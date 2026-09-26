"""End-to-end preprocessing pipeline.

    load (counts-verified) -> QC gates -> filters -> strategy -> features -> outputs

Every run emits ``processed.h5ad``, ``qc_metrics.json``, ``manifest.json``
(full provenance: steps, params, environment, input fingerprint) and a
human-readable ``report.md``.
"""
from __future__ import annotations

import dataclasses
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import anndata as ad
import scipy.sparse as sp

from grnbench.features import select_features
from grnbench.ingest import load_any
from grnbench.qc import QCReport, Thresholds, run_qc
from grnbench.recommend import recommend_strategy
from grnbench.report import (build_manifest, file_fingerprint, render_report, write_json,
                             sequence_identity, require_empty_output)
from grnbench.strategies import REGISTRY, apply_strategy, get_strategy
from grnbench.tfs import match_tfs_in_data, resolve_tf_names


class PrepError(RuntimeError):
    """QC gate failure or unrecoverable preprocessing condition."""


@dataclass
class PrepResult:
    out_dir: Path
    processed_path: Path
    manifest_path: Path
    report_path: Path
    qc_metrics_path: Path
    strategy: str
    manifest: dict


def _nnz_per_gene(X) -> np.ndarray:
    return X.getnnz(axis=0) if sp.issparse(X) else (np.asarray(X) > 0).sum(axis=0)


def _stratified_indices(labels: np.ndarray, n_keep: int,
                        rng: np.random.Generator) -> tuple[np.ndarray, dict]:
    """Largest-remainder proportional allocation, >=1 per class when possible.

    Returns ``(indices, info)``; ``info`` records what the allocation could not
    honour, because a stratified subsample that quietly drops strata is worse
    than one that refuses.

    The previous version bumped every zero-quota class to 1 and then sliced
    ``picked[:n_keep]``. Since ``picked`` is built in ``np.unique`` (alphabetical)
    class order, that truncation deleted whole classes from the tail: 200
    singleton classes plus one class holding half the cells kept **zero** cells
    of the majority class, while the manifest still read "stratified by 'ct'".
    """
    classes, counts = np.unique(labels, return_counts=True)
    n_keep = int(n_keep)

    # Largest-remainder allocation already sums to exactly n_keep and keeps each
    # class's share proportional. The ">=1 per class" nicety is applied only
    # while there is room for it: bumping unconditionally is what used to
    # overshoot and force the truncation that deleted whole classes.
    quota = counts * (n_keep / len(labels))
    take = np.floor(quota).astype(int)
    remainder = int(n_keep - take.sum())
    order = np.argsort(-(quota - take), kind="stable")
    for i in range(remainder):
        take[order[i % len(classes)]] += 1
    take = np.minimum(take, counts)

    # spend anything still unallocated on classes that would otherwise get none
    spare = n_keep - int(take.sum())
    if spare > 0:
        for i in np.argsort(-counts, kind="stable"):
            if spare == 0:
                break
            if take[i] == 0 and counts[i] > 0:
                take[i] = 1
                spare -= 1

    zero = np.flatnonzero(take == 0)
    info: dict = {
        "n_classes": int(len(classes)),
        "classes_dropped": int(len(zero)),
        "classes_dropped_examples": [str(c) for c in classes[zero][:5]],
        "min_per_class_honoured": bool(len(zero) == 0),
        "policy": "largest-remainder proportional allocation",
    }

    picked: list[int] = []
    for cls, k in zip(classes, take):
        if k <= 0:
            continue
        pool = np.flatnonzero(labels == cls)
        picked.extend(rng.choice(pool, size=int(min(k, len(pool))),
                                 replace=False).tolist())
    info["n_kept"] = len(picked)
    return np.sort(np.asarray(picked, dtype=int)), info


def run_prep(
    input_path: str | Path,
    out_dir: str | Path,
    *,
    strategy: str = "auto",
    tf_list: list[str] | None = None,
    n_hvg: int | None = None,
    group_key: str | None = None,
    label_key: str | None = None,
    label_mapping: dict[str, str] | None = None,
    stratify_key: str | None = None,
    batch_key: str | None = None,
    min_group_cells: int | None = None,
    downcast: bool = True,
    thresholds: Thresholds | None = None,
    subsample_cells: int | None = None,
    orientation: str = "auto",
    drop_doublets: bool = False,
    force: bool = False,
    seed: int = 42,
    subset_key: str | None = None,
    subset_value: str | None = None,
    cell_calling: str | None = None,
    exclude_qc_genes: bool = False,
    per_batch_bounds: dict | None = None,
    cells_per_metacell: int | None = None,
) -> PrepResult:
    """Preprocess a dataset for GRN inference. See module docstring."""
    t = thresholds or Thresholds()
    out = require_empty_output(out_dir)
    input_identity = file_fingerprint(input_path)
    execution_params = {
        'strategy': strategy, 'n_hvg': n_hvg, 'group_key': group_key,
        'label_key': label_key, 'label_mapping': label_mapping,
        'stratify_key': stratify_key, 'batch_key': batch_key,
        'min_group_cells': min_group_cells, 'downcast': downcast,
        'thresholds': dataclasses.asdict(t), 'subsample_cells': subsample_cells,
        'orientation': orientation, 'drop_doublets': drop_doublets,
        'force': force, 'seed': seed,
        'subset_key': subset_key, 'subset_value': subset_value, 'cell_calling': cell_calling,
        'exclude_qc_genes': exclude_qc_genes, 'per_batch_bounds': per_batch_bounds,
        'cells_per_metacell': cells_per_metacell,
        'requested_tfs': sequence_identity(tf_list) if tf_list is not None else None,
    }
    rng = np.random.default_rng(seed)
    steps: list[dict] = []
    warnings_list: list[str] = []

    # -- 1. ingest (counts-verified) ---------------------------------------- #
    adata, ingest_rep = load_any(input_path, orientation=orientation)
    steps.append({"step": "ingest",
                  "detail": f"format={ingest_rep.fmt}, counts_source={ingest_rep.counts_source}"})
    if subset_key is not None and subset_key not in adata.obs.columns:
        raise PrepError(f"subset key '{subset_key}' not in obs")
    if cell_calling is not None:
        if cell_calling != "knee":
            raise PrepError(f"unknown cell calling method {cell_calling!r} (available: knee)")
        from grnbench.hygiene import knee_cell_threshold
        totals = np.asarray(adata.X.sum(axis=1)).ravel()
        call = knee_cell_threshold(totals)
        if call is None:
            raise PrepError("knee cell calling found no clear cell/empty-droplet cliff in the "
                            "barcode-rank curve; is this an unfiltered droplet matrix?")
        adata = adata[totals >= call["threshold"]].copy()
        steps.append({"step": "cell_calling", "n_cells_after": int(adata.n_obs),
                      "n_genes_after": int(adata.n_vars), **call,
                      "detail": f"knee: kept {call['n_called']} of {call['n_barcodes']} barcodes with "
                                f">= {call['threshold']:.0f} counts"})

    # -- 1b. memory hygiene + metadata context ------------------------------- #
    from grnbench import hygiene as hyg
    if hyg.sparsify_if_dense(adata):
        steps.append({"step": "sparsify",
                      "detail": "dense matrix with >=80% zeros converted to CSR"})
    if downcast and hyg.downcast_floats(adata):
        steps.append({"step": "downcast", "detail": "float64 -> float32"})

    id_verdict = hyg.detect_id_type(adata.var_names)
    organism = hyg.organism_hint(adata.var_names)
    n_dup_barcodes = hyg.duplicate_barcodes_report(adata)
    if n_dup_barcodes:
        # Cells are never matched to anything by name, so a suffix costs nothing;
        # duplicated names, on the other hand, make every obs lookup ambiguous.
        adata.obs["barcode_as_loaded"] = adata.obs_names.astype(str).to_numpy()
        adata.obs_names_make_unique()
        steps.append({"step": "dedup_barcodes", "n_duplicates": int(n_dup_barcodes),
                      "detail": f"{n_dup_barcodes} duplicated obs_names made unique with "
                                "numeric suffixes; originals kept in obs['barcode_as_loaded']"})
    batch_col, batch_how = hyg.detect_batch_column(adata, explicit=batch_key)
    if batch_col is not None:
        steps.append({"step": "batch_detect",
                      "detail": f"column '{batch_col}' ({batch_how})"})
    if organism == "mouse-like":
        warnings_list.append(
            "gene symbols look mouse-like: make sure the TF list matches the "
            "organism (human TF lists will silently under-match)")

    label_report = None
    if label_key is not None:
        if label_key not in adata.obs.columns:
            raise PrepError(f"label key '{label_key}' not in obs")
        # label_mapping arrives already parsed (dict old->new)
        adata.obs[label_key], label_report = hyg.normalize_labels(
            adata.obs[label_key], mapping=label_mapping)
        steps.append({"step": "label_hygiene",
                      "detail": ", ".join(f"{k}={v}" for k, v in label_report.items()
                                          if k.startswith("n_"))})
        if label_report["n_missing"]:
            warnings_list.append(
                f"{label_report['n_missing']} cells have no usable '{label_key}' label "
                f"(tokens seen: {label_report['missing_tokens_seen']}); they are kept in "
                "the matrix and excluded from stratified subsampling - pass the column as "
                "--group-key instead if unlabelled cells should be dropped")

    if group_key is not None and group_key in adata.obs.columns:
        adata.obs[group_key], _grp_rep = hyg.normalize_labels(adata.obs[group_key])
        missing = adata.obs[group_key].isna()
        if missing.any():
            adata = adata[~missing].copy()
            steps.append({"step": "drop_unlabelled_groups",
                          "n_cells_after": int(adata.n_obs),
                          "n_genes_after": int(adata.n_vars),
                          "detail": f"removed {int(missing.sum())} cells with "
                                    f"missing '{group_key}'"})

    # -- 2. QC on the (hygiene-passed) matrix -------------------------------- #
    qc_report: QCReport = run_qc(adata, t, group_key=group_key, label_key=label_key,
                                 batch_key=batch_col, id_verdict=id_verdict,
                                 n_dup_barcodes=n_dup_barcodes)
    steps.append({"step": "qc", "detail": f"{len(qc_report.gates)} gates evaluated"})

    # -- 3. recommendations + strategy choice -------------------------------- #
    med_cells_group = n_groups = None
    if group_key is not None:
        if group_key not in adata.obs.columns:
            cols = ", ".join(map(str, adata.obs.columns[:25]))
            raise PrepError(f"group key '{group_key}' not in obs (columns: {cols})")
        counts = adata.obs[group_key].astype(str).value_counts()
        med_cells_group, n_groups = float(counts.median()), int(len(counts))
    recommendations = recommend_strategy(
        qc_report.metrics, group_key=group_key, n_groups=n_groups,
        median_cells_per_group=med_cells_group,
    )

    if strategy == "auto":
        chosen = "cpm_log1p"
        for rec in recommendations:
            spec = REGISTRY.get(rec["strategy"])
            if spec is None or spec.unavailable is not None:
                continue
            if rec["strategy"] == "pseudobulk" and group_key is None:
                continue  # not applicable without a group key
            chosen = rec["strategy"]
            break
    else:
        chosen = strategy
    get_strategy(chosen)  # validates availability early

    # -- 4. gate enforcement ------------------------------------------------- #
    if qc_report.has_failures and not force:
        failed = [g for g in qc_report.gates if g.status == "fail"]
        detail = "; ".join(f"{g.name}: {g.detail}" for g in failed)
        raise PrepError(f"QC gate failure(s) - {detail}")

    # -- 5. cell / gene filters (logged) -------------------------------------- #
    low_genes = adata.obs["_n_genes"].to_numpy() < t.min_genes_per_cell
    zero_total = adata.obs["_total_counts"].to_numpy() <= 0
    upper = {}
    for label, col, bound in (("> max genes", "_n_genes", t.max_genes_per_cell),
                              ("> max counts", "_total_counts", t.max_counts_per_cell),
                              ("> max % mito", "_pct_mito", t.max_pct_mito_cell)):
        if bound is not None:
            upper[f"{label} ({bound:g})"] = adata.obs[col].to_numpy() > bound
    if per_batch_bounds:
        if batch_col is None:
            raise PrepError("per_batch_bounds needs a batch column (--batch-key)")
        labels = adata.obs[batch_col].astype(str).to_numpy()
        unknown = sorted(set(map(str, per_batch_bounds)) - set(labels))
        if unknown:
            raise PrepError(f"per_batch_bounds names batches not in '{batch_col}': {unknown[:5]}")
        hit = np.zeros(adata.n_obs, dtype=bool)
        for b, bounds in per_batch_bounds.items():
            rows = labels == str(b)
            for key, col, below in (("min_genes_per_cell", "_n_genes", True),
                                    ("max_genes_per_cell", "_n_genes", False),
                                    ("max_counts_per_cell", "_total_counts", False),
                                    ("max_pct_mito_cell", "_pct_mito", False)):
                if bounds.get(key) is None:
                    continue
                v = adata.obs[col].to_numpy()
                hit |= rows & ((v < bounds[key]) if below else (v > bounds[key]))
        upper[f"outside per-batch bounds ('{batch_col}')"] = hit
    mask_cells = ~(low_genes | zero_total)
    for hit in upper.values():
        mask_cells &= ~hit
    if not mask_cells.all():
        adata = adata[mask_cells].copy()
    # A zero-count cell is necessarily also a low-gene cell, so reporting both
    # counts side by side double-counted every one of them and the numbers did
    # not add up to the cells actually removed. Report the union, then the
    # zero-count cells as the subset they are.
    n_removed = int((~mask_cells).sum())
    steps.append({"step": "filter_cells",
                  "n_cells_after": int(adata.n_obs), "n_genes_after": int(adata.n_vars),
                  "n_cells_removed": n_removed,
                  "detail": f"removed {n_removed} cells with "
                            f"< {t.min_genes_per_cell} detected genes "
                            f"(of which {int(zero_total.sum())} had zero counts)"
                            + "".join(f"; {int(hit.sum())} flagged {label}"
                                      for label, hit in upper.items())})

    mask_genes = _nnz_per_gene(adata.X) >= t.min_cells_per_gene
    if not mask_genes.all():
        adata = adata[:, mask_genes].copy()
    steps.append({"step": "filter_genes",
                  "n_cells_after": int(adata.n_obs), "n_genes_after": int(adata.n_vars),
                  "detail": f"removed {int((~mask_genes).sum())} genes "
                            f"detected in < {t.min_cells_per_gene} cells"})
    if exclude_qc_genes:
        from grnbench.qc import MITO_PREFIXES, RIBO_PREFIXES, pd_index_upper
        upper_names = pd_index_upper(adata.var_names)
        qc_genes = np.zeros(adata.n_vars, dtype=bool)
        for pref in tuple(MITO_PREFIXES) + tuple(RIBO_PREFIXES):
            qc_genes |= upper_names.str.startswith(pref)
        if qc_genes.any():
            adata = adata[:, ~qc_genes].copy()
        steps.append({"step": "exclude_qc_genes",
                      "n_cells_after": int(adata.n_obs), "n_genes_after": int(adata.n_vars),
                      "detail": f"removed {int(qc_genes.sum())} mitochondrial/ribosomal genes "
                                "from the features (they track cell quality, not regulation)"})

    if drop_doublets:
        method, is_dbl, dbl_warning = hyg.detect_doublets(adata, batch_key=batch_col, seed=seed)
        if dbl_warning:
            warnings_list.append(dbl_warning)
        n_dbl = int(is_dbl.sum())
        rate = n_dbl / max(adata.n_obs, 1)
        if rate > t.doublet_rate_max:
            # On data without cluster structure scrublet's simulated doublets look
            # like everything else and the auto-threshold lands inside the singlet
            # mode: 83% of a homogeneous test set came back "doublet". That is a
            # failed threshold, not a doublet population; dropping it is data loss.
            raise PrepError(
                f"{method} called {n_dbl}/{adata.n_obs} cells ({rate:.0%}) doublets, above "
                f"the {t.doublet_rate_max:.0%} ceiling. Real doublet rates are ~1-10%; a call "
                "rate like this means the threshold failed (no cluster structure, or batches "
                "too small), not that most cells are doublets. Nothing was dropped. Inspect "
                "the scores without --drop-doublets, or raise --max-doublet-rate deliberately.")
        if n_dbl:
            adata = adata[~is_dbl].copy()
        steps.append({"step": "drop_doublets", "method": method,
                      "fallback": method != "scrublet",
                      "fallback_reason": dbl_warning if method != "scrublet" else None,
                      "n_removed": n_dbl,
                      "n_cells_after": int(adata.n_obs), "n_genes_after": int(adata.n_vars),
                      "detail": f"removed {n_dbl} predicted doublets ({method}"
                                + (f", per batch '{batch_col}'" if batch_col and method == "scrublet" else "")
                                + "); scores kept in obs['doublet_score']"})

    # Context subsetting comes after cell QC and doublet calling on ALL cells:
    # inside one cell type there is no cluster structure, and doublet callers
    # over-call (scrublet flagged 22-28% of a single demo cell type). Genes are
    # then re-filtered within the context.
    if subset_key is not None:
        keep = adata.obs[subset_key].astype(str).to_numpy() == str(subset_value)
        if not keep.any():
            raise PrepError(f"no cells with {subset_key} == {subset_value!r} after QC")
        adata = adata[keep].copy()
        steps.append({"step": "subset_context", "n_cells_after": int(adata.n_obs),
                      "n_genes_after": int(adata.n_vars),
                      "detail": f"kept {int(keep.sum())} cells with {subset_key} == {subset_value!r}"})
        mask_ctx = _nnz_per_gene(adata.X) >= t.min_cells_per_gene
        if not mask_ctx.all():
            adata = adata[:, mask_ctx].copy()
        steps.append({"step": "filter_genes_in_context",
                      "n_cells_after": int(adata.n_obs), "n_genes_after": int(adata.n_vars),
                      "detail": f"removed {int((~mask_ctx).sum())} genes detected in < "
                                f"{t.min_cells_per_gene} cells of this context"})

    if subsample_cells is not None and int(subsample_cells) < adata.n_obs:
        strat_col = stratify_key or group_key or label_key
        n_unlabelled = 0
        if strat_col is not None and strat_col in adata.obs.columns:
            # Missing labels are not a stratum. Folding them into a literal "nan"
            # class gave them a proportional quota like any real class.
            unlabelled = hyg.missing_label_mask(adata.obs[strat_col])
            n_unlabelled = int(unlabelled.sum())
            pool = np.flatnonzero(~unlabelled)
            if n_unlabelled:
                warnings_list.append(
                    f"{n_unlabelled} cells with missing '{strat_col}' were excluded from "
                    f"the stratified subsample")
            rel, strat_info = _stratified_indices(
                adata.obs[strat_col].iloc[pool].astype(str).to_numpy(),
                min(int(subsample_cells), len(pool)), rng)
            idx = pool[rel]
            kind = f"stratified by '{strat_col}'"
            if strat_info["classes_dropped"]:
                ex = ", ".join(strat_info["classes_dropped_examples"])
                warn = (f"subsample_cells={int(subsample_cells)} is smaller than the "
                        f"{strat_info['n_classes']} classes in '{strat_col}': "
                        f"{strat_info['classes_dropped']} class(es) are absent from "
                        f"the subsample (e.g. {ex}). The largest classes were kept.")
                warnings_list.append(warn)
                kind += f", {strat_info['classes_dropped']} class(es) dropped"
        else:
            strat_info = None
            idx = np.sort(rng.choice(adata.n_obs, size=int(subsample_cells),
                                     replace=False))
            kind = "random (no stratify key)"
        adata = adata[idx].copy()
        steps.append({"step": "subsample_cells",
                      "n_cells_after": int(adata.n_obs), "n_genes_after": int(adata.n_vars),
                      "stratification": strat_info, "n_unlabelled_excluded": n_unlabelled,
                      "detail": f"{kind}, to {int(subsample_cells)} cells (seed={seed})"})

    # -- 6. normalization strategy ------------------------------------------- #
    params: dict = {}
    if chosen == "pseudobulk":
        if group_key is None:
            raise PrepError("strategy 'pseudobulk' requires --group-key")
        params["group_key"] = group_key
        if min_group_cells is not None:
            params["min_group_cells"] = int(min_group_cells)
    if chosen == "cpm_log1p_center_batch":
        if batch_col is None:
            raise PrepError("strategy 'cpm_log1p_center_batch' needs a batch column")
        params["batch_key"] = batch_col
    if chosen == "metacells":
        params.update(batch_key=batch_col, label_key=stratify_key or label_key, seed=seed)
        if cells_per_metacell is not None:
            params["cells_per_metacell"] = int(cells_per_metacell)
    if chosen == "cpm_log1p_combat":
        if batch_col is None:
            raise PrepError("strategy 'cpm_log1p_combat' needs a batch column: pass "
                            "--batch-key, or have a batch/donor/sample-like obs column")
        params["batch_key"] = batch_col
    # A strategy warns via warnings.warn (e.g. pseudobulk dropping tiny groups);
    # that is invisible to anyone who only reads manifest.json/report.md, so it
    # is captured here rather than left to whoever happens to see stderr.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        processed, prov = apply_strategy(chosen, adata, **params)
    warnings_list.extend(str(w.message) for w in caught)
    if processed.n_obs == 0 or processed.n_vars == 0:
        raise PrepError(
            f"strategy '{chosen}' produced an empty matrix ({processed.n_obs} cells x "
            f"{processed.n_vars} genes) - check --min-group-cells and group sizes; "
            "an empty processed.h5ad cannot feed inference"
        )
    steps.append({"step": "strategy", **prov,
                  "n_cells_after": int(processed.n_obs), "n_genes_after": int(processed.n_vars)})

    # -- 7. optional feature selection (HVG union TFs) ------------------------ #
    feature_report = None
    if n_hvg is not None:
        processed, feature_report = select_features(processed, n_hvg=int(n_hvg),
                                                    tf_list=tf_list)
        steps.append({"step": "feature_selection",
                      "n_cells_after": int(processed.n_obs),
                      "n_genes_after": int(processed.n_vars),
                      "detail": ", ".join(f"{k}={v}" for k, v in feature_report.items())})
    tf_coverage = None
    if tf_list:
        matched, requested = match_tfs_in_data(tf_list, processed.var_names)
        tf_coverage = {"requested": requested, "matched_in_final_matrix": len(matched)}
        rate = (len(matched) / requested) if requested else 0.0
        if rate < 0.3:
            warnings_list.append(
                f"TF coverage low: only {len(matched)}/{requested} "
                f"({rate:.0%}) TFs found in the final matrix - check organism "
                "and gene-ID type before trusting TF-restricted inference")

    # -- 8. outputs ----------------------------------------------------------- #
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    processed_path = out / "processed.h5ad"
    qc_path = out / "qc_metrics.json"
    manifest_path = out / "manifest.json"
    report_path = out / "report.md"

    processed.write_h5ad(processed_path, compression="gzip")
    write_json(qc_path, {"qc": qc_report.to_dict(), "ingest": ingest_rep.to_dict(),
                         "thresholds": dataclasses.asdict(t)})

    outputs = {"processed_h5ad": str(processed_path), "qc_metrics": str(qc_path),
               "manifest": str(manifest_path), "report": str(report_path)}
    manifest = build_manifest(
        input_fp=input_identity,
        ingest=ingest_rep.to_dict(),
        qc=qc_report.to_dict(),
        steps=steps,
        strategy=chosen,
        recommendations=recommendations,
        outputs=outputs,
        extra={
            'args': execution_params,
            'effective_genes': sequence_identity(processed.var_names),
            'effective_tfs': sequence_identity(resolve_tf_names(processed.var_names, tf_list)),
            "seed": seed,
            "forced_past_gates": bool(force and qc_report.has_failures),
            "thresholds": dataclasses.asdict(t),
            "feature_selection": feature_report,
            "tf_coverage": tf_coverage,
            "warnings": warnings_list,
            "label_audit": label_report,
            "gene_id_type": id_verdict,
            "organism_hint": organism,
        },
    )
    write_json(manifest_path, manifest)
    report_path.write_text(render_report(manifest), encoding="utf-8")

    return PrepResult(out_dir=out, processed_path=processed_path,
                      manifest_path=manifest_path, report_path=report_path,
                      qc_metrics_path=qc_path, strategy=chosen, manifest=manifest)
