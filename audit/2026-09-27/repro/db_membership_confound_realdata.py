"""Membership-only effect of the motif prior ON A REAL DATASET of EXP-20260925-LM-2.

Real code path: diagnose -> baseline plan (engine_bench.BASELINE) -> execute_plan ->
real_bench.score_models (pearson_signed, pearson_motif, motif_prior) with the evaluation half
from screen.split_reference(pairs, 0), exactly as real_bench does for the baseline arm.

The cisTarget DB is FAKE with RANDOM ranks (no sequence information). Two DBs per draw, built
from the same rank matrices:
  actual : absent genes are not columns (ENSG-named genes, optionally also non-protein-coding
           symbol genes)
  ctrl   : every universe gene is a column (the absent genes take the slots that held filler)
actual - ctrl = lift from DB membership alone on real expression and the real reference.
"""
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.feather as pf

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "benchmarks"))
from engine_bench import BASELINE  # noqa: E402
from real_bench import _reference, score_models  # noqa: E402
from grnbench.autoprep import diagnose, execute_plan, plan_preprocessing  # noqa: E402
from grnbench.profile import build_profile  # noqa: E402
from grnbench.screen import split_reference  # noqa: E402
from grnbench.selector import select  # noqa: E402
from grnbench.infer import cistarget as ct  # noqa: E402
import anndata as ad  # noqa: E402

HERE = Path(__file__).resolve().parent
N_DB = 5000
METHODS = ["pearson_signed", "pearson_motif", "motif_prior"]


def build_db(root: Path, cols: list[str], R: list[np.ndarray], motifs: list[str], ann: list[tuple]):
    root.mkdir(parents=True, exist_ok=True)
    for d, Rd in enumerate(R):
        arrays = [pa.array(motifs)] + [pa.array(Rd[:, j]) for j in range(Rd.shape[1])]
        tab = pa.Table.from_arrays(arrays, names=["motifs"] + cols)
        name = ["hg38_10kbp_up_10kbp_down", "hg38_500bp_up_100bp_down"][d]
        pf.write_feather(tab, str(root / f"{name}_fake.genes_vs_motifs.rankings.feather"))
    pd.DataFrame(ann, columns=["#motif_id", "gene_name"]).to_csv(
        root / "motifs-fake.hgnc-m0.001-o0.0.tbl", sep="\t", index=False)


def main(ds: str, n_draws: int, motif_mode: str, absent_mode: str):
    path = HERE / "data" / f"{ds}.h5ad"
    diag = diagnose(str(path))
    pairs = _reference(HERE / "refs", "human", diag.adata.var_names)
    tfs = sorted({s for s, _ in pairs})
    _sel, eval_half = split_reference(pairs, 0)
    prof = build_profile(diag.adata, diag.qc.metrics, diag.context, id_verdict=diag.id_verdict, tf_list=tfs)
    selection = select(prof)
    plan = plan_preprocessing(diag, tf_list=tfs, selection=selection,
                              choices={**selection.chosen(), **BASELINE})
    recs = execute_plan(plan, str(path), HERE / "work_rm" / f"{ds[:8]}_{motif_mode}_{absent_mode}", tf_list=tfs, seed=0, force=True)
    if not any(r["processed"] for r in recs):
        raise SystemExit(f"no processed context: {[r['error'] for r in recs]}")
    proc = ad.read_h5ad([r for r in recs if r["processed"]][0]["processed"], backed="r")
    genes = [str(g) for g in proc.var_names]
    gencode = pd.read_csv(HERE / "cxs/cellxgene_schema/gencode_files/genes_homo_sapiens.csv.gz",
                          header=None, names=["id", "name", "ver", "len", "type"])
    id2type = dict(zip(gencode.id, gencode.type))
    ids = proc.var["gene_id_as_loaded"].astype(str).tolist() if "gene_id_as_loaded" in proc.var else genes
    btype = {g: id2type.get(i.split(".")[0], "?") for g, i in zip(genes, ids)}
    ensg = {g for g in genes if g.startswith("ENSG")}
    tfset = {t.upper() for t in tfs}
    if absent_mode == "ensg":
        absent = sorted(ensg)
    elif absent_mode == "ensg+nonpc":
        absent = sorted(g for g in genes if (g in ensg or btype[g] != "protein_coding") and g.upper() not in tfset)
    else:
        raise SystemExit(absent_mode)
    present = [g for g in genes if g not in set(absent)]
    n_fill = N_DB - len(genes)
    fill_abs = [f"FILLA{i}" for i in range(len(absent))]
    fill = [f"FILL{i}" for i in range(n_fill)]
    cols_actual = present + fill_abs + fill
    cols_ctrl = present + absent + fill
    print(f"{ds[:8]} genes={len(genes)} tfs={len(tfs)} absent={len(absent)} ({len(absent)/len(genes):.3f}) "
          f"eval_pairs={len(eval_half)} motif_mode={motif_mode} absent_mode={absent_mode}", flush=True)
    rows = []
    for draw in range(n_draws):
        rng = np.random.default_rng(1000 + draw)
        if motif_mode == "one":
            k = np.ones(len(tfs), dtype=int)
        else:  # ~10% motif-less, the rest geometric with mean ~15, capped at 60
            k = np.where(rng.random(len(tfs)) < 0.1, 0, np.minimum(rng.geometric(1 / 15, len(tfs)), 60))
        motifs, ann = [], []
        for t, kk in zip(tfs, k):
            for j in range(kk):
                motifs.append(f"{t}__m{j}")
                ann.append((f"{t}__m{j}", t))
        R = [np.stack([rng.permutation(N_DB).astype(np.int16) for _ in motifs]) for _ in range(2)]
        res = {}
        for tag, cols in (("actual", cols_actual), ("ctrl", cols_ctrl)):
            root = Path(tempfile.mkdtemp(dir=HERE / "tmpd"))
            build_db(root, cols, R, motifs, ann)
            os.environ["GRNBENCH_CISTARGET_DIR"] = str(root)
            ct.load_annotation.cache_clear()
            out = score_models({"baseline": recs}, {"baseline": plan}, {"baseline": "baseline"},
                               METHODS, tfs, eval_half)
            for m in METHODS:
                res[f"{m}_{tag}"] = out[f"m_{m}_baseline_auprc_ratio"]
                if f"m_{m}_baseline_error" in out:
                    print("ERROR", m, tag, out[f"m_{m}_baseline_error"], flush=True)
        rows.append(res)
        print(f"  draw {draw}: " + " ".join(f"{k}={v:.4f}" for k, v in res.items()), flush=True)
    df = pd.DataFrame(rows)
    mem = df.pearson_motif_actual - df.pearson_motif_ctrl
    print(f"SUMMARY {ds[:8]} {motif_mode} {absent_mode} absent_share={len(absent)/len(genes):.3f} "
          f"pearson={df.pearson_signed_actual.mean():.4f} "
          f"pearson_motif(random) actual={df.pearson_motif_actual.mean():.4f} ctrl={df.pearson_motif_ctrl.mean():.4f} "
          f"membership_lift={mem.mean():+.4f} [{mem.min():+.4f},{mem.max():+.4f}] "
          f"motif_prior(random) actual={df.motif_prior_actual.mean():.4f} ctrl={df.motif_prior_ctrl.mean():.4f}",
          flush=True)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4])
