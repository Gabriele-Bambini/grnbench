"""Scenario simulator for testing the cleaning engine against a known network.

Each :class:`Scenario` generates raw counts from a planted TF -> target
network (shared by every cell type) and injects the problems real data have,
one knob each: low depth, batches (balanced or confounded with cell type),
damaged cells, near-empty droplets, doublets, ambient RNA, full-length read
counts, Ensembl or mouse-style identifiers, an unfiltered droplet matrix, and a
README that is correct, missing or misleading.

``obs['simulated_state']`` records what every barcode really is (``ok``,
``damaged``, ``empty``, ``doublet``, ``raw_empty``) so cleaning accuracy can be
measured, and :meth:`Scenario.expected` states what a correct engine should
detect and choose.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

DROPLET_README = ("Single-cell RNA-seq of {tissue} ({organism}), 10x Genomics Chromium v3, raw UMI "
                  "counts, {n_cells} cells. Batches in column `batch`, cell types in `cell_type`.")
PLATE_README = ("Single-cell RNA-seq of {tissue} ({organism}) with Smart-seq2 full-length protocol, "
                "read counts (not normalized), {n_cells} cells.")


@dataclass(frozen=True)
class Scenario:
    name: str = "clean"
    n_cells: int = 1500
    n_genes: int = 1000
    n_tfs: int = 30
    targets_per_tf: int = 12
    n_cell_types: int = 3
    type_program_sd: float = 0.8
    effect: float = 0.7
    n_batches: int = 1
    batch_shift_sd: float = 0.0
    batch_depth_ratio: float = 1.0
    batch_confounding: float = 0.0
    median_library: float = 3000.0
    library_sd: float = 0.35
    damaged_frac: float = 0.0
    empty_frac: float = 0.0
    doublet_frac: float = 0.0
    ambient_frac: float = 0.0
    protocol: str = "umi"            # umi | full_length
    gene_ids: str = "symbols"        # symbols | ensembl | mouse
    raw_droplets: int = 0            # extra near-empty barcodes (unfiltered matrix)
    readme: str = "correct"          # correct | none | misleading
    seed: int = 0

    def with_seed(self, seed: int) -> "Scenario":
        return replace(self, seed=seed)

    def expected(self) -> dict:
        """What a correct engine should detect / choose on this scenario."""
        exp = {"technology": self.protocol,
               "organism": "mouse" if self.gene_ids == "mouse" else "human",
               "batch_detected": self.n_batches > 1}
        if self.protocol == "full_length":
            exp["doublets"] = "none"
            exp["normalization"] = "cpm_log1p"
        elif self.n_cells >= 1000:
            exp["doublets_removed"] = True
        if self.n_batches > 1 and self.batch_confounding >= 0.8:
            exp["batch"] = "none"
        if self.gene_ids == "ensembl":
            exp["identifier_mapping_best"] = "ensembl_to_symbol"
        if self.raw_droplets:
            exp["cell_calling_best"] = {"emptydrops", "knee"}
        types_big = self.n_cell_types >= 2 and self.n_cells / self.n_cell_types >= 250
        if types_big:
            exp["context_scope"] = "per_cell_type"
        if self.readme == "misleading":
            exp["conflict"] = "technology"
        return exp


def _gene_names(scn: Scenario, n_mito: int):
    tfs = [f"TF{i:02d}" for i in range(1, scn.n_tfs + 1)]
    others = [f"GENE{i:04d}" for i in range(1, scn.n_genes - scn.n_tfs - n_mito + 1)]
    mito = [f"MT-GENE{i:02d}" for i in range(1, n_mito + 1)]
    names = tfs + others + mito
    if scn.gene_ids == "mouse":
        names = [n[:1] + n[1:].lower() if not n.startswith("MT-") else "mt-" + n[3:].capitalize()
                 for n in names]
    elif scn.gene_ids == "ensembl":
        names = [f"ENSG{i + 1:011d}" for i in range(len(names))]
    return names, names[: scn.n_tfs], names[scn.n_tfs: scn.n_genes - n_mito]


def simulate(scn: Scenario):
    """Return ``(adata, tf_names, gold)``; gold columns: source, target, sign."""
    rng = np.random.default_rng(scn.seed)
    # ~1.3% of genes are mitochondrial, as in a 1,000-gene panel with 13 MT genes; a fixed
    # 13 on a 300-gene panel pushed the median mito fraction above the 25% QC gate.
    n_mito = max(1, round(13 * scn.n_genes / 1000)) if scn.gene_ids != "ensembl" else 0
    names, tfs, others = _gene_names(scn, n_mito)
    n_c, n_g, n_t = scn.n_cells, len(names), scn.n_tfs

    pool = others[: min(len(others), n_t * scn.targets_per_tf * 2)]
    edges = []
    for tf in tfs:
        for tgt in rng.choice(pool, scn.targets_per_tf, replace=False):
            edges.append((tf, str(tgt), int(rng.choice([-1, 1], p=[0.3, 0.7]))))
    gold = pd.DataFrame(edges, columns=["source", "target", "sign"]).drop_duplicates(["source", "target"])
    index = {g: i for i, g in enumerate(names)}

    types = [f"type{k + 1}" for k in range(scn.n_cell_types)]
    weights = np.linspace(1.0, 0.5, scn.n_cell_types)
    cell_type = rng.choice(types, n_c, p=weights / weights.sum())
    batches = [f"run{b + 1}" for b in range(scn.n_batches)]
    if scn.n_batches > 1 and scn.batch_confounding > 0:
        home = {t: batches[i % scn.n_batches] for i, t in enumerate(types)}
        batch = np.array([home[t] if rng.random() < scn.batch_confounding else rng.choice(batches)
                          for t in cell_type])
    else:
        batch = rng.choice(batches, n_c)

    program = {t: rng.normal(0, scn.type_program_sd, n_t) for t in types}
    tf_log = np.stack([program[t] for t in cell_type]) + rng.normal(0, 0.6, (n_c, n_t))
    z = (tf_log - tf_log.mean(0)) / tf_log.std(0)
    base = rng.normal(-1.0, 1.2, n_g)
    base[:n_t] = rng.normal(0.8, 0.3, n_t)
    if n_mito:
        base[-n_mito:] = rng.normal(2.5, 0.3, n_mito)
    log_mu = np.tile(base, (n_c, 1))
    log_mu[:, :n_t] += tf_log - tf_log.mean(0)
    for s, t, sign in gold.itertuples(index=False):
        log_mu[:, index[t]] += sign * scn.effect * z[:, tfs.index(s)]
    lib = np.exp(rng.normal(np.log(scn.median_library), scn.library_sd, n_c))
    for b_i, b in enumerate(batches[1:], start=1):
        rows = batch == b
        log_mu[rows] += rng.normal(0, scn.batch_shift_sd, n_g) if scn.batch_shift_sd else 0
        lib[rows] *= scn.batch_depth_ratio ** (b_i / max(1, scn.n_batches - 1))

    state = rng.choice(["ok", "damaged", "empty", "doublet"], n_c,
                       p=_state_p(scn.damaged_frac, scn.empty_frac, scn.doublet_frac))
    lib[state == "empty"] *= 0.03
    if n_mito:
        log_mu[state == "damaged", -n_mito:] += np.log(8.0)
    lib[state == "damaged"] *= 0.5

    mu = np.exp(log_mu)
    mu = mu / mu.sum(1, keepdims=True)
    soup = mu.mean(0)
    if scn.ambient_frac:
        mu = (1 - scn.ambient_frac) * mu + scn.ambient_frac * soup
    X = rng.poisson(rng.gamma(5.0, mu * lib[:, None] / 5.0)).astype(np.float32)
    dbl = np.flatnonzero(state == "doublet")
    X[dbl] = X[dbl] + X[rng.integers(0, n_c, len(dbl))]
    if scn.protocol == "full_length":
        amp = np.exp(rng.normal(np.log(40), 0.6, n_g))          # length x amplification per gene
        X = rng.poisson(X * amp).astype(np.float32)

    obs = pd.DataFrame({"cell_type": cell_type, "simulated_state": state},
                       index=[f"cell{i:05d}" for i in range(n_c)])
    if scn.n_batches > 1:
        obs.insert(0, "batch", batch)
    if scn.raw_droplets:
        k = scn.raw_droplets
        lib_e = np.exp(rng.normal(np.log(20), 0.8, k))
        E = rng.poisson(soup[None, :] * lib_e[:, None]).astype(np.float32)
        X = np.vstack([X, E])
        extra = pd.DataFrame({"cell_type": "none", "simulated_state": "raw_empty"},
                             index=[f"bc{i:06d}" for i in range(k)])
        if "batch" in obs:
            extra.insert(0, "batch", rng.choice(batches, k))
        obs = pd.concat([obs, extra])
    adata = ad.AnnData(X=sparse.csr_matrix(X), obs=obs, var=pd.DataFrame(index=names))
    return adata, tfs, gold


def _state_p(damaged, empty, doublet):
    ok = 1.0 - damaged - empty - doublet
    if ok <= 0:
        raise ValueError("damaged + empty + doublet fractions must be < 1")
    return [ok, damaged, empty, doublet]


def readme_text(scn: Scenario, n_cells: int) -> str | None:
    if scn.readme == "none":
        return None
    organism = "Mus musculus" if scn.gene_ids == "mouse" else "Homo sapiens"
    plate = scn.protocol == "full_length"
    if scn.readme == "misleading":
        plate = not plate
    tpl = PLATE_README if plate else DROPLET_README
    return tpl.format(tissue="bone marrow", organism=organism, n_cells=n_cells)


def write_scenario(out_dir, scn: Scenario) -> dict[str, Path]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata, tfs, gold = simulate(scn)
    paths = {"counts": out / "counts.h5ad", "tfs": out / "tfs.txt", "gold": out / "gold.tsv",
             "scenario": out / "scenario.json"}
    adata.write_h5ad(paths["counts"])
    paths["tfs"].write_text("\n".join(tfs) + "\n", encoding="utf-8")
    gold.to_csv(paths["gold"], sep="\t", index=False)
    paths["scenario"].write_text(json.dumps(asdict(scn), indent=1), encoding="utf-8")
    text = readme_text(scn, scn.n_cells)
    if text:
        paths["readme"] = out / "README.md"
        paths["readme"].write_text(text + "\n", encoding="utf-8")
    return paths


PRESETS: dict[str, Scenario] = {s.name: s for s in [
    Scenario("clean"),
    Scenario("low_depth", median_library=800),
    Scenario("very_sparse", median_library=400, n_cells=3000),
    Scenario("batch_balanced", n_batches=2, batch_shift_sd=0.4, batch_depth_ratio=1.8),
    Scenario("batch_confounded", n_batches=3, batch_shift_sd=0.4, batch_confounding=1.0),
    Scenario("damaged", damaged_frac=0.10),
    Scenario("doublets", doublet_frac=0.08),
    Scenario("ambient", ambient_frac=0.15),
    Scenario("empty_droplets", empty_frac=0.08),
    Scenario("strong_composition", n_cell_types=4, type_program_sd=1.5),
    Scenario("single_type", n_cell_types=1),
    Scenario("full_length", protocol="full_length", n_cells=800),
    Scenario("mixed_problems", n_batches=2, batch_shift_sd=0.3, damaged_frac=0.05,
             doublet_frac=0.05, ambient_frac=0.08, empty_frac=0.03, median_library=1500),
    Scenario("mouse_ids", gene_ids="mouse"),
    Scenario("ensembl_ids", gene_ids="ensembl"),
    Scenario("misleading_readme", readme="misleading"),
    Scenario("no_readme", readme="none"),
    Scenario("small", n_cells=300, n_genes=500, n_tfs=15, targets_per_tf=8),
    Scenario("raw_droplets", raw_droplets=22000, n_cells=1200),
]}


def random_scenario(seed: int) -> Scenario:
    """A random mixture of every knob (for massive testing beyond the presets)."""
    r = np.random.default_rng(10_000 + seed)
    n_b = int(r.choice([1, 1, 2, 3]))
    return Scenario(
        name=f"random_{seed}", seed=seed,
        n_cells=int(r.choice([400, 800, 1500, 3000])), n_genes=int(r.choice([600, 1000, 1500])),
        n_tfs=int(r.choice([15, 30, 45])), targets_per_tf=int(r.choice([8, 12, 16])),
        n_cell_types=int(r.choice([1, 2, 3, 4])), type_program_sd=float(r.uniform(0.2, 1.6)),
        effect=float(r.uniform(0.3, 1.0)), n_batches=n_b,
        batch_shift_sd=float(r.uniform(0, 0.6)) if n_b > 1 else 0.0,
        batch_depth_ratio=float(r.uniform(1, 2.5)) if n_b > 1 else 1.0,
        batch_confounding=float(r.choice([0, 0, 0.5, 1.0])) if n_b > 1 else 0.0,
        median_library=float(r.choice([500, 1000, 2500, 5000])),
        damaged_frac=float(r.choice([0, 0.03, 0.08])), empty_frac=float(r.choice([0, 0.03, 0.06])),
        doublet_frac=float(r.choice([0, 0.03, 0.08])), ambient_frac=float(r.choice([0, 0.05, 0.15])),
        protocol=str(r.choice(["umi", "umi", "umi", "full_length"])),
        readme=str(r.choice(["correct", "correct", "none", "misleading"])),
    )
