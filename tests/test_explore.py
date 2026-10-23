"""grnbench.explore: atlas, TF filter, detection, card and CLI on a tiny synthetic artifacts folder."""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from grnbench import explore

GENES = ["TFA", "TFB"] + [f"G{i}" for i in range(1, 13)]


def _dataset(root: Path, name: str, ref: list[tuple[str, str]], seed: int) -> None:
    rng = np.random.default_rng(seed)
    X = rng.poisson(2.0, (50, len(GENES))).astype(np.float32)
    X[:, 0] = 0
    X[:20, 0] = 1.0                      # TFA detected in exactly 20/50 = 0.4 of cells
    ctx = root / name / "auto" / "ctx0"
    ctx.mkdir(parents=True)
    ad.AnnData(X=X, var=pd.DataFrame(index=GENES)).write_h5ad(ctx / "processed.h5ad")
    pd.DataFrame(ref, columns=["source", "target"]).to_csv(root / name / "reference_evaluation.tsv", sep="\t", index=False)
    for m in explore.MODELS:
        rows = [(tf, g, (10.0 if (tf, g) == ("TFA", "G1") else float(rng.normal())))
                for tf in ("TFA", "TFB") for g in GENES[2:]]
        pd.DataFrame(rows, columns=["source", "target", "weight"]).to_csv(ctx / f"{m}.tsv.gz", sep="\t", index=False)


def _artifacts(tmp_path: Path) -> Path:
    root = tmp_path / "art"
    ref1 = [("TFA", f"G{i}") for i in range(1, 7)] + [("TFB", f"G{i}") for i in range(1, 4)]
    ref2 = [("TFA", f"G{i}") for i in range(1, 6)] + [("TFA", "NOTINDATA")] + [("TFB", f"G{i}") for i in range(1, 6)]
    _dataset(root, "cxg__a__homo-sapiens__t1__aaaa", ref1, 1)
    _dataset(root, "cxg__a__homo-sapiens__t2__bbbb", ref2, 2)
    _dataset(root, "cxg__a__homo-sapiens__t3__cccc", ref1, 3)
    return root


def test_atlas_filter_and_detection(tmp_path):
    atlas = explore.build_atlas(explore.load_all(_artifacts(tmp_path)))
    got = {(r.dataset[-4:], r.tf): r for r in atlas.itertuples()}
    assert set(got) == {("aaaa", "TFA"), ("bbbb", "TFA"), ("bbbb", "TFB"), ("cccc", "TFA")}   # TFB has 3 targets in 1 and 3
    assert got[("bbbb", "TFA")].n_ref_targets == 5                            # NOTINDATA is outside the universe
    assert abs(got[("aaaa", "TFA")].detection - 0.4) < 1e-9
    assert got[("aaaa", "TFA")].auprc_ratio_pearson_motif > 1.0               # G1 ranked first and is a reference target
    assert (atlas.status == "EXPLORATORY").all()
    assert (atlas.n_stored_pearson_signed == 12).all()                        # stored edges per TF are recorded


def test_ranked_in_both_excludes_unranked(tmp_path):
    atlas = explore.build_atlas(explore.load_all(_artifacts(tmp_path)))
    atlas.loc[0, "n_stored_pearson_motif"] = 0                                # TF unranked with the prior
    s = explore.detectability_stats(atlas)["pearson"]
    assert s["n_rows"] == 4 and s["excluded"] == 1 and s["n_ranked_both"] == 3


def test_card_and_cli(tmp_path):
    art, out = _artifacts(tmp_path), tmp_path / "out"
    text = explore.evidence_card(explore.load_all(art), "TFA", out)
    assert "TFA" in text and "EXPLORATORY" in text
    top_row = [ln for ln in text.splitlines() if ln.startswith("| G")][0]
    assert top_row.startswith("| G1 | 3 | 3 | 1.00 |") and top_row.rstrip().endswith("| y |")   # hit in 3 of 3 datasets
    assert "## Known targets" in text and "NOTINDATA" not in text              # absent genes are never "present"
    assert (out / "card_TFA.md").exists() and (out / "card_TFA.png").exists()
    assert explore.main(["--artifacts", str(art), "--out", str(out / "cli")]) == 0
    for f in ("tf_atlas.tsv", "fig_tf_detectability.png", "fig_tf_detectability.pdf",
              "fig_promoted_edges.png", "README.md", "card_TFA.md"):
        assert (out / "cli" / f).exists(), f
