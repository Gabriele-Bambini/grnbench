"""subgroup_breakdown.py: descriptive breakdown of a paired difference by organism, assay family and
suspension (declared in EXP-20261022-GB-1 and -GB-2 before the runs). Fake rows and manifests."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks"))
import subgroup_breakdown as sb  # noqa: E402


@pytest.mark.parametrize("label, family", [
    ("10x 3' v3", "10x 3'"), ("10x 3' transcription profiling", "10x 3'"), ("10x 5' v2", "10x 5'"),
    ("10x multiome", "10x multiome"), ("Smart-seq2", "full-length"), ("Patch-seq", "full-length"),
    ("Smart-seq v4", "full-length"), ("10x Next GEM Flex v1", "other"), ("Drop-seq", "other")])
def test_assay_families(label, family):
    assert sb.assay_family(label) == family


def _manifest(tmp_path):
    rows = [("10x-3-v3/homo-sapiens/lung__aaaaaaa1", "Homo sapiens", "10x 3' v3", "cell"),
            ("10x-3-v3/homo-sapiens/gut__aaaaaaa2", "Homo sapiens", "10x 3' v3", "nucleus"),
            ("smart-seq2/mus-musculus/brain__aaaaaaa3", "Mus musculus", "Smart-seq2", "cell"),
            ("10x-5-v1/homo-sapiens/blood__aaaaaaa4", "Homo sapiens", "10x 5' v1", "cell;nucleus")]
    text = "folder\tdataset_id\torganism\tassay\tsuspension_type\n" + "".join(
        f"{f}\t{f[-8:]}-0000\t{o}\t{a}\t{s}\n" for f, o, a, s in rows)
    (tmp_path / "d.tsv").write_text(text, encoding="utf-8")
    return tmp_path / "d.tsv"


def test_real_rows_give_the_mean_and_t_interval_per_subgroup(tmp_path):
    diffs = {"aaaaaaa1": 0.3, "aaaaaaa2": 0.1, "aaaaaaa3": -0.2, "aaaaaaa4": 0.5}
    rows = pd.DataFrame({"task": [f"cxg__x__{k}" for k in diffs],
                         "m_pearson_motif_auto_epr": [1.0 + v for v in diffs.values()],
                         "m_pearson_signed_auto_epr": [1.0] * 4})
    table = sb.breakdown(sb.real_differences(rows, ("pearson_motif", "pearson_signed"), "epr", "auto"),
                         sb.load_meta(_manifest(tmp_path)))
    human = table[(table.factor == "organism") & (table.level == "human")].iloc[0]
    x = np.array([0.3, 0.1, 0.5])
    half = stats.t.ppf(0.975, 2) * x.std(ddof=1) / np.sqrt(3)
    assert human.n == 3 and human["mean"] == pytest.approx(x.mean())
    assert (human.lo, human.hi) == pytest.approx((x.mean() - half, x.mean() + half))
    mouse = table[(table.factor == "organism") & (table.level == "mouse")].iloc[0]
    assert mouse.n == 1 and np.isnan(mouse.lo)               # one dataset: no interval
    both = table[(table.factor == "suspension") & (table.level == "both")].iloc[0]
    assert both.n == 1


def test_a_scored_dataset_missing_from_the_manifest_is_an_error(tmp_path):
    rows = pd.DataFrame({"task": ["cxg__x__bbbbbbb9"], "m_pearson_motif_auto_epr": [1.2],
                         "m_pearson_signed_auto_epr": [1.0]})
    with pytest.raises(ValueError, match="bbbbbbb9"):
        sb.breakdown(sb.real_differences(rows, ("pearson_motif", "pearson_signed"), "epr", "auto"),
                     sb.load_meta(_manifest(tmp_path)))


def test_rows_without_both_scores_are_left_out_and_counted(tmp_path):
    rows = pd.DataFrame({"task": ["cxg__x__aaaaaaa1", "cxg__x__aaaaaaa2"],
                         "m_pearson_motif_auto_epr": [1.2, np.nan], "m_pearson_signed_auto_epr": [1.0, 1.0]})
    diffs = sb.real_differences(rows, ("pearson_motif", "pearson_signed"), "epr", "auto")
    assert diffs == {"aaaaaaa1": pytest.approx(0.2)}
    md = sb.render(sb.breakdown(diffs, sb.load_meta(_manifest(tmp_path))), "EPR, pearson_motif - pearson_signed, auto",
                   left_out=1)
    assert "1 dataset(s) without both scores left out" in md and "no test" in md.lower()
