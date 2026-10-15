"""Interface: model picker, data and results explorers, quick run (logic and pages)."""
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from grnbench.demo import make_demo, write_demo
from grnbench.ui import catalog, logic

ROOT = Path(__file__).resolve().parents[1]


def _goto(at, step):
    at.session_state["step"] = step
    at.run()


def _methods():
    base = {"origin": "builtin", "default": True, "available": True, "unavailable_reason": None, "description": "d"}
    return [{**base, "name": "pearson_signed"},
            {**base, "name": "bicor", "default": False},
            {**base, "name": "clr"},
            {**base, "name": "pearson_motif", "available": False, "unavailable_reason": "no cisTarget files"},
            {**base, "name": "mine", "origin": "python"},
            {**base, "name": "mrnet", "default": False}]


# --------------------------------------------------------------------------- #
# model picker
# --------------------------------------------------------------------------- #

def test_model_table_marks_tier_speed_and_availability():
    t = catalog.model_table(_methods(), recommended=["pearson_signed", "pearson_motif"], user_methods=["mine"])
    row = t.set_index("method")
    assert row.loc["pearson_signed", "tier"] == "default" and row.loc["bicor", "tier"] == "extra"
    assert row.loc["mine", "tier"] == "yours" and row.loc["mine", "speed"] == "unknown speed"
    assert row.loc["clr", "speed"] == "about a minute" and row.loc["pearson_signed", "speed"] == "seconds"
    assert not row.loc["pearson_motif", "available"] and row.loc["pearson_motif", "why not available"]
    assert row.loc["pearson_signed", "run"] and not row.loc["pearson_motif", "run"], "unavailable is never ticked"
    order = t["tier"].map({"default": 0, "extra": 1, "yours": 2})
    assert list(order) == sorted(order) and set(t["tier"]) <= set(catalog.TIERS)


def test_presets_tick_only_available_models():
    t = catalog.model_table(_methods(), recommended=["pearson_signed"])
    assert catalog.preset_selection(t, "Fast (seconds)") == ["pearson_signed"], "extras are not ticked by Fast"
    assert set(catalog.preset_selection(t, "Fast, including extras")) == {"pearson_signed", "bicor", "mrnet"}
    assert "clr" not in catalog.preset_selection(t, "Fast, including extras"), "clr is about a minute"
    assert "pearson_motif" not in catalog.preset_selection(t, "Everything available")
    assert catalog.preset_selection(t, "Recommended by the engine", ["pearson_signed", "pearson_motif"]) == ["pearson_signed"]
    assert catalog.preset_selection(t, "None") == []
    assert set(catalog.PRESETS) >= {"Fast (seconds)", "Fast, including extras", "None"}


def test_the_default_tier_is_what_the_engine_calls_default():
    """``list_methods()`` marks the default set with ``default`` once the extra models exist; before that
    every method is in the default set."""
    from grnbench.infer import list_methods
    builtin = [m for m in list_methods() if m["origin"] == "builtin"]
    t = catalog.model_table(builtin).set_index("method")
    for m in builtin:
        assert (t.loc[m["name"], "tier"] == "default") == bool(m.get("default", True)), m["name"]


def test_every_registered_method_is_described_or_listed_as_unknown():
    from grnbench.infer import list_methods
    builtin = [m for m in list_methods() if m["origin"] == "builtin"]     # a test may have registered a plugin
    t = catalog.model_table(builtin)
    assert set(t["method"]) == {m["name"] for m in builtin}
    assert (t["what it is"].str.len() > 0).all()
    unknown = t.loc[t["speed"] == "unknown speed", "method"].tolist()
    assert unknown == [], f"add these to MODEL_INFO: {unknown}"


# --------------------------------------------------------------------------- #
# looking at the data
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def demo_adata():
    adata, tfs, gold = make_demo(n_cells=400, n_genes=300, n_tfs=8, targets_per_tf=8)
    return adata, tfs, gold


def test_embedding_frame_keeps_annotations_and_subsamples(demo_adata):
    adata, _, _ = demo_adata
    f = logic.embedding_frame(adata, "PCA", max_cells=150)
    assert len(f) == 150 and f.attrs == {"method": "PCA", "n_cells_total": 400}
    assert {"x", "y", "batch", "cell_type", "detected genes"} <= set(f.columns)
    assert f.index.is_unique and set(f.index) <= set(adata.obs_names)
    again = logic.embedding_frame(adata, "PCA", max_cells=150)
    pd.testing.assert_frame_equal(f, again)           # same seed, same cells, same map
    opts = logic.color_choices(f)
    assert opts.index("batch") < opts.index("detected genes"), "categories come before numbers"
    assert "x" not in opts and "y" not in opts


def test_embedding_errors_are_clear(demo_adata, monkeypatch):
    adata, _, _ = demo_adata
    with pytest.raises(ValueError, match="unknown embedding"):
        logic.embedding_frame(adata, "TSNE")
    monkeypatch.setitem(sys.modules, "umap", None)
    with pytest.raises(RuntimeError, match="umap-learn"):
        logic.embedding_frame(adata, "UMAP", max_cells=100)
    tiny = ad.AnnData(np.ones((2, 5), dtype="float32"))
    with pytest.raises(ValueError, match="too small"):
        logic.embedding_frame(tiny)


def test_embedding_figure_categorical_and_continuous(demo_adata):
    adata, _, _ = demo_adata
    f = logic.embedding_frame(adata, "PCA", max_cells=200)
    cat = logic.embedding_figure(f, "batch")
    assert len(cat.data) == f["batch"].nunique()
    cont = logic.embedding_figure(f, "detected genes")
    assert len(cont.data) == 1 and cont.layout.coloraxis.colorscale is not None
    assert len(logic.embedding_figure(f, None).data) == 1


def test_gene_summary_and_expression_frame(demo_adata):
    adata, tfs, _ = demo_adata
    s = logic.gene_summary(adata, tfs)
    assert len(s) == adata.n_vars and s["detection"].between(0, 1).all()
    assert s["detection"].is_monotonic_decreasing and s.attrs["detection_source"] == "X"
    assert s.loc[s["TF"], "gene"].str.upper().isin([t.upper() for t in tfs]).all() and s["TF"].sum() == len(tfs)
    genes = [tfs[0], adata.var_names[40], "NOT_A_GENE"]
    f = logic.expression_frame(adata, genes, "cell_type", max_cells=100)
    assert set(f["gene"]) == {genes[0], genes[1]} and len(f) == 2 * 100
    assert {"cell", "gene", "expression", "cell_type"} <= set(f.columns)
    assert len(logic.expression_frame(adata, ["NOPE"]).index) == 0
    fig = logic.expression_figure(f, "cell_type")
    assert len(fig.data) >= 1


def test_sparse_summaries_give_the_numbers_of_the_dense_computation(demo_adata):
    """B2 (PR #71 review): the summaries stopped making the matrix dense; the numbers are the ones the dense code gave."""
    from sklearn.decomposition import PCA
    from grnbench.infer import prepare_adata

    adata, tfs, _ = demo_adata

    def dense(a):
        return prepare_adata(a)[0].X.toarray().astype(np.float32)

    X = dense(adata)
    mean = X.mean(axis=0)
    s = logic.gene_summary(adata, tfs).set_index("gene").reindex(adata.var_names.astype(str))
    assert s["mean"].to_numpy() == pytest.approx(mean, rel=1e-5, abs=1e-6)
    assert s["variance"].to_numpy() == pytest.approx((X ** 2).mean(axis=0) - mean ** 2, rel=1e-4, abs=1e-5)
    X = dense(logic.subsample_cells(adata, 150, 0))
    top = np.argsort(-X.var(axis=0), kind="stable")[:100]
    D = X[:, top] - X[:, top].mean(axis=0)
    xy = PCA(n_components=20, random_state=0).fit_transform(D)[:, :2]
    f = logic.embedding_frame(adata, "PCA", max_cells=150, n_top=100)
    # same genes in the same order; float32 PCA moves by ~1e-4 with the memory layout alone (the old code run on a
    # C-ordered copy of its own input differs from itself by that much)
    assert f[["x", "y"]].to_numpy() == pytest.approx(xy, abs=1e-3)
    X = dense(logic.subsample_cells(adata, 100, 0))
    genes = [str(adata.var_names[5]), str(adata.var_names[40])]
    e = logic.expression_frame(adata, genes, max_cells=100)
    assert e["expression"].to_numpy() == pytest.approx(np.r_[X[:, 5], X[:, 40]])


def test_the_explorer_never_makes_a_whole_sparse_matrix_dense(monkeypatch):
    """B2 (PR #71 review): opening Explore data called toarray() on the whole matrix, 8.7 GB of float32 for a
    247,914 x 8,749 Perturb-seq dataset, and the map and the violins made every gene of the sampled cells dense.
    Here any conversion that spans every cell or every gene fails."""
    from scipy import sparse

    n_cells, n_genes = 20000, 3000
    X = sparse.random(n_cells, n_genes, density=0.01, format="csr", dtype=np.float32, random_state=0)
    X.data = np.ceil(X.data * 5)              # integer counts, so the raw-count normalisation runs too
    a = ad.AnnData(X)
    a.var_names = [f"G{i}" for i in range(n_genes)]
    a.obs["group"] = np.where(np.arange(n_cells) % 2, "a", "b")
    made_dense = []
    for cls in (sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix, sparse.csr_array, sparse.csc_array,
                sparse.coo_array):      # todense() goes through toarray()
        def guard(self, *args, _toarray=cls.toarray, **kwargs):
            if self.shape[0] == n_cells or self.shape[1] == n_genes:
                raise AssertionError(f"a whole axis made dense: {self.shape}")
            made_dense.append(self.shape)
            return _toarray(self, *args, **kwargs)
        monkeypatch.setattr(cls, "toarray", guard)
    s = logic.gene_summary(a)
    assert len(s) == n_genes and made_dense == []
    f = logic.embedding_frame(a, "PCA")
    assert len(f) == logic.MAX_EMBED_CELLS and f.attrs["n_cells_total"] == n_cells
    assert made_dense == [(logic.MAX_EMBED_CELLS, 2000)], "the sampled cells x the 2,000 most variable genes"
    made_dense.clear()
    e = logic.expression_frame(a, ["G1", "G2", "G3"], "group")
    assert len(e) == 3 * 5000 and made_dense == [(5000, 3)]


def _counts_adata(n=60, g=12, seed=0):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(rng.gamma(1.0, 1.0, size=g), size=(n, g)).astype("float32")
    a = ad.AnnData(counts.copy())
    a.var_names = [f"G{i}" for i in range(g)]
    return a, counts


def test_detection_is_read_from_the_counts_layer_not_from_a_centred_matrix():
    a, counts = _counts_adata()
    expected = (counts > 0).mean(axis=0)
    # cleaned data as the pipeline leaves it: centred X (no zero means "not seen") plus the counts layer
    a.layers["counts"] = counts
    a.X = (np.log1p(counts) - np.log1p(counts).mean(axis=0)).astype("float32")
    det, source = logic.gene_detection(a)
    assert source == "counts layer" and det == pytest.approx(expected)
    s = logic.gene_summary(a)
    assert s.attrs["detection_source"] == "counts layer"
    assert s.set_index("gene")["detection"].reindex(a.var_names).to_numpy() == pytest.approx(expected)
    assert logic.cell_detected_genes(a) == pytest.approx((counts > 0).sum(axis=1))
    # X = Pearson residuals are negative too
    a.X = (counts - counts.mean(axis=0)).astype("float32")
    assert logic.gene_detection(a)[0] == pytest.approx(expected)


def test_detection_is_nan_not_a_guess_when_nothing_says_which_entries_are_zero():
    a, counts = _counts_adata()
    a.X = (counts - counts.mean(axis=0)).astype("float32")            # centred, no counts layer
    det, source = logic.gene_detection(a)
    assert np.isnan(det).all() and "negative" in source
    s = logic.gene_summary(a)
    assert s["detection"].isna().all() and len(s) == a.n_vars          # still lists genes and the other columns
    assert np.isnan(logic.cell_detected_genes(a)).all()
    frame = logic.embedding_frame(a, "PCA", max_cells=40)
    assert frame["detected genes"].isna().all()
    assert "detected genes" not in logic.color_choices(frame), "an all-NaN column is not offered"


def test_stage_choices_and_loading(demo_adata, tmp_path):
    adata, _, _ = demo_adata
    path = tmp_path / "clean.h5ad"
    adata[:50].copy().write_h5ad(path)

    class D:  # the only attribute load_stage reads
        pass
    d = D()
    d.adata = adata
    recs = [{"context": "all", "processed": str(path)}, {"context": "bad", "processed": None}]
    assert logic.stage_choices(recs) == ["As loaded (raw counts)", "Cleaned: all"]
    assert logic.load_stage(d, recs, "As loaded (raw counts)").n_obs == 400
    assert logic.load_stage(d, recs, "Cleaned: all").n_obs == 50
    with pytest.raises(KeyError):
        logic.load_stage(d, recs, "Cleaned: bad")


# --------------------------------------------------------------------------- #
# looking at the results
# --------------------------------------------------------------------------- #

def _edges(rows):
    return pd.DataFrame(rows, columns=["source", "target", "weight"])


def test_edges_table_ranks_filters_and_flags_reference():
    e = _edges([("tfa", "G1", 0.9), ("TFA", "G2", 0.8), ("TFB", "G1", 0.7), ("TFB", "G3", 0.6)])
    ref = {("TFA", "G1"), ("TFB", "G3")}
    t = logic.edges_table(e, ref)
    assert t["rank"].tolist() == [1, 2, 3, 4] and t["in reference"].tolist() == [True, False, False, True]
    assert logic.edges_table(e, ref, tf="tfb")["target"].tolist() == ["G1", "G3"]
    assert set(logic.edges_table(e, ref, gene="g1")["rank"]) == {1, 3}
    assert len(logic.edges_table(e, top=2)) == 2 and "in reference" not in logic.edges_table(e)


def test_edges_table_gives_tied_edges_one_rank_whatever_the_row_order():
    rows = [("T", "G1", 5.0), ("T", "G2", 3.0), ("T", "G3", 3.0), ("T", "G4", 3.0), ("T", "G5", 1.0)]
    a = logic.edges_table(_edges(rows))
    b = logic.edges_table(_edges(rows[::-1]))
    assert a["rank"].tolist() == [1, 2, 2, 2, 5] and b["rank"].tolist() == a["rank"].tolist()


def test_network_figure_draws_at_most_top_n_edges_and_marks_reference():
    rng = np.random.default_rng(0)
    e = _edges([(f"TF{i % 4}", f"G{i}", 1 - i / 100) for i in range(60)])
    ref = {("TF0", "G0"), ("TF1", "G1")}
    fig = logic.network_figure(e, ref, top_n=20)
    lines = [t for t in fig.data if t.mode == "lines"]
    n_segments = sum(len([v for v in t.x if v is not None]) // 2 for t in lines)
    assert n_segments == 20 and len([a for a in fig.layout.annotations if a.showarrow]) == 20, "one arrow per edge"
    assert {t.name for t in fig.data if t.mode and "markers" in t.mode} == {"TF", "gene"}
    assert any(t.name == "in reference" for t in lines)
    only = logic.network_figure(e, ref, top_n=50, tf="TF2")
    assert all(t.name == "gene" or t.name == "TF" or t.mode == "lines" for t in only.data)
    empty = logic.network_figure(e.iloc[0:0], ref)
    assert empty.layout.annotations[0].text == "no edges to draw"
    del rng


def test_network_figure_keeps_the_direction_and_both_edges_of_a_reciprocal_pair():
    """B4 (PR #71 review): the network was an undirected graph, so TF1 -> TF2 and TF2 -> TF1 were one line with
    no direction (2 input edges, 1 line segment, 0 arrows)."""
    fig = logic.network_figure(_edges([("TF1", "TF2", 0.9), ("TF2", "TF1", 0.8)]), {("TF1", "TF2")})
    segments = [((t.x[i], t.y[i]), (t.x[i + 1], t.y[i + 1])) for t in fig.data if t.mode == "lines"
                for i in range(0, len(t.x), 3)]
    assert len(segments) == 2 and segments[0] != segments[1] and segments[0] != segments[1][::-1]
    arrows = [a for a in fig.layout.annotations if a.showarrow]
    assert len(arrows) == 2
    node = {h.split(":")[0]: np.array([x, y]) for t in fig.data if t.mode and "markers" in t.mode
            for x, y, h in zip(t.x, t.y, t.hovertext)}

    def nearest(x, y):
        return min(node, key=lambda n: np.hypot(*(node[n] - (x, y))))
    drawn = set()
    for a in arrows:
        head = nearest(a.x, a.y)
        tail = next(n for n in node if n != head)
        assert np.dot(np.array([a.x - a.ax, a.y - a.ay]), node[head] - node[tail]) > 0, "the arrow points at its target"
        drawn.add((tail, head, a.arrowcolor))
    assert drawn == {("TF1", "TF2", "#009E73"), ("TF2", "TF1", "#B9B9B9")}, "only TF1 -> TF2 is in the reference"


def test_every_arrow_of_a_dense_network_keeps_its_head_wherever_two_markers_do_not_overlap():
    """B4 follow-up (verifier on the PR #71 fix): each arrow ran from the middle of its edge to the target, half the
    edge long. Plotly drops an arrow, line and head, shorter than ``standoff + startstandoff + head backoff`` (1.3 x
    arrowsize x arrowwidth for arrowhead 2, plotly.js ``draw_arrow_head.js``) and first cuts its tail to the edge of
    the empty label box (4 x 4 px, so at most 2*sqrt(2) px). Two close TF squares (size 30, standoff 15) then showed a
    plain undirected segment: 2 of 150 heads lost at 1,400 px, 5 at 700 px, among them the only reference edge.

    The test is width-free: an arrow spanning a fraction ``f`` of an edge whose markers are ``d`` px apart is drawn iff
    ``f * d - trim >= standoff + startstandoff + backoff``. Two markers do not overlap when ``d >= (s_u + s_v) / 2``,
    so at that ``d`` every arrow must still be drawn, whatever the plot width."""
    tfs = [f"TF{i}" for i in range(6)]
    tf_tf = [("TF0", "TF1"), ("TF1", "TF0"), ("TF2", "TF3"), ("TF3", "TF2"),           # two reciprocal pairs
             ("TF0", "TF2"), ("TF0", "TF3"), ("TF0", "TF4"), ("TF0", "TF5"), ("TF1", "TF2"), ("TF1", "TF4"),
             ("TF2", "TF5"), ("TF3", "TF4"), ("TF4", "TF5"), ("TF5", "TF1")]
    tf_gene = [(tfs[i % 6], f"G{(7 * i) % 97}") for i in range(136)]       # 97 genes, some shared by two TFs
    assert len(set(tf_gene)) == 136
    rows = [(s, t, 2.0 - i / 1000) for i, (s, t) in enumerate(tf_tf)] + \
           [(s, t, 1.0 - i / 1000) for i, (s, t) in enumerate(tf_gene)]
    fig = logic.network_figure(_edges(rows), {("TF0", "TF1"), *tf_gene[:10]}, top_n=150)
    arrows = [a for a in fig.layout.annotations if a.showarrow]
    assert len(arrows) == 150
    node, size = {}, {}
    for t in fig.data:
        if t.mode and "markers" in t.mode:
            for x, y, h, s in zip(t.x, t.y, t.hovertext, t.marker.size):
                node[h.split(":")[0]], size[h.split(":")[0]] = np.array([x, y]), s

    def nearest(x, y):
        return min(node, key=lambda n: np.hypot(*(node[n] - (x, y))))
    backoff = {2: 1.3}          # plotly.js arrow_paths.js, arrowhead 2
    seen = set()
    for a in arrows:
        u, v = nearest(a.ax, a.ay), nearest(a.x, a.y)
        assert u in tfs and u != v, "the tail is at the regulator"
        seen.add((u, v))
        edge = node[v] - node[u]
        f = np.hypot(a.x - a.ax, a.y - a.ay) / np.hypot(*edge)
        assert np.dot([a.x - a.ax, a.y - a.ay], edge) > 0, "the arrow points at its target"
        need = (a.standoff + (a.startstandoff or 0) + backoff[a.arrowhead] * a.arrowsize * a.arrowwidth
                + 2 * np.sqrt(2))
        assert f * (size[u] + size[v]) / 2 >= need, (u, v, f, size[u], size[v], need)
        assert a.standoff >= size[v] / 2 - 1e-9, "the head stops at the rim of the target, not inside it"
    assert seen == set(tf_tf) | set(tf_gene), "one arrow per edge, reciprocal pairs included"


def test_overlap_matrix_is_jaccard_of_the_top_k():
    a = _edges([("T", "G1", 3), ("T", "G2", 2), ("T", "G3", 1)])
    b = _edges([("T", "G2", 3), ("T", "G3", 2), ("T", "G1", 1)])
    m = logic.overlap_matrix({"a": a, "b": b, "same": a}, k=2)
    assert m.loc["a", "a"] == 1 and m.loc["a", "same"] == 1
    assert m.loc["a", "b"] == pytest.approx(1 / 3) and m.loc["b", "a"] == pytest.approx(1 / 3)
    assert logic.overlap_matrix({"a": a, "b": b}, k=3).loc["a", "b"] == 1.0
    assert len(logic.overlap_figure(m, 2).data) == 1


def test_overlap_keeps_a_tie_block_at_the_cut_whole_and_ignores_row_order():
    rows = [("T", "G1", 9.0), ("T", "G2", 5.0), ("T", "G3", 5.0), ("T", "G4", 5.0), ("T", "G5", 1.0)]
    other = _edges([("T", "G1", 9.0), ("T", "G2", 5.0), ("T", "G3", 5.0), ("T", "G4", 5.0), ("T", "G6", 1.0)])
    forward = logic.overlap_matrix({"a": _edges(rows), "b": other}, k=2)
    backward = logic.overlap_matrix({"a": _edges(rows[::-1]), "b": other.iloc[::-1]}, k=2)
    assert forward.loc["a", "b"] == 1.0, "k=2 cuts inside a block of three ties: all of the block is kept"
    pd.testing.assert_frame_equal(forward, backward)
    assert forward.attrs["sizes"] == {"a": 4, "b": 4}, "the realised set size is reported, not k"


def test_overlap_skips_a_model_it_cannot_compare_and_says_why():
    good = _edges([("T", "G1", 3.0), ("T", "G2", 2.0), ("T", "G3", 1.0)])
    bad = _edges([("T", "G1", np.nan), ("T", "G2", 2.0), ("T", "G3", 1.0)])
    m = logic.overlap_matrix({"good": good, "bad": bad, "empty": good.iloc[0:0]}, k=2)
    assert list(m.index) == ["good"] and "non-finite" in m.attrs["skipped"]["bad"]


def test_tf_recovery_and_detail_separate_a_good_ranking_from_a_bad_one():
    ref = {("T1", f"G{i}") for i in range(4)} | {("T2", f"G{i}") for i in range(2)}
    good = _edges([("T1", f"G{i}", 10 - i) for i in range(4)] + [("T1", f"X{i}", 1) for i in range(5)]
                  + [("T2", f"G{i}", 10 - i) for i in range(2)])
    bad = _edges([("T1", f"X{i}", 10 - i) for i in range(5)] + [("T1", f"G{i}", 0.1) for i in range(4)]
                 + [("T2", "X1", 5), ("T2", "G0", 0.1), ("T2", "G1", 0.1)])
    t = logic.tf_recovery_table({"good": good, "bad": bad}, ref, k_per_tf=4)
    assert t.loc["T1", "good"] == 1.0 and t.loc["T1", "bad"] == 0.0
    assert t.loc["T1", "reference targets"] == 4 and t.index[0] == "T1"
    assert len(logic.tf_recovery_figure(t, 4).data) == 1
    top, missed, n_out = logic.tf_detail(bad, ref, "T1", k=3)
    assert top["in reference"].tolist() == [False, False, False] and missed == ["G0", "G1", "G2", "G3"] and n_out == 0
    assert top["rank for this TF"].tolist() == [1, 2, 3]
    top, missed, _ = logic.tf_detail(good, ref, "t1", k=2)
    assert top["in reference"].all() and missed == ["G2", "G3"]


def test_tf_recovery_counts_only_what_is_in_the_data_and_marks_no_edges_as_nan():
    ref = {("T1", "G1"), ("T1", "G2"), ("T1", "G3"), ("T1", "OUT"), ("GONE", "G1"), ("T2", "G1")}
    genes = ["T1", "T2", "G1", "G2", "G3"]
    full = _edges([("T1", "G1", 3), ("T1", "G2", 2), ("T1", "G3", 1), ("T1", "T2", 0.5)])
    no_t2 = full.copy()
    t = logic.tf_recovery_table({"m": full, "n": no_t2}, ref, genes, k_per_tf=2)
    assert list(t.index) == ["T1"], "TF absent from the data and TF with no ranked edges are not rows"
    assert t.loc["T1", "reference targets"] == 3, "OUT is not a gene of the data: not in the denominator"
    assert t.loc["T1", "m"] == pytest.approx(2 / 3)
    only_t1 = _edges([("T1", "G1", 3), ("T1", "G2", 2)])
    only_t2 = _edges([("T2", "G1", 3)])
    t = logic.tf_recovery_table({"a": only_t1, "b": only_t2}, ref, genes, k_per_tf=2)
    assert np.isnan(t.loc["T2", "a"]) and np.isnan(t.loc["T1", "b"]), "no stored edge is NaN, not a zero recovery"
    assert t.loc["T2", "b"] == 1.0
    top, missed, n_out = logic.tf_detail(full, ref, "T1", k=1, genes=genes)
    assert missed == ["G2", "G3"] and n_out == 1


def test_tf_recovery_folds_case_and_does_not_depend_on_row_order_at_ties():
    ref = {("tf1", "g1"), ("tf1", "g2")}
    rows = [("TF1", "G1", 2.0), ("TF1", "G2", 1.0), ("TF1", "G3", 1.0), ("TF1", "G4", 1.0)]
    genes = ["TF1", "G1", "G2", "G3", "G4"]
    a = logic.tf_recovery_table({"m": _edges(rows)}, ref, genes, k_per_tf=2)
    b = logic.tf_recovery_table({"m": _edges(rows[::-1])}, ref, genes, k_per_tf=2)
    pd.testing.assert_frame_equal(a, b)
    # G1 is certain; one slot is left in a block of three tied edges of which one is a hit: 1 + 1/3
    assert a.loc["TF1", "m"] == pytest.approx((1 + 1 / 3) / 2)


def test_reference_for_scoring_is_one_rule_for_every_page():
    full = {("A", "B"), ("A", "C"), ("D", "E")}
    sel, ev = {("A", "B")}, {("A", "C"), ("D", "E")}
    screen = {"reference": {"selection": sel, "evaluation": ev}}
    pairs, label, warn = logic.reference_for_scoring(None, full)
    assert pairs == full and "full reference" in label and warn is None
    pairs, label, warn = logic.reference_for_scoring(screen, full)
    assert pairs == ev and "evaluation half" in label and warn is None
    pairs, label, warn = logic.reference_for_scoring(None, None)
    assert pairs is None and warn is None
    # an upload that contains the selection half is flagged; one that does not is not
    pairs, label, warn = logic.reference_for_scoring(screen, full, {("a", "b"), ("x", "y")}, "mine.tsv")
    assert pairs == {("a", "b"), ("x", "y")} and "mine.tsv" in label and "1 of the uploaded" in warn
    assert logic.reference_for_scoring(screen, full, {("x", "y")}, "mine.tsv")[2] is None
    # the cleaning was chosen on a reference that was too small to hold anything out: nothing fair to score on
    pairs, _, warn = logic.reference_for_scoring({"reference": {"selection": sel, "evaluation": set()}}, full)
    assert pairs is None and "too small" in warn


def test_check_edges_rejects_what_the_pages_cannot_draw_or_score():
    genes = ["T", "G1", "G2"]
    logic._check_edges(_edges([("T", "G1", 1.0), ("T", "G2", 0.5)]), [g.lower() for g in genes])   # names fold case
    for bad, why in [(_edges([("T", "G1", np.nan)]), "finite"), (_edges([("T", "G1", 1), ("T", "G1", 2)]), "duplicate"),
                     (_edges([("T", "T", 1.0)]), "self-loop"), (_edges([("T", "OTHER", 1.0)]), "not genes of the data"),
                     (_edges([]), "empty")]:
        with pytest.raises(ValueError, match=why):
            logic._check_edges(bad, genes)


def test_pr_figure_has_one_curve_per_model_and_the_random_line():
    genes = [f"G{i}" for i in range(12)] + ["T1", "T2"]
    ref = {("T1", "G0"), ("T1", "G1"), ("T2", "G2")}
    good = _edges([("T1", "G0", 9), ("T1", "G1", 8), ("T2", "G2", 7)]
                  + [(s, g, 1 - j / 100) for j, (s, g) in enumerate((s, g) for s in ("T1", "T2") for g in genes[3:12])])
    rand = good.sample(frac=1, random_state=0).assign(weight=lambda d: np.linspace(1, 0, len(d)))
    fig = logic.pr_figure({"good": good, "random": rand}, ref, genes, tf_list=["T1", "T2"])
    assert [t.name for t in fig.data] == ["good", "random"]
    assert fig.layout.shapes, "the random baseline line"
    assert fig.data[0].y[0] >= fig.data[1].y[0]


def test_pr_curve_without_listing_the_universe_is_the_one_on_the_listed_universe():
    """B2 (PR #71 review): the curve was scored on a frame of every regulator x gene pair (76 million at 8,749 genes
    without a TF list). It is now computed from the scored edges and two counts; it must equal ``evaluate_ranking``
    on the listed universe exactly: ties, reference edges left unscored, more than 200 points, with and without TFs."""
    from grnbench.eval.metrics import GoldIndex, evaluate_ranking

    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(100)] + ["T1", "T2", "t3"]
    tfs = ["T1", "T2", "T3"]
    ref = {("T1", f"G{i}") for i in range(0, 40, 3)} | {("t2", f"g{i}") for i in range(41, 100, 4)} | {("T1", "T1"),
                                                                                                   ("T3", "OUT")}
    scored = [(s, g) for s in ("T1", "T2", "t3") for g in genes if g != s]
    ties = _edges([(s, t, float(rng.integers(0, 5))) for s, t in scored[:120]])          # few score levels
    many = _edges([(s, t, float(w)) for (s, t), w in zip(scored, rng.permutation(len(scored)))])   # every pair, > 200
    gi = GoldIndex.coerce(ref)
    for tf_list in (tfs, None):
        universe = logic._universe(genes, tf_list)
        sources = [g for g in genes if g.upper() in {t.upper() for t in tf_list}] if tf_list else genes
        listed = pd.DataFrame([(s, g) for s in sources for g in genes if s != g], columns=["source", "target"])
        assert universe["n"] == len(listed)
        for e in (ties, many):
            expected = evaluate_ranking(e, gi, with_curve=True, universe=listed)["pr_curve"]
            assert logic._pr_curve_on_universe(e, gi, universe) == expected
    assert len(evaluate_ranking(many, gi, universe=listed)["pr_curve"]["recall"]) == 200, "the thinning ran"
    with pytest.raises(ValueError, match="outside the explicit universe"):
        logic._pr_curve_on_universe(_edges([("T1", "NOT_A_GENE", 1.0)]), gi, logic._universe(genes, tfs))
    with pytest.raises(ValueError, match="outside the explicit universe"):
        logic._pr_curve_on_universe(_edges([("G1", "G2", 1.0)]), gi, logic._universe(genes, tfs))


# --------------------------------------------------------------------------- #
# one-click run
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def quick(tmp_path_factory):
    d = tmp_path_factory.mktemp("quick")
    paths = write_demo(d / "demo", n_cells=900, n_genes=700, n_tfs=12, targets_per_tf=10)
    tfs = paths["tfs"].read_text().split()
    gold = pd.read_csv(paths["gold"], sep="\t")
    ref = set(zip(gold.source, gold.target))
    ticks = []
    out = logic.quick_run(paths["counts"], d / "runs", ["pearson_signed", "spearman_signed"], tf_list=tfs,
                          reference=ref, progress=lambda f, m: ticks.append((f, m)))
    return d, paths, tfs, ref, out, ticks


def test_quick_run_cleans_with_the_baseline_runs_and_scores(quick):
    d, _, _, ref, out, ticks = quick
    from grnbench.screen import BASELINE_CHOICES
    choices = out["plan"].choices
    assert {k: choices[k] for k in ("normalization", "context_scope", "batch", "doublets")} == {
        k: BASELINE_CHOICES[k] for k in ("normalization", "context_scope", "batch", "doublets")}
    assert out["screen"] is None
    assert all(r["error"] is None for r in out["records"]) and out["records"]
    assert set(out["runs"]["method"]) == {"pearson_signed", "spearman_signed"} and out["runs"]["edges"].notna().all()
    assert set(out["bench"]["method"]) == {"pearson_signed", "spearman_signed"} and (out["bench"]["AUPRC ratio"] > 1).all()
    assert out["reference_used"] == ref and Path(out["workdir"], "clean").exists() and out["reference_warning"] is None
    import json
    for row in out["runs"].itertuples():
        prov = json.loads(Path(row.edges).with_name(f"{row.method}.provenance.json").read_text())
        assert prov["method"] == row.method and prov["context"] == row.context and prov["n_edges"] == row.n_edges
    fracs = [f for f, _ in ticks]
    assert fracs == sorted(fracs) and fracs[-1] == 1.0 and ticks[0][1].startswith("Reading")


def test_quick_run_without_a_reference_runs_but_does_not_score(quick):
    d, paths, tfs, _, _, _ = quick
    out = logic.quick_run(paths["counts"], d / "runs2", ["pearson_signed"], tf_list=tfs)
    assert out["bench"] is None and out["runs"]["edges"].notna().all()


# --------------------------------------------------------------------------- #
# pages
# --------------------------------------------------------------------------- #

def _app_with_results(quick):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    _, paths, tfs, ref, out, _ = quick
    at = AppTest.from_file(str(ROOT / "src" / "grnbench" / "ui" / "app.py"), default_timeout=120)
    at.session_state["workdir"] = out["workdir"]
    at.session_state["data_path"] = str(paths["counts"])
    at.session_state["tfs"] = tfs
    at.session_state["reference"] = ref
    at.session_state["diag"] = out["diag"]
    at.session_state["plan"] = out["plan"]
    at.session_state["records"] = out["records"]
    at.session_state["runs"] = out["runs"]
    at.session_state["bench"] = out["bench"]
    at.session_state["user_choices"] = {}
    return at


def test_quick_run_page_renders_with_a_result(quick):
    at = _app_with_results(quick)
    at.session_state["quick_done"] = True
    at.run()
    _goto(at, "Quick run")
    assert not at.exception and not at.error
    assert any("Analysis complete" in s.value for s in at.success)


def test_explore_pages_render_on_a_finished_run(quick):
    at = _app_with_results(quick)
    at.run()
    for step in ("6 · Models", "Explore data", "Explore results"):
        _goto(at, step)
        assert not at.exception, (step, [e.value for e in at.exception])
        assert not at.error, (step, [e.value for e in at.error])
        if step == "Explore data":      # the before/after maps are sampled and embedded separately
            captions = " ".join(c.value for c in at.caption)
            assert "own random sample" in captions and "Same cells" not in captions, captions
    assert any("Explore the networks" in m.value for m in at.markdown), "the screen has its header"
    assert any("not a result" in i.value for i in at.info), "the exploratory banner is on the page"
    assert len(at.tabs) == 5


def test_explore_results_page_survives_a_model_with_unreadable_edges(quick, tmp_path):
    at = _app_with_results(quick)
    broken = quick[4]["runs"].copy()
    bad = tmp_path / "broken.edges.tsv.gz"
    bad.write_bytes(b"not a gzip file")
    broken.loc[broken.index[0], "edges"] = str(bad)
    at.session_state["runs"] = broken
    at.run()
    _goto(at, "Explore results")
    assert not at.exception and not at.error
    assert any("could not be read" in w.value for w in at.warning)


def test_no_tf_list_takes_per_target_regression_out_of_the_fast_presets():
    """UI-07: lasso is 'seconds' with a short TF list (3 s at 1,400 cells x 1,000 genes with 30 TFs) but
    about 27 minutes when every gene is a regulator, so without a TF list it must not be in Fast."""
    base = {"origin": "builtin", "default": True, "available": True, "unavailable_reason": None, "description": "d"}
    ms = [{**base, "name": n} for n in ("pearson_signed", "lasso", "community", "tigress")]
    with_tfs = catalog.model_table(ms).set_index("method")
    assert with_tfs.loc["lasso", "speed"] == "seconds"
    assert catalog.preset_selection(catalog.model_table(ms), "Fast (seconds)") == ["pearson_signed", "lasso"]
    no_tfs = catalog.model_table(ms, all_regulators=True)
    assert catalog.preset_selection(no_tfs, "Fast (seconds)") == ["pearson_signed"]
    assert set(no_tfs.set_index("method").loc[["lasso", "community", "tigress"], "speed"]) == {"minutes or more"}


def test_models_page_starts_from_fast_and_warns_when_a_slow_model_is_ticked(quick):
    """Step 6 must not pre-tick the engine's slow tree ensembles (one GENIE3 job is tens of minutes on a
    631-cell context, and the page is busy until every job ends): it opens like the quick run, and ticking a
    slow model says so before Run."""
    at = _app_with_results(quick)
    at.run()
    _goto(at, "6 · Models")
    assert not at.exception
    assert [r.value for r in at.radio if r.label == "Start from"] == ["Fast (seconds)"]
    assert not at.warning
    at.radio(key="models_preset").set_value("Everything available").run()
    assert any("minutes or more" in w.value for w in at.warning)
