"""Independent-review regressions for replay-safe standalone evaluation."""
import hashlib
import json
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from grnbench.cli_eval import main
from grnbench.diagnostics import new_document, write_diagnostics


def inputs(tmp_path, method='method'):
    edges = tmp_path / f'{method}.tsv'
    edges.write_text('source\ttarget\tweight\nTF\tG1\t2\nTF\tG2\t1\n')
    gold = tmp_path / 'gold.tsv'
    gold.write_text('source\ttarget\nTF\tG1\n')
    universe = tmp_path / 'universe.tsv'
    universe.write_text('source\ttarget\nTF\tG1\nTF\tG2\nTF\tG3\n')
    return [str(edges), '--gold', str(gold), '--universe', str(universe),
            '--out', str(tmp_path / 'out'), '--skip-nulls']


def test_eval_refuses_nonempty_output_before_loading_or_mutating(tmp_path, monkeypatch):
    args = inputs(tmp_path)
    out = tmp_path / 'out'
    out.mkdir()
    (out / 'eval_report.html').write_text('previous artifact')
    (out / 'eval_report.json').write_text('previous report')
    before = {p.name: p.read_bytes() for p in out.iterdir()}
    assert main(args + ['--no-html']) != 0
    assert {p.name: p.read_bytes() for p in out.iterdir()} == before


def test_standalone_eval_records_full_inputs_arguments_and_effective_gold(tmp_path):
    args = inputs(tmp_path)
    assert main(args + ['--seed', '17', '--top-k', '1', '--no-html']) == 0
    out = tmp_path / 'out'
    report = json.loads((out / 'eval_report.json').read_text())
    prov = report['provenance']
    assert prov['arguments']['seed'] == 17 and prov['arguments']['top_k'] == 1
    assert prov['arguments']['no_html'] is True
    for key, path in [('gold_input', tmp_path / 'gold.tsv'), ('universe_input', tmp_path / 'universe.tsv')]:
        assert prov[key]['sha256'] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert prov['edge_inputs']['method']['sha256'] == hashlib.sha256((tmp_path / 'method.tsv').read_bytes()).hexdigest()
    snapshot = out / 'effective_gold.tsv'
    assert prov['effective_gold']['snapshot']['sha256'] == hashlib.sha256(snapshot.read_bytes()).hexdigest()
    assert pd.read_csv(snapshot, sep='\t').to_records(index=False).tolist() == [('TF', 'G1')]
    assert prov['effective_gold']['identity']['count'] == 1
    assert prov['candidate_universe']['identity']['count'] == 3


def test_diagnostics_directory_fingerprints_only_the_portable_file(tmp_path):
    diagnostic_dir = tmp_path / "run"
    diagnostic_path = write_diagnostics(
        diagnostic_dir / "diagnostics.json", new_document(outcome="ok"))
    (diagnostic_dir / "processed.h5ad").write_bytes(b"not part of the diagnostic fingerprint")
    args = inputs(tmp_path)
    assert main(args + ["--diagnostics", str(diagnostic_dir), "--no-html"]) == 0
    report = json.loads((tmp_path / "out/eval_report.json").read_text())
    fingerprint = report["provenance"]["additional_inputs"]["diagnostics"]
    assert fingerprint["path"] == "diagnostics.json"
    assert fingerprint["sha256"] == hashlib.sha256(diagnostic_path.read_bytes()).hexdigest()
    assert report["provenance"]["arguments"]["diagnostics"] == "diagnostics.json"
    assert str(diagnostic_dir) not in json.dumps(report["provenance"])


@pytest.mark.parametrize('reason', ['not-scored', 'unavailable', 'unknown'])
def test_explicit_unrunnable_stability_is_failure(tmp_path, monkeypatch, reason):
    import grnbench.infer
    requested = 'absent' if reason == 'not-scored' else 'method'
    methods = [] if reason == 'unknown' else [{'name': requested, 'available': reason == 'not-scored'}]
    monkeypatch.setattr(grnbench.infer, 'list_methods', lambda: methods)
    data = tmp_path / 'data.h5ad'
    ad.AnnData(np.ones((3, 3))).write_h5ad(data)
    args = inputs(tmp_path)
    assert main(args + ['--stability-data', str(data), '--stability-methods', requested, '--no-html']) != 0
    report = json.loads((tmp_path / 'out/eval_report.json').read_text())
    assert report['status'] == 'partial-failure'
    assert any(f['stage'] == 'stability' and f['method'] == requested for f in report['failures'])
