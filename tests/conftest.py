"""Shared fixtures: synthetic count datasets with known properties, and a fake Rscript."""
from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy import sparse


@pytest.fixture(scope="session")
def _r_probe_cache(tmp_path_factory) -> Path:
    return tmp_path_factory.mktemp("r_probe_cache")


@pytest.fixture(autouse=True)
def _keep_r_probe_cache_out_of_the_user_cache(_r_probe_cache, monkeypatch):
    """No test writes R probe results into the person's real cache folder."""
    monkeypatch.setenv("GRNBENCH_CACHE_DIR", str(_r_probe_cache))


# A stand-in for Rscript (plain Python, no R needed). It answers the bridge's two R files:
# _probe.R prints the JSON line built from FAKE_R_PACKAGES / FAKE_R_BUILT / FAKE_R_LIBRARY (folders joined
# with os.pathsep) / FAKE_R_VERSION; _wrapper.R reads params.json from its working folder and acts on
# FAKE_R_MODE = ok | fail | sleep | no_output | bad_session (sleep starts a sleeping child too and writes its
# pid to FAKE_R_CHILD_PID); like the real wrapper it takes no arguments of its own. Every call appends one JSON
# line to the file named by FAKE_R_CALLS: the R file, the arguments, the first PATH entries, R_LIBS_USER and
# the working folder, so tests can count R starts.
FAKE_R = r'''
import json, os, subprocess, sys, time
from pathlib import Path

args = sys.argv[1:]
vanilla = args[:1] == ["--vanilla"]
if vanilla:
    args = args[1:]
script = Path(args[0]).name if args else ""
calls = os.environ.get("FAKE_R_CALLS")
if calls:
    with open(calls, "a", encoding="utf-8") as fh:
        fh.write(json.dumps({"script": script, "argv": sys.argv[1:], "vanilla": vanilla,
                             "path": os.environ.get("PATH", "").split(os.pathsep)[:8],
                             "r_libs_user": os.environ.get("R_LIBS_USER"), "cwd": os.getcwd()}) + "\n")
version = os.environ.get("FAKE_R_VERSION", "4.5.3")
if script == "_probe.R":
    have = json.loads(os.environ.get("FAKE_R_PACKAGES", "{}"))
    built = json.loads(os.environ.get("FAKE_R_BUILT", "{}"))
    libs = [p for p in os.environ.get("FAKE_R_LIBRARY", "").split(os.pathsep) if p]
    print("a chatty site profile line before the JSON")
    print(json.dumps({"version": version, "version_string": f"R version {version} (fake)", "library": libs,
                      "packages": {p: have.get(p) for p in args[1:]},
                      "built": {p: built.get(p, version) if have.get(p) else None for p in args[1:]}}))
    sys.exit(0)
if script == "_wrapper.R" and len(args) == 1:      # like the real wrapper: everything is in ./params.json
    p = json.loads(Path("params.json").read_text(encoding="utf-8"))
    workdir = Path(p["workdir"])
    mode = os.environ.get("FAKE_R_MODE", "ok")
    print(f"fake R runs {Path(p['script']).name} ({mode})")
    if mode == "fail":
        for i in range(30):
            print(f"log line {i}", file=sys.stderr)
        print("Error: boom", file=sys.stderr)
        sys.exit(3)
    if mode == "sleep":
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
        if os.environ.get("FAKE_R_CHILD_PID"):
            Path(os.environ["FAKE_R_CHILD_PID"]).write_text(str(child.pid), encoding="utf-8")
        time.sleep(120)
    if mode != "no_output":
        for name, rel in p["expect"].items():
            (workdir / rel).write_text(f"{name} seed={p['seed']}\n", encoding="utf-8")
    kind = ["Mersenne-Twister", "Inversion", "Rejection"]
    session = {"r_version": f"R version {version} (fake)", "platform": "fake",
               "lib_paths": [p for p in os.environ.get("FAKE_R_LIBRARY", "").split(os.pathsep) if p],
               "rng_kind": kind, "rng_kind_after": kind, "set_seed_calls": [], "seed": p["seed"],
               "packages": {"jsonlite": "2.0.0"}, "seconds": 0.01}
    text = "not json {" if mode == "bad_session" else json.dumps(session)
    (workdir / "session.json").write_text(text, encoding="utf-8")
    sys.exit(0)
print(f"fake Rscript: unexpected arguments {sys.argv[1:]}", file=sys.stderr)
sys.exit(9)
'''


@dataclass
class FakeR:
    path: Path          # the launcher to pass as an Rscript
    script: Path        # the Python code behind every launcher
    calls_file: Path
    library: Path

    def launcher(self, where: Path, name: str = "Rscript") -> Path:
        """Another launcher of the same fake at ``where/name`` (``.cmd`` added on Windows)."""
        where.mkdir(parents=True, exist_ok=True)
        if os.name == "nt":
            exe = where / f"{name}.cmd"
            exe.write_text(f'@"{sys.executable}" "{self.script}" %*\r\n', encoding="utf-8")
        else:
            exe = where / name
            exe.write_text(f"#!{sys.executable}\n" + FAKE_R, encoding="utf-8")
            exe.chmod(0o755)
        return exe

    def calls(self) -> list[dict]:
        """One dict per start of the fake: script, argv, vanilla, path (first entries), r_libs_user, cwd."""
        if not self.calls_file.exists():
            return []
        return [json.loads(line) for line in self.calls_file.read_text(encoding="utf-8").splitlines()]


@pytest.fixture()
def fake_rscript(tmp_path, monkeypatch) -> FakeR:
    """A fake Rscript in ``tmp_path/fake_r`` that has jsonlite and SoupX; it counts its own calls."""
    home = tmp_path / "fake_r"
    home.mkdir()
    lib = home / "library"
    lib.mkdir()
    script = home / "fake_r.py"
    script.write_text(FAKE_R, encoding="utf-8")
    fake = FakeR(path=Path(), script=script, calls_file=home / "calls.txt", library=lib)
    fake.path = fake.launcher(home / "bin")
    monkeypatch.setenv("FAKE_R_CALLS", str(fake.calls_file))
    monkeypatch.setenv("FAKE_R_LIBRARY", str(lib))
    monkeypatch.setenv("FAKE_R_PACKAGES", json.dumps({"jsonlite": "2.0.0", "SoupX": "1.6.2"}))
    for name in ("FAKE_R_MODE", "FAKE_R_BUILT", "FAKE_R_VERSION", "FAKE_R_CHILD_PID"):
        monkeypatch.delenv(name, raising=False)
    from grnbench import rbridge
    rbridge.clear_memo()
    yield fake
    rbridge.clear_memo()


@pytest.fixture()
def synthetic_adata() -> ad.AnnData:
    """400 cells x 600 genes raw counts.

    - 12 MT- genes (2% of space) boosted 30x in 100 cells -> elevated but
      non-failing mito fraction.
    - log-normal library sizes (CV ~1) -> healthy gates.
    - three groups in obs['group']: A(200), B(195), C(5 tiny).
    """
    rng = np.random.default_rng(0)
    n_cells, n_genes = 400, 600
    lib = np.exp(rng.normal(7.5, 0.8, n_cells)).astype(np.float32)
    base = rng.negative_binomial(2, 0.2, size=(n_cells, n_genes)).astype(np.float32)
    X = np.round(base * lib[:, None])  # integer-valued counts (float32 exact)
    genes = [f"GENE{i:03d}" for i in range(n_genes)]
    genes[:12] = [f"MT-GENE{j}" for j in range(12)]
    X[:100, :12] *= 30.0  # still integer
    obs = pd.DataFrame({"group": ["A"] * 200 + ["B"] * 195 + ["C"] * 5})
    Xs = sparse.csr_matrix(X.astype(np.float32))
    adata = ad.AnnData(X=Xs.copy(), obs=obs, var=pd.DataFrame(index=genes))
    adata.layers["counts"] = Xs
    return adata


@pytest.fixture()
def extreme_mito_adata() -> ad.AnnData:
    """50 cells where mitochondrial genes dominate -> mito gate must FAIL."""
    rng = np.random.default_rng(7)
    X = rng.poisson(50.0, size=(50, 9)).astype(np.float32)  # MT- genes huge
    normal = rng.poisson(0.5, size=(50, 1)).astype(np.float32)
    X = np.hstack([X * 20.0, normal])
    return ad.AnnData(X=X,
                      var=pd.DataFrame(index=[f"MT-G{i}" for i in range(9)] + ["ACTB"]))


@pytest.fixture()
def cascade_adata() -> ad.AnnData:
    """200 cells x 20 genes with a planted linear cascade TF_A -> G_B -> G_C.

    Remaining genes are independent Poisson noise; expression is integer
    counts so the ingest contract holds.
    """
    rng = np.random.default_rng(11)
    n_cells = 200
    A = rng.poisson(8.0, n_cells).astype(np.float64)
    B = np.round(np.maximum(0.0, 0.9 * A - 2.0 + rng.normal(0, 1.5, n_cells)))
    C = np.round(np.maximum(0.0, 0.9 * B - 2.0 + rng.normal(0, 1.5, n_cells)))
    genes = {"TF_A": A, "G_B": B, "G_C": C}
    for i in range(17):
        genes[f"N{i:02d}"] = rng.poisson(5.0, n_cells).astype(np.float64)
    cols = list(genes)
    X = np.column_stack([genes[g] for g in cols]).astype(np.float32)
    return ad.AnnData(X=X.astype(np.float32),
                      var=pd.DataFrame(index=cols))
