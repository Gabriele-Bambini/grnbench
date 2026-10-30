"""The HPC job scripts pass what a registered experiment needs, checked without Slurm or conda:
conda and python are shell stubs that record their arguments."""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _bash():
    bash = shutil.which("bash") or "C:/Program Files/Git/bin/bash.exe"
    if not Path(bash).exists():
        pytest.skip("Bash unavailable")
    return bash


def _stubs(tmp_path: Path) -> Path:
    log = tmp_path / "calls.log"
    profile = tmp_path / "conda.sh"
    profile.write_text(
        'conda() { if [ "$1" = "env" ]; then echo test-env; fi; return 0; }\n'
        "grnbench() { return 0; }\n"
        # ${PYTHONPATH:-}: setup_env.sh runs under `set -u`, and PYTHONPATH is usually unset
        f'python() {{ echo "python $* PYTHONPATH=${{PYTHONPATH:-}}" >> "{log.as_posix()}"; return 0; }}\n',
        encoding="utf-8")
    return profile


def test_sbatch_passes_models_artifacts_and_records_the_commit(tmp_path):
    bash = _bash()
    repo, out = tmp_path / "repo", tmp_path / "out"
    repo.mkdir()
    (repo / "GIT_COMMIT").write_text("abc1234\n")
    env = {**os.environ, "CONDA_SH": _stubs(tmp_path).as_posix(), "REPO_DIR": repo.as_posix(),
           "OUT": out.as_posix(), "SCRIPT": "benchmarks/real_bench.py", "REFS": "/refs",
           "MODELS": "regdiffusion+regdiffusion_motif+grnboost2", "KEEP_ARTIFACTS": "1",
           "CHUNK": "100000", "SLURM_CPUS_PER_TASK": "8", "TMPDIR": tmp_path.as_posix()}
    r = subprocess.run([bash, str(ROOT / "hpc" / "engine_bench.sbatch")], env=env,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    call = [line for line in (tmp_path / "calls.log").read_text().splitlines() if " run " in line][0]
    assert "--models regdiffusion,regdiffusion_motif,grnboost2" in call
    assert "--keep-artifacts" in call and "--refs /refs" in call and "--chunk 100000" in call
    assert "--max-attempts 2" in call
    assert (out / "CODE_COMMIT").read_text().strip() == "abc1234"
    # the job imports the code of the checkout CODE_COMMIT names, not the env's editable install
    assert f"PYTHONPATH={repo.as_posix()}/src" in call


def test_two_jobs_sharing_tmpdir_get_their_own_work_folders(tmp_path):
    """Without --array both jobs had work folder gb_engine_0_0, and the first to end removed it under
    the other (the 2026-10-23 pool pilot: FileNotFoundError writing a manifest.json, #139)."""
    bash = _bash()
    repo = tmp_path / "repo"
    repo.mkdir()
    (repo / "GIT_COMMIT").write_text("abc1234\n")
    roots = []
    for job in ("715093", "715094"):
        env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_ARRAY")}
        env.update({"CONDA_SH": _stubs(tmp_path).as_posix(), "REPO_DIR": repo.as_posix(),
                    "OUT": (tmp_path / f"out{job}").as_posix(), "SCRIPT": "benchmarks/real_bench.py",
                    "CHUNK": "100000", "SLURM_CPUS_PER_TASK": "2", "TMPDIR": tmp_path.as_posix(),
                    "SLURM_JOB_ID": job})
        r = subprocess.run([bash, str(ROOT / "hpc" / "engine_bench.sbatch")], env=env,
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
        call = [line for line in (tmp_path / "calls.log").read_text().splitlines() if " run " in line][-1]
        roots.append(call.split("--workroot ")[1].split()[0])
    assert roots[0] != roots[1] and all(job in root for job, root in zip(("715093", "715094"), roots))


def test_sbatch_extra_arguments_are_safe_empty_and_keep_spaces(tmp_path):
    """An empty EXTRA under set -u stops Bash 3.2 (macOS) with 'EXTRA[@]: unbound variable', which a
    newer Bash does not show, so the expansion must be the guarded form; a path with a space stays
    one argument."""
    text = (ROOT / "hpc" / "engine_bench.sbatch").read_text(encoding="utf-8")
    assert text.count('"${EXTRA[@]}"') == text.count('${EXTRA[@]+"${EXTRA[@]}"}') == 1
    bash = _bash()
    repo = tmp_path / "repo"
    repo.mkdir()
    log = tmp_path / "args.log"
    profile = tmp_path / "conda.sh"
    profile.write_text(
        'conda() { if [ "$1" = "env" ]; then echo test-env; fi; return 0; }\n'
        f'python() {{ for a in "$@"; do printf "[%s]" "$a"; done >> "{log.as_posix()}"; echo >> "{log.as_posix()}"; }}\n',
        encoding="utf-8")
    env = {**os.environ, "CONDA_SH": profile.as_posix(), "REPO_DIR": repo.as_posix(),
           "OUT": (tmp_path / "out").as_posix(), "SCRIPT": "benchmarks/real_bench.py",
           "REFS": "/data/my refs", "CHUNK": "100000", "SLURM_CPUS_PER_TASK": "2",
           "TMPDIR": tmp_path.as_posix()}
    r = subprocess.run([bash, str(ROOT / "hpc" / "engine_bench.sbatch")], env=env,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    call = [line for line in log.read_text().splitlines() if "[run]" in line][0]
    assert "[--refs][/data/my refs]" in call


@pytest.mark.parametrize("script", ["grnbench.sbatch", "beeline_backtest.sbatch", "engine_bench.sbatch"])
def test_spooled_sbatch_resolves_the_submit_dir_checkout(tmp_path, script):
    """sbatch runs a copy of the script from the node's spool directory, so a REPO_DIR derived
    from BASH_SOURCE points at /var/spool/slurmd, not at the checkout `hpc.sh submit` cd'ed into."""
    bash = _bash()
    line = next(x for x in (ROOT / "hpc" / script).read_text(encoding="utf-8").splitlines()
                if x.startswith("REPO_DIR="))
    spool = tmp_path / "var" / "spool" / "slurmd" / "job4242"
    spool.mkdir(parents=True)
    (spool / "slurm_script").write_text(line + '\necho "$REPO_DIR"\n', encoding="utf-8")
    checkout = (tmp_path / "checkout").as_posix()
    env = {k: v for k, v in os.environ.items() if k != "REPO_DIR"}
    env["SLURM_SUBMIT_DIR"] = checkout
    r = subprocess.run([bash, str(spool / "slurm_script")], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == checkout


@pytest.mark.parametrize("script", ["grnbench.sbatch", "beeline_backtest.sbatch"])
@pytest.mark.parametrize("is_checkout", [False, True])
def test_sbatch_refuses_a_repo_dir_without_pyproject(tmp_path, script, is_checkout):
    """A submit from the wrong directory used to fail late, after conda and data checks."""
    bash = _bash()
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    if is_checkout:
        (checkout / "pyproject.toml").write_text("[project]\n", encoding="utf-8")
    env = {**os.environ, "REPO_DIR": checkout.as_posix(), "CONDA_SH": _stubs(tmp_path).as_posix(),
           "BEELINE": (tmp_path / "beeline").as_posix()}
    env.pop("DATA", None)
    r = subprocess.run([bash, str(ROOT / "hpc" / script)], env=env, capture_output=True, text=True)
    assert r.returncode != 0  # with a checkout it stops later, at the missing DATA/BEELINE input
    assert ("no pyproject.toml" in r.stderr) is not is_checkout, r.stderr


def test_setup_env_with_torch_runs_the_lock_header_line(tmp_path):
    bash = _bash()
    env = {**os.environ, "CONDA_SH": _stubs(tmp_path).as_posix(), "WITH_TORCH": "1"}
    r = subprocess.run([bash, str(ROOT / "hpc" / "setup_env.sh"), "test-env", str(ROOT)], env=env,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    calls = (tmp_path / "calls.log").read_text()
    lock_header = (ROOT / "requirements" / "lock-py311.txt").read_text(encoding="utf-8")
    pinned = next(line for line in lock_header.splitlines() if "pip install regdiffusion" in line)
    assert pinned.lstrip("# ").strip() in calls
    assert "import torch, regdiffusion" in calls


def test_setup_env_torch_line_survives_a_crlf_lock(tmp_path):
    """A Windows checkout (text=auto) syncs the lock with CRLF; the install line read from its
    header must not end in '\r', or pip gets 'https://pypi.org/simple\r' and finds nothing
    (what happened on the cluster on 2026-09-27)."""
    bash = _bash()
    repo = tmp_path / "repo"
    (repo / "requirements").mkdir(parents=True)
    header = [line for line in (ROOT / "requirements" / "lock-py311.txt").read_text(encoding="utf-8").splitlines()
              if line.startswith("#")]
    (repo / "requirements" / "lock-py311.txt").write_bytes(("\r\n".join(header) + "\r\n").encode())
    env = {**os.environ, "CONDA_SH": _stubs(tmp_path).as_posix(), "WITH_TORCH": "1"}
    r = subprocess.run([bash, str(ROOT / "hpc" / "setup_env.sh"), "test-env", str(repo)], env=env,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    calls = (tmp_path / "calls.log").read_bytes()
    assert b"\r" not in calls
    assert b"--extra-index-url https://pypi.org/simple PYTHONPATH=" in calls


def test_identifier_calibration_stages_only_the_registered_snapshot(tmp_path):
    bash = _bash()
    repo, data, out, scratch = (tmp_path / name for name in ("repo", "data", "out", "scratch"))
    (repo / "experiments").mkdir(parents=True)
    (data / "cohort").mkdir(parents=True)
    scratch.mkdir()
    (repo / "GIT_COMMIT").write_text("abc1234\n", encoding="utf-8")
    (data / "cohort" / "data.h5ad").write_bytes(b"registered")
    (data / "unrelated.h5ad").write_bytes(b"not in the snapshot")
    snapshot = repo / "experiments" / "EXP-TEST.snapshot.tsv"
    snapshot.write_text(
        "path\tbytes\tsha256\tfile_date\tsource\tversion\n"
        "cohort/data.h5ad\t10\tdigest\t2026-10-10\ttest\tv1\n",
        encoding="utf-8",
    )
    env = {
        **os.environ,
        "CONDA_SH": _stubs(tmp_path).as_posix(),
        "REPO_DIR": repo.as_posix(),
        "DATA_ROOT": data.as_posix(),
        "OUT": out.as_posix(),
        "EXP_ID": "EXP-TEST",
        "SNAPSHOT": snapshot.as_posix(),
        "SLURM_TMPDIR": scratch.as_posix(),
        "SLURM_JOB_ID": "123",
    }

    r = subprocess.run(
        [bash, str(ROOT / "hpc" / "identifier_policy_calibration.sbatch")],
        env=env,
        capture_output=True,
        text=True,
    )

    assert r.returncode == 0, r.stdout + r.stderr
    staged = scratch / "grnbench-EXP-TEST-123-inputs"
    staged_data = staged / "cohort" / "data.h5ad"
    assert staged_data.read_bytes() == b"registered"
    if os.name != "nt":
        assert staged_data.is_symlink()
    assert not (staged / "unrelated.h5ad").exists()
    assert (out / "input-snapshot.tsv").read_text(encoding="utf-8") == snapshot.read_text(encoding="utf-8")
    calls = (tmp_path / "calls.log").read_text(encoding="utf-8")
    assert f"data_manifest.py {staged.as_posix()} --check {snapshot.as_posix()} --require-all" in calls
    assert f"identifier_policy_calibration.py {staged.as_posix()}" in calls


def test_setup_pyscenic_env_uses_conda_forge_and_installs_requirements(tmp_path):
    bash = _bash()
    log = tmp_path / "calls.log"
    profile = tmp_path / "conda.sh"
    profile.write_text(
        f'conda() {{ echo "conda $*" >> "{log.as_posix()}"; if [ "$1" = "env" ] && [ "$2" = "list" ]; then echo "base"; fi; return 0; }}\n'
        f'python() {{ echo "python $*" >> "{log.as_posix()}"; return 0; }}\n'
        'pyscenic() { return 0; }\n',
        encoding="utf-8",
    )
    env = {**os.environ, "CONDA_SH": profile.as_posix()}
    r = subprocess.run([bash, str(ROOT / "hpc" / "setup_pyscenic_env.sh"), "test-pyscenic-env", ROOT.as_posix()],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    calls = log.read_text(encoding="utf-8")
    assert "conda create -y -c conda-forge -n test-pyscenic-env python=3.10 pip=23.3.2" in calls
    req = (ROOT / "requirements" / "pyscenic-0.12.1-py310.txt").as_posix()
    assert f"python -m pip install -r {req}" in calls


@pytest.mark.parametrize("script", ["grnbench.sbatch", "beeline_backtest.sbatch", "engine_bench.sbatch"])
def test_sbatch_imports_the_checkout_after_activating_the_env(script):
    lines = (ROOT / "hpc" / script).read_text(encoding="utf-8").splitlines()
    export = next(i for i, x in enumerate(lines)
                  if x == 'export PYTHONPATH="$REPO_DIR/src${PYTHONPATH:+:$PYTHONPATH}"')
    assert export > next(i for i, x in enumerate(lines) if "conda activate" in x)
