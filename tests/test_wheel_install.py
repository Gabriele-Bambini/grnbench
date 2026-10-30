import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def installed_wheel(tmp_path_factory):
    """Build a non-editable wheel and install it into an isolated target folder."""
    dist_dir = tmp_path_factory.mktemp("dist")
    target_dir = tmp_path_factory.mktemp("target")

    pip_env = dict(os.environ)
    pip_env["PIP_NO_INDEX"] = "1"
    pip_env["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"

    build_cmd = [
        sys.executable,
        "-m",
        "pip",
        "wheel",
        "--no-deps",
        "--no-build-isolation",
        "-w",
        str(dist_dir),
        str(ROOT),
    ]
    build_res = subprocess.run(build_cmd, env=pip_env, capture_output=True, text=True)
    assert build_res.returncode == 0, f"pip wheel failed:\n{build_res.stderr}\n{build_res.stdout}"

    wheels = list(dist_dir.glob("*.whl"))
    assert len(wheels) == 1, f"Expected 1 wheel, found {wheels}"
    wheel_path = wheels[0]

    install_cmd = [
        sys.executable,
        "-m",
        "pip",
        "install",
        "--target",
        str(target_dir),
        "--no-deps",
        str(wheel_path),
    ]
    install_res = subprocess.run(install_cmd, env=pip_env, capture_output=True, text=True)
    assert install_res.returncode == 0, f"pip install --target failed:\n{install_res.stderr}\n{install_res.stdout}"

    return target_dir


def _run_isolated(
    target_dir: Path,
    args: list[str],
    cwd: Path | None = None,
    extra_env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    if extra_env:
        env.update(extra_env)
    # Prepend target_dir to PYTHONPATH so grnbench is loaded exclusively from the installed wheel,
    # and strip any checkout paths if present.
    checkout_dirs = {ROOT.resolve(), (ROOT / "src").resolve(), (ROOT / "benchmarks").resolve()}
    existing_pythonpath = [
        p for p in env.get("PYTHONPATH", "").split(os.pathsep)
        if p and Path(p).resolve() not in checkout_dirs
    ]
    env["PYTHONPATH"] = os.pathsep.join([str(target_dir), *existing_pythonpath])
    cmd = [sys.executable, "-m", "grnbench", *args]
    return subprocess.run(cmd, cwd=str(cwd) if cwd is not None else None, env=env, capture_output=True, text=True)


def test_wheel_import_provenance(installed_wheel):
    """Ensure grnbench imports from the installed target directory, not the editable checkout."""
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{installed_wheel}{os.pathsep}{env.get('PYTHONPATH', '')}"
    cmd = [sys.executable, "-c", "import grnbench; print(grnbench.__file__)"]
    res = subprocess.run(cmd, env=env, capture_output=True, text=True)
    assert res.returncode == 0
    imported_path = Path(res.stdout.strip()).resolve()
    assert str(installed_wheel.resolve()) in str(imported_path), f"{imported_path} not in {installed_wheel}"


def test_wheel_board_outside_checkout_fails_cleanly(installed_wheel, tmp_path):
    """grnbench board outside a repository checkout returns exit code 2 and a clean error."""
    res = _run_isolated(installed_wheel, ["board"], cwd=tmp_path)
    assert res.returncode == 2
    assert "grnbench board: repository root not found" in res.stderr
    assert res.stdout == ""
    assert "Traceback" not in res.stderr


@pytest.mark.parametrize("flag", ["--digest", "--claims"])
def test_wheel_board_modes_outside_checkout_fail_cleanly(installed_wheel, tmp_path, flag):
    """Every repository-backed board mode has the same controlled failure."""
    res = _run_isolated(installed_wheel, ["board", flag], cwd=tmp_path)
    assert res.returncode == 2
    assert "grnbench board: repository root not found" in res.stderr
    assert res.stdout == ""
    assert "Traceback" not in res.stderr


def test_wheel_board_in_checkout_success(installed_wheel, tmp_path):
    """grnbench board from a nested directory inside a checkout succeeds with code 0."""
    fake_repo = tmp_path / "fake_repo"
    fake_repo.mkdir()
    (fake_repo / ".projectos").mkdir()
    (fake_repo / ".projectos" / "STATE.md").write_text(
        "# State - fake repo\n\n## In progress\n- fake wheel task in progress\n\n## Next steps\n- fake next step\n\n## Blockers\n- fake blocker\n",
        encoding="utf-8",
    )
    (fake_repo / "AGENTS.md").write_text("# Agents for fake repo\n", encoding="utf-8")
    (fake_repo / "pyproject.toml").write_text(
        '[project]\nname = "grnbench"\n', encoding="utf-8"
    )

    nested = fake_repo / "nested" / "subfolder"
    nested.mkdir(parents=True)

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_program = fake_bin / "fake_gh.py"
    fake_program.write_text(f"""
import json, sys
from pathlib import Path

# Verify that gh is executed with cwd at the fake repository root
cwd = Path.cwd().resolve()
expected = Path({repr(str(fake_repo.resolve()))})
if cwd != expected:
    sys.stderr.write(f"fake gh cwd mismatch: {{cwd}} != {{expected}}\\n")
    sys.exit(1)

if len(sys.argv) > 1 and sys.argv[1] == "issue":
    if "--label" in sys.argv:
        label_idx = sys.argv.index("--label") + 1
        label = sys.argv[label_idx] if label_idx < len(sys.argv) else ""
        if label == "running":
            print(json.dumps([{{
                "number": 101,
                "title": "fake running issue",
                "assignees": [{{"login": "fake_user"}}],
                "updatedAt": "2026-10-16T12:00:00Z",
                "url": "https://example.com/101",
            }}]))
            sys.exit(0)
    print("[]")
    sys.exit(0)
elif len(sys.argv) > 1 and sys.argv[1] == "pr":
    print(json.dumps([{{
        "number": 202,
        "title": "fake open PR",
        "author": {{"login": "fake_author"}},
        "headRefName": "fake/branch",
        "isDraft": False,
    }}]))
    sys.exit(0)

print("[]")
""")
    if os.name == "nt":
        fake_gh = fake_bin / "gh.cmd"
        fake_gh.write_text(
            f'@echo off\r\n"{sys.executable}" "{fake_program}" %*\r\n',
            encoding="utf-8",
        )
    else:
        fake_gh = fake_bin / "gh"
        fake_gh.write_text(
            f"#!{sys.executable}\n" + fake_program.read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        fake_gh.chmod(0o755)

    extra_env = {
        "PATH": f"{fake_bin}{os.pathsep}{os.environ.get('PATH', '')}",
        "GRNBENCH_GH": str(fake_gh),
    }

    res = _run_isolated(installed_wheel, ["board"], cwd=nested, extra_env=extra_env)
    assert res.returncode == 0, f"board failed:\n{res.stderr}\n{res.stdout}"
    assert res.stderr == ""
    assert "#101 fake running issue" in res.stdout
    assert "#202 fake open PR" in res.stdout
    assert "- fake wheel task in progress" in res.stdout
    assert "- fake next step" in res.stdout
    assert "- fake blocker" in res.stdout
    assert "Packaging and wheel resource safety" not in res.stdout
    assert "Traceback" not in res.stderr


def test_wheel_direct_run_directory_report(installed_wheel, tmp_path):
    """grnbench report on an explicit run directory succeeds outside a checkout without real_bench."""
    run_dir = tmp_path / "standalone_run"
    rows_dir = run_dir / "rows"
    rows_dir.mkdir(parents=True)
    row = {
        "task": "cxg__d0",
        "error": None,
        "m_pearson_signed_auto_auprc_ratio": 1.25,
        "m_pearson_motif_auto_auprc_ratio": 1.45,
        "m_motif_prior_auto_auprc_ratio": 1.35,
        "auto_doublets": "scrublet",
    }
    (rows_dir / "d0.json").write_text(json.dumps(row))

    out_dir = tmp_path / "latex_out"
    res = _run_isolated(installed_wheel, ["report", str(run_dir), "--no-compile", "--out", str(out_dir)], cwd=tmp_path)
    assert res.returncode == 0, f"Report failed:\n{res.stderr}\n{res.stdout}"
    assert (out_dir / "report.tex").exists()
    assert (out_dir / "values.tex").exists()
    assert (out_dir / "tables" / "configuration.tex").exists()
    assert "Traceback" not in res.stderr


def test_wheel_exp_id_outside_checkout_fails_cleanly(installed_wheel, tmp_path):
    """grnbench report with an EXP ID outside a checkout fails cleanly with code 2."""
    res = _run_isolated(installed_wheel, ["report", "EXP-20261001-LD-1", "--no-compile"], cwd=tmp_path)
    assert res.returncode == 2
    assert "not a folder and grnbench repository root not found" in res.stderr
    assert "Traceback" not in res.stderr


def test_wheel_paper_outside_checkout_fails_cleanly(installed_wheel, tmp_path):
    """grnbench report --paper outside a checkout fails cleanly with code 2."""
    res = _run_isolated(installed_wheel, ["report", "--paper", "--no-compile"], cwd=tmp_path)
    assert res.returncode == 2
    assert "requires a grnbench repository checkout" in res.stderr
    assert "Traceback" not in res.stderr


def test_wheel_paper_in_checkout_success(installed_wheel, tmp_path):
    """grnbench report --paper --no-compile in a checkout succeeds and writes paper LaTeX output."""
    fake_repo = tmp_path / "fake_paper_repo"
    fake_repo.mkdir()
    (fake_repo / ".projectos").mkdir()
    (fake_repo / ".projectos" / "STATE.md").write_text("# State\n", encoding="utf-8")
    (fake_repo / "AGENTS.md").write_text("# Agents\n", encoding="utf-8")
    (fake_repo / "pyproject.toml").write_text(
        '[project]\nname = "grnbench"\n', encoding="utf-8"
    )

    (fake_repo / "paper").mkdir()
    claims_text = (
        "| # | claim | experiments | files | status | reproduced by |\n"
        "|---|---|---|---|---|---|\n"
        "| C1 | fake confirmatory claim | EXP-20261016-TST-1 | results | confirmed | - |\n"
    )
    (fake_repo / "paper" / "CLAIMS.md").write_text(claims_text, encoding="utf-8")

    exp = "EXP-20261016-TST-1"
    (fake_repo / "experiments").mkdir()
    record_text = (
        f"# {exp} Fake Confirmatory Experiment\n\n"
        "- type: confirmatory\n"
        "- primary comparison: pearson_motif > pearson_signed | auprc_ratio | auto\n"
        f"- report input: `benchmarks/results/{exp}`\n"
    )
    (fake_repo / "experiments" / f"{exp}.md").write_text(record_text, encoding="utf-8")

    run_dir = fake_repo / "benchmarks" / "results" / exp
    rows_dir = run_dir / "rows"
    rows_dir.mkdir(parents=True)
    for i in range(4):
        row = {
            "task": f"cxg__d{i}",
            "error": None,
            "m_pearson_signed_auto_auprc_ratio": 1.20 + i * 0.02,
            "m_pearson_motif_auto_auprc_ratio": 1.40 + i * 0.02,
            "m_motif_prior_auto_auprc_ratio": 1.30 + i * 0.02,
            "auto_doublets": "scrublet",
        }
        (rows_dir / f"d{i}.json").write_text(json.dumps(row))

    res = _run_isolated(installed_wheel, ["report", "--paper", "--no-compile"], cwd=fake_repo)
    assert res.returncode == 0, f"report --paper failed:\n{res.stderr}\n{res.stdout}"
    assert "Traceback" not in res.stderr

    paper_tex = fake_repo / "paper" / "build" / "paper.tex"
    assert paper_tex.exists()
    paper_content = paper_tex.read_text(encoding="utf-8")
    assert "\\begin{document}" in paper_content
    assert f"\\section{{{exp}}}" in paper_content
    assert (fake_repo / "paper" / "build" / exp / "values.tex").exists()
    assert (fake_repo / "paper" / "build" / "references.bib").exists()
    assert "EXP-20260925-LM-1" not in paper_content
