import sys
import textwrap

import pytest

from grnbench.board import _split_command, board


def test_split_command_preserves_windows_paths_and_removes_outer_quotes():
    assert _split_command(r"C:\Users\x\bin\gh.cmd", posix=False) == [r"C:\Users\x\bin\gh.cmd"]
    assert _split_command(
        r'"C:\Program Files\GitHub CLI\gh.exe" --repo owner/project',
        posix=False,
    ) == [r"C:\Program Files\GitHub CLI\gh.exe", "--repo", "owner/project"]


def test_board_lists_running_issues_prs_and_state(tmp_path):
    fake = tmp_path / "fake_gh.py"
    fake.write_text(textwrap.dedent('''
        import json, sys
        run = {"number": 7, "title": "EXP-20260926-GB-1 synthetic v3", "assignees": [{"login": "someone"}],
               "updatedAt": "2026-09-26T10:00:00Z", "labels": [{"name": "experiment"}, {"name": "running"}]}
        claim = {"number": 9, "title": "EXP-20260927-LM-3 HPC motif", "assignees": [{"login": "other"}],
                 "labels": [{"name": "experiment"}]}
        if sys.argv[1] == "issue":
            label = sys.argv[sys.argv.index("--label") + 1]
            print(json.dumps([run] if label == "running" else [run, claim]))
        else:
            print(json.dumps([{"number": 3, "title": "scrublet import fix", "author": {"login": "a"},
                               "headRefName": "lm/scrublet-import-fix", "isDraft": False}]))
    '''))
    (tmp_path / ".projectos").mkdir()
    (tmp_path / ".projectos" / "STATE.md").write_text(
        "# State\n## In progress\n- job 123 on HPC\n## Next steps (x)\n1. do the thing\n## Notes\n- n\n")
    out = board(gh=[sys.executable, str(fake)], root=tmp_path)
    assert "#7 EXP-20260926-GB-1 synthetic v3 (someone" in out
    assert "#3 scrublet import fix" in out
    claimed = out.split("## Claimed experiments")[1].split("## Open pull requests")[0]
    assert "#9 EXP-20260927-LM-3 HPC motif (other)" in claimed and "#7" not in claimed
    assert "- job 123 on HPC" in out and "1. do the thing" in out
    assert "- n" not in out.split("STATE.md: Next steps")[1]


def test_board_survives_missing_gh(tmp_path):
    out = board(gh=[str(tmp_path / "no-such-gh")], root=tmp_path)
    assert "gh not usable" in out


def test_project_root_discovery(tmp_path):
    from grnbench.project import find_project_root

    assert find_project_root(tmp_path) is None

    # Bare .git directory alone is intentionally not accepted
    (tmp_path / ".git").mkdir()
    assert find_project_root(tmp_path) is None

    # Generic project markers are not enough: another repository may use the
    # same team files without being grnbench.
    (tmp_path / ".projectos").mkdir()
    (tmp_path / ".projectos" / "STATE.md").write_text("# State\n")
    (tmp_path / "AGENTS.md").write_text("# Agents\n")
    (tmp_path / "pyproject.toml").write_text('[project]\nname = "another-project"\n')
    assert find_project_root(tmp_path) is None

    (tmp_path / "pyproject.toml").write_text('[project]\nname = "grnbench"\n')
    assert find_project_root(tmp_path) == tmp_path
    nested = tmp_path / "sub" / "dir"
    nested.mkdir(parents=True)
    assert find_project_root(nested) == tmp_path


@pytest.mark.parametrize("args", [[], ["--digest"], ["--claims"]])
def test_board_outside_checkout_returns_code_2(tmp_path, monkeypatch, capsys, args):
    from grnbench.board import main

    monkeypatch.chdir(tmp_path)
    code = main(args)
    assert code == 2
    captured = capsys.readouterr()
    assert "grnbench board: repository root not found" in captured.err
    assert captured.out == ""


def test_board_from_subdirectory_finds_root_and_runs_gh_in_root(tmp_path, monkeypatch):
    import json
    from grnbench.board import board

    (tmp_path / ".projectos").mkdir()
    (tmp_path / ".projectos" / "STATE.md").write_text("# State\n## In progress\n- test item\n")
    (tmp_path / "AGENTS.md").write_text("# Agents\n")
    (tmp_path / "pyproject.toml").write_text('[project]\nname = "grnbench"\n')

    fake = tmp_path / "fake_gh.py"
    fake.write_text(textwrap.dedent(f'''
        import json, os, sys
        from pathlib import Path
        cwd = Path.cwd().resolve()
        expected = Path({repr(str(tmp_path.resolve()))})
        if cwd != expected:
            sys.stderr.write(f"wrong cwd: {{cwd}} != {{expected}}\\n")
            sys.exit(1)
        print("[]")
    '''))

    sub = tmp_path / "nested" / "subfolder"
    sub.mkdir(parents=True)
    monkeypatch.chdir(sub)

    out = board(gh=[sys.executable, str(fake)])
    assert "## In flight" in out
    assert "- test item" in out
