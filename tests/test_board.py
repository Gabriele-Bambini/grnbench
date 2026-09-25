import sys
import textwrap

from grnbench.board import board


def test_board_lists_running_issues_prs_and_state(tmp_path):
    fake = tmp_path / "fake_gh.py"
    fake.write_text(textwrap.dedent('''
        import json, sys
        if sys.argv[1] == "issue":
            print(json.dumps([{"number": 7, "title": "EXP-20260926-GB-1 synthetic v3",
                               "assignees": [{"login": "someone"}], "updatedAt": "2026-09-26T10:00:00Z"}]))
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
    assert "- job 123 on HPC" in out and "1. do the thing" in out
    assert "- n" not in out.split("STATE.md: Next steps")[1]


def test_board_survives_missing_gh(tmp_path):
    out = board(gh=[str(tmp_path / "no-such-gh")], root=tmp_path)
    assert "gh not usable" in out
