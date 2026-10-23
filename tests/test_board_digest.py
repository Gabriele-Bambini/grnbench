import os
import subprocess
import sys
import textwrap

from grnbench.board import claims_table, digest, main

FAKE_GH = '''
import json, sys
a = sys.argv[1:]
if a[:2] == ["pr", "list"] and "merged" in a:
    print(json.dumps(MERGED))
elif a[:2] == ["pr", "list"]:
    print(json.dumps(OPEN))
elif a[:2] == ["issue", "list"]:
    print(json.dumps(ISSUES))
else:
    sys.exit("unexpected gh call: %r" % a)
'''


def fake_gh(tmp_path, merged=(), opened=(), issues=()):
    f = tmp_path / "fake_gh.py"
    f.write_text(f"MERGED={list(merged)!r}\nOPEN={list(opened)!r}\nISSUES={list(issues)!r}\n" + FAKE_GH)
    return [sys.executable, str(f)]


def commit(repo, msg, date):
    env = {**os.environ, "GIT_AUTHOR_DATE": date, "GIT_COMMITTER_DATE": date}
    subprocess.run(["git", "-C", str(repo), "add", "-A"], check=True)
    subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", msg],
                   check=True, env=env)


CLAIMS = "| # | claim | experiments | files | status | reproduced by |\n|---|---|---|---|---|---|\n" \
         "| C1 | text | EXP-20260101-AA-1 (confirmatory) | `x` | {} | - |\n"
DECISIONS = "### D-001 | 2026-01-01 | a | ACTIVE | First\n"


def make_repo(tmp_path):
    repo = tmp_path / "repo"
    for d in ("paper/tables", "experiments", ".projectos"):
        (repo / d).mkdir(parents=True)
    subprocess.run(["git", "-C", str(repo), "init", "-q"], check=True)
    (repo / "paper/CLAIMS.md").write_text(CLAIMS.format("draft"))
    (repo / ".projectos/DECISIONS.md").write_text(DECISIONS)
    commit(repo, "old", "2026-01-01T12:00:00")
    return repo


def test_empty_week(tmp_path):
    repo = make_repo(tmp_path)
    out = digest("2026-02-01", gh=fake_gh(tmp_path), root=repo)
    for title in ("PRs merged", "PRs open", "Experiment records added or changed", "Claims whose status changed",
                  "Decisions added", "Open issues labelled running/experiment"):
        assert out.split(f"## {title}")[1].split("\n\n")[0].strip() == "- none"


def test_merged_open_and_issues(tmp_path):
    repo = make_repo(tmp_path)
    gh = fake_gh(tmp_path,
                 merged=[{"number": 5, "title": "fix x", "author": {"login": "lm"}}],
                 opened=[{"number": 6, "title": "wip", "author": {"login": "ld"}, "isDraft": False,
                          "reviewDecision": "APPROVED"}],
                 issues=[{"number": 9, "title": "EXP run", "assignees": [{"login": "lm"}]}])
    out = digest("2026-02-01", gh=gh, root=repo)
    assert "- #5 fix x (lm)" in out and "- #6 wip (ld, approved)" in out and "- #9 EXP run (lm)" in out


def test_claim_change_decision_and_experiment(tmp_path):
    repo = make_repo(tmp_path)
    (repo / "paper/CLAIMS.md").write_text(CLAIMS.format("confirmed"))
    (repo / ".projectos/DECISIONS.md").write_text(DECISIONS + "### D-002 | 2026-02-10 | a | ACTIVE | Second\n")
    (repo / "experiments/EXP-20260210-AA-1.md").write_text(
        "# EXP\n- type: **confirmatory** (x)\n- status: done (2026-02-10)\n")
    commit(repo, "new", "2026-02-10T12:00:00")
    out = digest("2026-02-01", gh=fake_gh(tmp_path), root=repo)
    assert "- C1: draft -> confirmed" in out
    assert "- D-002: Second" in out and "D-001" not in out
    assert "- EXP-20260210-AA-1 (confirmatory): status done (2026-02-10)" in out


def test_claims_table_joins_experiments_and_outputs(tmp_path):
    repo = make_repo(tmp_path)
    (repo / "paper/CLAIMS.md").write_text(CLAIMS.replace("| - |\n", "| rerun, issue #15 |\n").format("reproduced"))
    (repo / "paper/tables/T_a.md").write_text("from EXP-20260101-AA-1")
    (repo / "paper/tables/T_b.md").write_text("unrelated")
    out = claims_table(repo)
    assert "| C1 | reproduced | rerun, issue #15 | EXP-20260101-AA-1 | #15 | paper/tables/T_a.md |" in out


def test_post_only_when_asked(tmp_path, monkeypatch, capsys):
    repo = make_repo(tmp_path)
    log = tmp_path / "calls.txt"
    f = tmp_path / "gh_rec.py"
    f.write_text(textwrap.dedent(f'''
        import sys
        open({str(log)!r}, "a").write(" ".join(sys.argv[1:3]) + "\\n")
        if sys.argv[1] == "issue" and sys.argv[2] == "comment":
            open({str(log)!r}, "a").write("BODY:" + sys.stdin.read()[:20] + "\\n")
        else:
            print("[]")
    '''))
    monkeypatch.setenv("GRNBENCH_GH", f'"{sys.executable}" "{f}"')
    monkeypatch.chdir(repo)
    main(["--digest", "--since", "2026-02-01"])
    assert "comment" not in log.read_text()
    main(["--digest", "--since", "2026-02-01", "--post", "12"])
    assert "issue comment" in log.read_text() and "BODY:# grnbench digest" in log.read_text()


def test_digest_reads_origin_main_else_main_and_says_which(tmp_path):
    """A checked-out feature branch must not leak into the digest."""
    repo = make_repo(tmp_path)
    git = lambda *a: subprocess.run(["git", "-C", str(repo), *a], check=True, capture_output=True, text=True)
    git("branch", "-M", "main")
    git("checkout", "-q", "-b", "feat")
    (repo / ".projectos/DECISIONS.md").write_text(DECISIONS + "### D-002 | 2026-02-10 | a | ACTIVE | Second\n")
    commit(repo, "on the branch only", "2026-02-10T12:00:00")
    out = digest("2026-02-01", gh=fake_gh(tmp_path), root=repo)
    assert "files from main)" in out and "D-002" not in out
    git("update-ref", "refs/remotes/origin/main", "feat")
    out = digest("2026-02-01", gh=fake_gh(tmp_path), root=repo)
    assert "files from origin/main)" in out and "- D-002: Second" in out
