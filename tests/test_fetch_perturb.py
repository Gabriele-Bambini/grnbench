"""benchmarks/fetch_perturb.py: pinned md5, parallel byte ranges, skip when present, resume, refuse a bad file."""
import hashlib
import importlib.util
import io
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("fetch_perturb", ROOT / "benchmarks" / "fetch_perturb.py")
fp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fp)

BODY = bytes(range(256)) * 97          # 24,832 bytes, not a multiple of the range count


class _Resp(io.BytesIO):
    status = 206

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


@pytest.fixture
def server(monkeypatch):
    """Serve BODY under a fake pinned key, honouring 'Range: bytes=a-b' like Zenodo does."""
    calls = []
    monkeypatch.setitem(fp.FILES, "fake", ("fake.h5ad", len(BODY), hashlib.md5(BODY).hexdigest()))

    def urlopen(req, timeout):
        a, b = map(int, re.fullmatch(r"bytes=(\d+)-(\d+)", req.headers["Range"]).groups())
        calls.append(a)
        return _Resp(BODY[a:b + 1])

    monkeypatch.setattr(fp, "urlopen", urlopen)
    return calls


def test_parallel_ranges_join_to_the_file_then_skip(tmp_path, server):
    assert fp.main(["--out", str(tmp_path), "--only", "fake", "--parallel", "3"]) == 0
    assert (tmp_path / "fake.h5ad").read_bytes() == BODY
    assert sorted(server) == [0, 8278, 16556]
    assert not list(tmp_path.glob("*.part*"))
    assert fp.main(["--out", str(tmp_path), "--only", "fake", "--parallel", "3"]) == 0
    assert len(server) == 3                       # second run fetched nothing


def test_resumes_each_range_where_it_stopped(tmp_path, server):
    (tmp_path / "fake.h5ad.part01").write_bytes(BODY[8278:8278 + 1000])
    assert fp.main(["--out", str(tmp_path), "--only", "fake", "--parallel", "3"]) == 0
    assert sorted(server) == [0, 9278, 16556]
    assert (tmp_path / "fake.h5ad").read_bytes() == BODY


def test_md5_mismatch_is_refused(tmp_path, server, monkeypatch):
    monkeypatch.setitem(fp.FILES, "fake", ("fake.h5ad", len(BODY), "0" * 32))
    assert fp.main(["--out", str(tmp_path), "--only", "fake", "--parallel", "2"]) == 1
    assert not (tmp_path / "fake.h5ad").exists()
    assert len(list(tmp_path.glob("fake.h5ad.part*"))) == 2


def test_bad_arguments_are_errors(tmp_path):
    with pytest.raises(SystemExit):
        fp.main(["--out", str(tmp_path), "--only", "nope"])
    with pytest.raises(SystemExit):
        fp.main(["--out", str(tmp_path), "--parallel", "0"])


def test_success_removes_every_part_and_joining_file_and_hints_on_parallel_mismatch(tmp_path, server, capsys):
    (tmp_path / "fake.h5ad.part05").write_bytes(b"left over from --parallel 8")
    (tmp_path / "fake.h5ad.part06").write_bytes(b"left over from --parallel 8")
    (tmp_path / "fake.h5ad.joining").write_bytes(b"half a join")
    assert fp.main(["--out", str(tmp_path), "--only", "fake", "--parallel", "3"]) == 0
    assert (tmp_path / "fake.h5ad").read_bytes() == BODY
    assert sorted(p.name for p in tmp_path.iterdir()) == ["fake.h5ad"]
    assert "same --parallel" in capsys.readouterr().err


def test_perturb_folder_has_a_source_in_the_data_manifest():
    spec = importlib.util.spec_from_file_location("data_manifest", ROOT / "benchmarks" / "data_manifest.py")
    dm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dm)
    assert "zenodo.org/records/7041849" in dm.SOURCES["perturb/"][0]


def test_join_appends_in_place_without_a_second_copy(tmp_path, server, monkeypatch):
    opened = []
    real_open = Path.open
    monkeypatch.setattr(Path, "open", lambda self, *a, **k: (opened.append(self.name), real_open(self, *a, **k))[1])
    assert fp.main(["--out", str(tmp_path), "--only", "fake", "--parallel", "3"]) == 0
    assert (tmp_path / "fake.h5ad").read_bytes() == BODY
    assert hashlib.md5(BODY).hexdigest() == fp.md5(tmp_path / "fake.h5ad")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["fake.h5ad"]
    assert not any(n.endswith(".joining") for n in opened)


def test_md5_mismatch_keeps_every_part_and_makes_no_final_file(tmp_path, server, monkeypatch):
    monkeypatch.setitem(fp.FILES, "fake", ("fake.h5ad", len(BODY), "0" * 32))
    written = []
    real_open = Path.open
    monkeypatch.setattr(Path, "open", lambda self, mode="r", *a, **k: (written.append(self.name) if "w" in mode else None, real_open(self, mode, *a, **k))[1])
    assert fp.main(["--out", str(tmp_path), "--only", "fake", "--parallel", "3"]) == 1
    assert not [n for n in written if not n.startswith("fake.h5ad.part")]   # the check writes nothing
    assert sorted(p.name for p in tmp_path.iterdir()) == [f"fake.h5ad.part{i:02d}" for i in range(3)]
    assert b"".join(p.read_bytes() for p in sorted(tmp_path.iterdir())) == BODY


def test_oserror_during_append_is_reported_and_leaves_no_final_file(tmp_path, server, monkeypatch, capsys):
    def boom(src, dst, *a):
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(fp, "shutil", type("S", (), {"copyfileobj": staticmethod(boom)}), raising=False)
    assert fp.main(["--out", str(tmp_path), "--only", "fake", "--parallel", "3"]) == 1
    assert not (tmp_path / "fake.h5ad").exists()
    err = capsys.readouterr().err
    assert "fake.h5ad" in err and "No space left" in err


def test_join_failing_after_some_parts_leaves_a_state_a_rerun_finishes(tmp_path, server, monkeypatch):
    real, n = fp.shutil.copyfileobj, []

    def second_fails(src, dst, *a):
        n.append(1)
        if len(n) == 2:
            raise OSError(28, "No space left on device")
        real(src, dst, *a)

    monkeypatch.setattr(fp.shutil, "copyfileobj", second_fails)
    assert fp.main(["--out", str(tmp_path), "--only", "fake", "--parallel", "3"]) == 1
    monkeypatch.setattr(fp.shutil, "copyfileobj", real)
    assert fp.main(["--out", str(tmp_path), "--only", "fake", "--parallel", "3"]) == 0
    assert (tmp_path / "fake.h5ad").read_bytes() == BODY
