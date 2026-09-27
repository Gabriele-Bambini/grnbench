"""hpc/hpc.sh maps ssh failures to the one action the person or agent should take next.
ssh is replaced by a stub that prints a given error and exits 255; no network is used."""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _check(tmp_path: Path, ssh_error: str, os_name: str) -> subprocess.CompletedProcess:
    bash = shutil.which("bash") or "C:/Program Files/Git/bin/bash.exe"
    if not Path(bash).exists():
        pytest.skip("Bash unavailable")
    stub = tmp_path / "ssh"
    stub.write_text(f"#!/bin/bash\necho '{ssh_error}' >&2\nexit 255\n", encoding="utf-8", newline="\n")
    stub.chmod(0o755)
    env = {**os.environ, "SSH": stub.as_posix(), "GRNBENCH_HPC_HOST": "alias", "GRNBENCH_OS": os_name}
    return subprocess.run([bash, str(ROOT / "hpc" / "hpc.sh"), "check"], env=env,
                          capture_output=True, text=True)


@pytest.mark.parametrize("os_name, expected", [
    ("MINGW64_NT-10.0-26200", "powershell -NoProfile -ExecutionPolicy Bypass -File hpc/vpn.ps1 reconnect"),
    ("Linux", "globalprotect connect --portal vpnstudents.unibocconi.it"),
    ("Darwin", "Open the GlobalProtect app"),
])
def test_unresolved_login_node_prints_the_vpn_command_first(tmp_path, os_name, expected):
    r = _check(tmp_path, "ssh: Could not resolve hostname slogin.example: Name or service not known", os_name)
    assert r.returncode == 255
    first = r.stderr.splitlines()[0]
    assert expected in first, r.stderr
    assert "Could not resolve hostname" in r.stderr


def test_timeout_also_points_at_the_vpn(tmp_path):
    r = _check(tmp_path, "ssh: connect to host slogin.example port 22: Connection timed out", "Linux")
    assert r.stderr.startswith("VPN down?"), r.stderr


def test_auth_failure_is_not_blamed_on_the_vpn(tmp_path):
    r = _check(tmp_path, "user@host: Permission denied (publickey).", "MINGW64_NT-10.0")
    assert "VPN down?" not in r.stderr
    assert "Auth failed" in r.stderr and "Do not retry in a loop" in r.stderr


def test_lockout_says_stop_and_wait(tmp_path):
    r = _check(tmp_path, "Not allowed at this time", "Linux")
    assert "LOCKED OUT" in r.stderr and "VPN down?" not in r.stderr


def test_remote_command_timeouts_are_not_blamed_on_the_vpn(tmp_path):
    """A remote command's own 'timed out' (e.g. pip) exits with its code, not ssh's 255."""
    bash = shutil.which("bash") or "C:/Program Files/Git/bin/bash.exe"
    stub = tmp_path / "ssh"
    stub.write_text("#!/bin/bash\necho 'ReadTimeoutError: Read timed out.' >&2\nexit 1\n",
                    encoding="utf-8", newline="\n")
    stub.chmod(0o755)
    env = {**os.environ, "SSH": stub.as_posix(), "GRNBENCH_HPC_HOST": "alias", "GRNBENCH_OS": "Linux"}
    r = subprocess.run([bash, str(ROOT / "hpc" / "hpc.sh"), "run", "pip install x"], env=env,
                       capture_output=True, text=True)
    assert r.returncode == 1 and "VPN down?" not in r.stderr and "Read timed out" in r.stderr
