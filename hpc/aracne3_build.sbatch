#!/bin/bash
#SBATCH --job-name=gb-aracne3-build
#SBATCH --partition=stud
#SBATCH --qos=stud
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.out
#
# Build ARACNe3 (native C++, CMake + OpenMP) inside a job: no compute on the login node.
# Download it first, on the login node (hpc/README.md, "ARACNe3"). The cluster has no cmake
# (no binary, no module): point CMAKE at the Kitware binary unpacked in ~/tools.
#
#   sbatch --export=ALL,CMAKE=$HOME/tools/cmake-3.31.8-linux-x86_64/bin/cmake hpc/aracne3_build.sbatch
#
# Prints the tool versions and the pinned commits, refuses an incomplete download, configures and
# builds a Release build, checks that build/src/app/ARACNe3_app_release exists, and calls it once
# without arguments: ARACNe3 then prints its usage line and exits 1 (src/app/ARACNe3.cpp).
# Home paths are printed as ~.
set -euo pipefail
SRC="${ARACNE3_DIR:-$HOME/tools/ARACNe3}"
CMAKE="${CMAKE:-cmake}"
JOBS="${SLURM_CPUS_PER_TASK:-4}"
REAL_HOME="$(readlink -f "$HOME")"
tilde() { sed -e "s#$REAL_HOME#~#g" -e "s#$HOME#~#g"; }

echo "cmake: $("$CMAKE" --version | head -1)"
echo "c++:   $(c++ --version | head -1)"
echo "ARACNe3 source: $(git -C "$SRC" rev-parse HEAD)"
git -C "$SRC" submodule status | tilde
incomplete=$(git -C "$SRC" submodule status --recursive | grep -c '^[-+U]' || true)
if [ "$incomplete" != 0 ]; then
  echo "$incomplete submodules missing or not at their pinned commit: finish the download on the login node" >&2
  exit 1
fi

"$CMAKE" -S "$SRC" -B "$SRC/build" -DCMAKE_BUILD_TYPE=Release 2>&1 | tilde
"$CMAKE" --build "$SRC/build" -j "$JOBS" 2>&1 | tilde

BIN="$SRC/build/src/app/ARACNe3_app_release"
test -x "$BIN" || { echo "no build/src/app/ARACNe3_app_release after the build" >&2; exit 1; }
echo "program sha256: $(sha256sum "$BIN" | cut -d' ' -f1)"
ldd "$BIN" | grep -E 'gomp|stdc\+\+' | tilde || true
rc=0
"$BIN" > "$SRC/build/usage.txt" 2>&1 || rc=$?
echo "called without arguments: exit status $rc"
tilde < "$SRC/build/usage.txt"
