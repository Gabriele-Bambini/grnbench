#!/bin/bash
# Fetch the BEELINE benchmark corpus (Pratapa et al., Nature Methods 2020).
#
#   bash hpc/fetch_beeline.sh [DEST]        # default: $HOME/beeline
#
# RUN THIS ON THE LOGIN NODE: compute nodes usually have no outbound network.
# ~280 MB total, pinned Zenodo record 19009603 (v4, BEELINE v1.1). Existing complete
# archives are not re-downloaded.
set -euo pipefail

DEST="${1:-$HOME/beeline}"
BASE="https://zenodo.org/records/19009603/files"
command -v md5sum >/dev/null 2>&1 || { echo "md5sum required for archive verification" >&2; exit 1; }
mkdir -p "$DEST" || exit 1
cd "$DEST" || exit 1

fetch() {
  local name="$1" md5="$2"
  if [ -s "$name" ]; then
    if command -v md5sum >/dev/null 2>&1 && \
       [ "$(md5sum "$name" | cut -d' ' -f1)" = "$md5" ]; then
      echo "  $name already present and intact"
      return 0
    fi
    echo "  $name present but incomplete/corrupt - refetching"
    rm -f "$name"
  fi
  echo "  downloading $name ..."
  # --retry-all-errors + resume: Zenodo has been observed truncating this
  # transfer while the client still reports success, so the md5 check below is
  # the only thing standing between you and a corrupt archive.
  if command -v curl >/dev/null 2>&1; then
    curl -fSL --retry 8 --retry-delay 5 --retry-all-errors -C - \
         -o "$name" "$BASE/$name?download=1" || return 1
  else
    wget -q --tries=8 --continue -O "$name" "$BASE/$name?download=1" || return 1
  fi
  if command -v md5sum >/dev/null 2>&1 && \
     [ "$(md5sum "$name" | cut -d' ' -f1)" != "$md5" ]; then
    echo "  ERROR: $name failed its md5 check after download (truncated?)" >&2
    return 1
  fi
}

echo "=== fetching BEELINE into $DEST ==="
fetch BEELINE-Networks.zip a7dca84d302bb127bad962ae52815662 || { echo FETCH_FAIL; exit 1; }
fetch BEELINE-data.zip     9ea93065004f2f12022adeeaf3728edd || { echo FETCH_FAIL; exit 1; }

echo "=== unpacking ==="
for z in BEELINE-Networks.zip BEELINE-data.zip; do
  d="${z%.zip}"
  [ -d "$d" ] || { mkdir -p "$d" && unzip -q -o "$z" -d "$d"; }
done

echo "=== layout ==="
find "$DEST" -maxdepth 3 -type d | sed -n '1,30p'
echo
echo "ground-truth networks:"
find "$DEST" -name "*-network.csv" | sed "s|$DEST/||" | sort
echo
echo "expression datasets (ExpressionData.csv + GeneOrdering.csv):"
find "$DEST" -name ExpressionData.csv -printf '%h\n' 2>/dev/null | sort | sed -n '1,20p'

echo
echo "Now submit the backtest:"
echo "  sbatch --export=ALL,BEELINE=$DEST hpc/beeline_backtest.sbatch"
echo FETCH_DONE
