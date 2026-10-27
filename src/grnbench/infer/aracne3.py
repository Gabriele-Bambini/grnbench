"""ARACNe3 adapter: drives the native ARACNe3 program, with honest availability probing.

ARACNe3 (github.com/califano-lab/ARACNe3) is a C++ program built with CMake and OpenMP. It is not
a jar: the Java tool is its predecessor ARACNe-AP, which this adapter does not drive. Written
against upstream commit 3d8791a (2025-01-16): the README, ``src/app/ARACNe3.cpp`` (the flags) and
``src/app/io.cpp`` (the file formats). How to build it: ``hpc/README.md``.

Where the program is: ``GRNBENCH_ARACNE3_BIN``, else ``ARACNe3_app_release`` (the README's Release
build) or ``ARACNe3_app`` (CMake's name without a build type) on PATH. Without it the method is
reported unavailable and running it raises :class:`MethodError` naming the missing piece.

One call is one ARACNe3 run in a fresh temporary folder::

    ARACNe3_app_release -e expr.tsv -r regulators.txt -o out \\
        --seed 42 --threads 1 --runid grnbench [extra_args]

- ``expr.tsv``: genes as rows, cells as columns, a header row of cell ids (ARACNe3 counts samples
  from its tabs). ARACNe3 rank-transforms every gene itself, so the normalized matrix goes as is.
- ``regulators.txt``: the TFs that are genes of the matrix, one per line (every gene without a list).
- The temporary folder is also the working folder: ARACNe3 caches its null model in
  ``./.ARACNe3_cached/`` and, when a cached one exists, skips the seeded draws that build it, so a
  shared folder would make a run depend on the runs before it.
- A network file system can refuse to delete that folder right after a finished run (``EBUSY`` on the
  cluster's BeeGFS home in EXP-20261019-GB-1, which lost a finished network that way). The network is
  then returned anyway, the folder's name goes into ``attrs["temporary_folder_left"]`` (the run's
  provenance) and a ``RuntimeWarning`` says it is left; a run that failed still reports its own error.
  Where to put ``TMPDIR`` on the cluster: ``hpc/README.md``.
- ARACNe3 runs with one thread unless asked for more: in ``src/app/apmi_nullmodel.cpp`` the null-model
  loop shuffles one vector shared by all threads and reads it outside the critical section, so with more
  threads two runs with the same seed can differ. Measured in EXP-20261018-GB-1 (500 cells x 300 genes):
  with 4 threads the cached null model differed between two same-seed runs, while the network came out
  byte-identical; with 1 thread both repeated. So the engine's generic ``threads`` (``grnbench infer`` and
  ``grnbench run`` pass 4, ``hpc/grnbench.sbatch`` the job's CPU count, to every method) is accepted and
  not passed on; only ``aracne3_threads`` sets ARACNe3's ``--threads``. ``attrs["threads"]`` (the run's
  provenance, ``method_info``) is the count ARACNe3 got; ``attrs["threads_not_used"]`` names a generic
  count set aside.
- ``extra_args`` may add flags ARACNe3 reads (``-x 30``, ``--alpha 0.01``, ``--noMaxEnt``, ...).
  ARACNe3 ignores flags it does not know without a word, so those are refused, and so are the flags
  the adapter sets itself. A keyword the adapter does not take is refused too, instead of dropped.
- A failed run's error carries the end of ARACNe3's output, which names the files it was given: the
  temporary folder is written ``<tmp>`` and the program's path its file name, because the error goes
  into run manifests.

Output: ``out/consolidated-net_grnbench.tsv``. weight = MI of the pair on the full matrix,
direction = sign of its Spearman correlation; the subnetwork count and p-value are not kept.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ENV_BIN = "GRNBENCH_ARACNE3_BIN"
BINARY_NAMES = ("ARACNe3_app_release", "ARACNe3_app")
RUN_ID = "grnbench"
CONSOLIDATED_COLUMNS = ["regulator.values", "target.values", "mi.values",
                        "scc.values", "count.values", "log.p.values"]
# ARACNe3 keeps sample counts and gene ids in 16 bits: more cells end in a row-length error, more
# genes wrap the gene id without a word. Both are refused before the run.
MAX_UINT16 = 65_535

# the flags main() reads (src/app/ARACNe3.cpp)
_SET_BY_ADAPTER = {"-e", "-r", "-o", "--seed", "--threads", "--runid", "--consolidate", "--noConsolidate",
                   "-h", "--help"}
_WITH_VALUE = {"-x", "--alpha", "--subsample", "--min-subnets", "--mithresh", "--numnulls"}
_SWITCHES = {"--noalpha", "--noMaxEnt", "--FDR", "--FWER", "--FPR", "--adaptive", "--verbose"}

# the keywords run_aracne3 takes; ``seed`` and ``threads`` are the ones the engine's callers pass to
# every method (cli_infer, beeline_backtest, stability, community), ``tf_names`` comes from run_method
_PARAMETERS = ("tf_names", "seed", "threads", "aracne3_threads", "timeout_s", "extra_args")


def _find_binary() -> tuple[str | None, str]:
    """(path of the ARACNe3 program or None, why not). The reason never repeats the path: it is
    written into run manifests."""
    raw = os.environ.get(ENV_BIN, "")
    if raw:
        if not Path(raw).is_file():
            return None, f"{ENV_BIN} does not point at a file; set it to the built ARACNe3_app_release"
        if not os.access(raw, os.X_OK):
            return None, f"{ENV_BIN} points at a file that is not executable (chmod +x it)"
        # absolute: the program runs inside its own temporary folder, where a relative path means nothing
        return os.path.abspath(raw), ""
    for name in BINARY_NAMES:
        found = shutil.which(name)
        if found:
            return found, ""
    why = (f"ARACNe3 program not found: set {ENV_BIN} to the built ARACNe3_app_release or put it on "
           "PATH (native C++ build of github.com/califano-lab/ARACNe3, see hpc/README.md)")
    if os.environ.get("ARACNE3_JAR"):
        why += "; ARACNE3_JAR is no longer used (ARACNe3 is a native program, not a jar)"
    return None, why


def aracne3_availability() -> tuple[bool, str]:
    binary, why = _find_binary()
    return binary is not None, why


def _checked_extra_args(extra_args: list[str] | None) -> list[str]:
    """``extra_args`` as strings, refusing what ARACNe3 would misread or silently ignore."""
    from grnbench.infer.base import MethodError

    args = [str(a) for a in (extra_args or [])]
    i = 0
    while i < len(args):
        flag = args[i]
        if flag in _SET_BY_ADAPTER:
            raise MethodError(f"aracne3: {flag} cannot be passed in extra_args (the adapter sets the input "
                              "files, output folder, run id, seed and threads, and always consolidates; "
                              "aracne3_threads sets the threads)")
        if flag in _WITH_VALUE:
            if i + 1 == len(args):
                raise MethodError(f"aracne3: {flag} needs a value")
            i += 2
            continue
        if flag not in _SWITCHES:
            known = sorted(_WITH_VALUE | _SWITCHES)
            near = [k for k in known if k.lstrip("-").lower() == flag.lstrip("-").lower()]
            hint = f"; upstream spells it {near[0]}" if near else ""
            raise MethodError(f"aracne3 does not read {flag!r} and would ignore it silently{hint}; "
                              f"flags it reads: {' '.join(known)}")
        i += 1
    return args


def _without_paths(text: str, folder: str, binary: str) -> str:
    """ARACNe3's output with the temporary folder written ``<tmp>`` and the program's path its file name."""
    for path in sorted({folder, os.path.realpath(folder)}, key=len, reverse=True):
        text = text.replace(path, "<tmp>")
    return text.replace(binary, Path(binary).name)


def run_aracne3(expr: pd.DataFrame,
                tf_names: list[str] | None = None,
                seed: int = 42,
                threads: int | None = None,
                timeout_s: int = 3600,
                extra_args: list[str] | None = None,
                aracne3_threads: int = 1,
                **unknown) -> pd.DataFrame:
    """Run ARACNe3 once on ``expr`` (cells x genes); return [source, target, weight, direction].

    ``threads`` is the engine's generic thread count: accepted, recorded, not passed to ARACNe3.
    ``aracne3_threads`` is ARACNe3's ``--threads`` (module docstring).
    """
    from grnbench.infer.base import MethodError

    if unknown:
        raise MethodError(f"aracne3 does not take {', '.join(sorted(unknown))}; its parameters: "
                          f"{', '.join(_PARAMETERS)} (threads is the engine's generic count and is not "
                          "passed to ARACNe3: aracne3_threads sets its --threads)")
    n_threads = int(aracne3_threads)
    if n_threads < 1:
        raise MethodError(f"aracne3_threads must be at least 1, got {aracne3_threads!r}")
    binary, why = _find_binary()
    if binary is None:
        raise MethodError(f"aracne3 unavailable: {why}")
    n_cells, n_genes = expr.shape
    if n_cells > MAX_UINT16 or n_genes > MAX_UINT16:
        raise MethodError(f"aracne3 reads at most {MAX_UINT16:,} cells and {MAX_UINT16:,} genes (16-bit "
                          f"counters); got {n_cells:,} cells x {n_genes:,} genes. Subsample the cells or "
                          "keep fewer genes first.")
    genes = [str(g) for g in expr.columns]
    wanted = [str(t) for t in tf_names] if tf_names else genes
    in_matrix = set(genes)
    regulators = [t for t in dict.fromkeys(wanted) if t in in_matrix]
    if not regulators:
        raise MethodError(f"aracne3: none of the {len(wanted)} regulators is a gene of the expression matrix")
    extra = _checked_extra_args(extra_args)

    def run_in(td: str) -> pd.DataFrame:
        expr_path, reg_path, out_dir = Path(td) / "expr.tsv", Path(td) / "regulators.txt", Path(td) / "out"
        # nine significant digits: rounding must not create ties that ARACNe3 would break at random
        expr.T.to_csv(expr_path, sep="\t", float_format="%.9g", lineterminator="\n")
        reg_path.write_text("".join(f"{r}\n" for r in regulators), encoding="utf-8", newline="\n")
        cmd = [binary, "-e", str(expr_path), "-r", str(reg_path), "-o", str(out_dir),
               "--seed", str(int(seed)), "--threads", str(n_threads), "--runid", RUN_ID, *extra]
        try:
            proc = subprocess.run(cmd, cwd=td, capture_output=True, text=True,
                                  timeout=int(timeout_s))  # noqa: S603 - fixed argv list
        except subprocess.TimeoutExpired:
            raise MethodError(f"aracne3 timed out after {int(timeout_s)} s") from None
        # paths hidden before the cut: a path cut in half would no longer match
        tail = _without_paths((proc.stderr or proc.stdout or "").strip(), td, binary)[-800:]
        if proc.returncode != 0:
            raise MethodError(f"aracne3 failed (rc={proc.returncode}): {tail}")
        net_path = out_dir / f"consolidated-net_{RUN_ID}.tsv"
        if not net_path.exists():
            raise MethodError(f"aracne3 exited 0 but wrote no {net_path.name}: {tail}")
        net = parse_aracne3_text(net_path.read_text(encoding="utf-8", errors="replace"))
        net.attrs["threads"] = n_threads
        if threads is not None and int(threads) != n_threads:
            net.attrs["threads_not_used"] = (f"threads={threads} (the engine's generic count) was not passed "
                                             "to ARACNe3: its null model repeats from the seed only with one "
                                             "thread; aracne3_threads sets its --threads")
        return net

    edges, run_error = None, None
    try:
        with tempfile.TemporaryDirectory(prefix="grnbench_aracne3_") as td:
            try:
                edges = run_in(td)
            except BaseException as err:  # kept, so a folder that cannot be deleted either never hides it
                run_error = err
                raise
    except OSError as exc:
        if run_error is not None:
            raise run_error
        if edges is None:                 # the temporary folder could not even be made
            raise
        # the folder could not be deleted after a finished run (module docstring): keep the network
        left = f"{Path(td).name} ({type(exc).__name__}: [Errno {exc.errno}] {exc.strerror})"
        edges.attrs["temporary_folder_left"] = left
        warnings.warn(f"ARACNE3 FOLDER LEFT: the network is kept, but the temporary folder {left} could not "
                      "be deleted; remove it from the temporary directory", RuntimeWarning, stacklevel=2)
    return edges


def parse_aracne3_text(text: str) -> pd.DataFrame:
    """Parse ARACNe3's consolidated network into [source, target, weight, direction].

    weight = ``mi.values`` (MI of the pair on the full matrix); direction = sign of ``scc.values``
    (its Spearman correlation, 0 when exactly 0). Any other table is refused: guessing columns is how
    a wrong file turns into a network.
    """
    from grnbench.infer.base import MethodError

    lines = text.splitlines()
    header = lines[0].split("\t") if lines else []
    if header != CONSOLIDATED_COLUMNS:
        raise MethodError("not an ARACNe3 consolidated network: expected the columns "
                          f"{' '.join(CONSOLIDATED_COLUMNS)}, got {' '.join(header)[:200] or 'an empty file'}")
    rows = [line.split("\t") for line in lines[1:] if line.strip()]
    if not rows:
        raise MethodError("ARACNe3 kept no edges: its pruning (FDR, FWER or FPR threshold, then MaxEnt) "
                          "removed every pair; relax it with extra_args, e.g. ['--alpha', '0.2'] or ['--noMaxEnt']")
    bad = next((n for n, row in enumerate(rows, start=2) if len(row) != len(CONSOLIDATED_COLUMNS)), None)
    if bad is not None:
        raise MethodError(f"ARACNe3 consolidated network: line {bad} does not have "
                          f"{len(CONSOLIDATED_COLUMNS)} tab-separated fields")
    table = pd.DataFrame(rows, columns=CONSOLIDATED_COLUMNS)
    return pd.DataFrame({
        "source": table["regulator.values"],
        "target": table["target.values"],
        "weight": table["mi.values"].astype(float),
        "direction": np.sign(table["scc.values"].astype(float)).astype(int),
    })
