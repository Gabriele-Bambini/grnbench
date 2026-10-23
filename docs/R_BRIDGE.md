# The optional R bridge

Several cleaning options of the knowledge base exist only as R packages: EmptyDrops (DropletUtils),
SoupX, DecontX (celda), scDblFinder, DoubletFinder, scran, sctransform, Slingshot. The R bridge
(`src/grnbench/rbridge.py`, decision D-032) lets grnbench run them without making R a dependency of
the Python package:

- R runs as `Rscript --vanilla` in a child process. Nothing R-related is imported into Python (no rpy2).
- R lives in a dedicated conda env, `grnbench-r`, built from `requirements/`. The Python lock and the
  registered results do not change.
- grnbench never installs anything. It finds R, asks it which packages it has, and tells you the
  command to run for what is missing.
- Every R run goes through `src/grnbench/r/_wrapper.R`, started in the run's work folder with nothing
  else on its command line. Everything it needs is in `params.json` in that folder (UTF-8, so paths with
  any characters work on Windows): the method script, `seed`, `params`, `inputs`, `expect` and
  `workdir`. It is read with `jsonlite::fromJSON(simplifyVector = FALSE)`: JSON arrays arrive in R as
  lists, scalars as scalars. Parameters must be plain JSON values: NaN and infinity are refused.
- The wrapper fixes the random number generator (`RNGkind("Mersenne-Twister", "Inversion",
  "Rejection")`) and the seed, and writes `session.json`: R version, platform, library folders, the RNG
  kind it set (`rng_kind`) and the one in force after the method (`rng_kind_after`), the method's own
  `set.seed(...)` calls (`set_seed_calls`), the seed and the version of every package that was loaded.
  Only calls written `set.seed` in the method script are seen: one written `base::set.seed`, or made
  inside a package, is not listed (a changed generator kind still shows in `rng_kind_after`). A run
  manifest (`RResult.manifest_block()`) reports `rng_changed` when the method changed the generator or
  reseeded as above, and `skipped`: the R candidates passed over before the one used, and why.
  Everything R prints lands in `r.log`.
- A run that times out, or that you interrupt (Ctrl+C), is killed with all its child processes; if the
  kill cannot be confirmed, the error says `process may still be running (pid N)`. Killing the Python
  process itself from outside (for example a SIGTERM to it) can still leave R running.

## Install

Run these from a clone of the repository: `requirements/` is not part of the installed package.
The env gets only the approved set: DropletUtils, scDblFinder, SoupX, celda, miQC and jsonlite, plus
their dependencies (scran and sctransform come in that way). DoubletFinder and Slingshot are not
installed: their options have no grnbench executor yet, and installing them needs the owner's OK.

**Linux (and the cluster).** Installing R on the cluster needs the owner's OK first.

```bash
conda env create -f requirements/r-bridge.yml
conda run -n grnbench-r Rscript --vanilla requirements/r-bridge.R
```

The yml pins R 4.5.3 and the Bioconductor packages from conda-forge and bioconda (no `defaults`
channel). `r-bridge.R` then adds SoupX 1.6.2 from CRAN (conda-forge builds `r-soupx` only for
R 4.2/4.3; SoupX is pure R, and its imports ggplot2, Matrix and Seurat come from the env).

**Windows.** bioconda has no win-64 builds, so the env holds R only and the packages come as CRAN and
Bioconductor binaries:

```bash
conda create -n grnbench-r -c conda-forge r-base=4.5.3 r-biocmanager r-jsonlite
conda run -n grnbench-r Rscript --vanilla requirements/r-bridge.R
```

Bioconductor 3.22 has no Windows binary of `bluster` yet (scDblFinder needs it). `r-bridge.R`
then installs the Bioconductor 3.21 build (1.17.0), which works with R 4.5.

`r-bridge.R` installs into the env's own R library (`.Library`) and nowhere else: run outside the
env (no `CONDA_PREFIX`, or an R whose library is not inside it) it exits 1 before installing anything.
It installs only what is missing and upgrades nothing, so running it again is safe. It prints every
package's version and exits 1 if one of them does not load from the env library.

## Which R is used

`grnbench.rbridge.find_rscript()` takes the first of:

1. an Rscript given explicitly (`grnbench env --rscript PATH`);
2. `GRNBENCH_RSCRIPT`: the path of an Rscript. The value `none` (any case) turns R off: grnbench then
   behaves as if no R existed and starts nothing;
3. the conda env named by `GRNBENCH_R_ENV` (default `grnbench-r`). Its folder is read from
   `~/.conda/environments.txt`; when that file does not list it, from `conda env list --json` (conda,
   mamba or micromamba, found through `CONDA_EXE` or PATH);
4. `$R_HOME/bin/Rscript`;
5. `Rscript` on PATH.

grnbench never quietly switches to another R:

- an Rscript you named (1 or 2) that does not work is an error;
- `GRNBENCH_R_ENV` naming an env that does not exist is an error;
- a conda env that exists (named, or the default `grnbench-r`) but whose R is missing or does not
  answer is an error, even when another R is on PATH.

Only when there is no env called `grnbench-r` at all does grnbench go on to 4 and 5. A candidate that
is passed over is listed, with the reason, by `grnbench env` and in the run manifest (`skipped`).
Relative paths are made absolute first (R runs in the work folder).

**Libraries.** A conda R child gets a clean library setting: `R_LIBS` and `R_LIBS_SITE` are removed
(another R's packages must not leak in) and `R_LIBS_USER` is the env's library. Any other R keeps its
own library settings, since grnbench does not know the layout of that installation. `R_HOME` is
removed unless that is where R was found. An Rscript counts as a conda R when a parent folder (up to
five levels) holds `conda-meta`, including one given through `GRNBENCH_RSCRIPT`.

**Windows PATH detail.** A conda R on Windows (`<env>\Scripts\Rscript.exe`, and also
`<env>\lib\R\bin\x64\Rscript.exe`) exits at once with code 0xC0000135 (a DLL was not found) unless
PATH starts with the env's own folders:
`<env>;<env>\Library\mingw-w64\bin;<env>\Library\usr\bin;<env>\Library\bin;<env>\Scripts;<env>\bin`.
grnbench puts them in front of PATH for every conda R child it starts. It does not use `conda run`,
which is slow and holds back R's output until the end.

## What R has: the probe and its cache

Asking R which packages it has is one `Rscript` call (`src/grnbench/r/_probe.R`, about 0.2 s; the
timeout is 30 s for a cold start behind an antivirus scan). It reads the packages' DESCRIPTION files
and loads nothing. A package built for another R major.minor than the running R counts as not usable
(`built for R x.y`), with the command to reinstall it.

The answer is remembered in the process and on disk, in `GRNBENCH_CACHE_DIR` if set, else
`%LOCALAPPDATA%\grnbench\cache` on Windows, else `$XDG_CACHE_HOME/grnbench` or `~/.cache/grnbench`.
Worker processes therefore do not start R again. An answer is keyed by the Rscript (path and
modification time), the package names and, for an R outside conda, your `R_LIBS_USER`, `R_LIBS` and
`R_LIBS_SITE`. It is reused while every library folder R reported keeps the modification time it had
before R was asked: installing a package changes its folder, so the next check asks R again. A failed
`conda env list` is never remembered. A cache that cannot be read or written is ignored; deleting it
is always safe.

## grnbench env

```bash
grnbench env            # table: tool, found, version, source, what needs it, install command
grnbench env --json     # the same rows as JSON
grnbench env --rscript ~/anaconda3/envs/grnbench-r/Scripts/Rscript.exe
```

One row covers R itself, one each `r:` package of the knowledge base (plus jsonlite, which the wrapper
needs), and one each `py:` module and `cli:` tool. It always asks R again (no cache): it is a
diagnostic. A missing R package shows the line to run in that R, from the repository it comes from:
`BiocManager::install("...")` (Bioconductor), `install.packages("...")` (CRAN) or
`remotes::install_github("...")` (DoubletFinder); for the conda env the line adds `lib = .Library`.
Missing R shows the commands above.

It exits 0 when R is absent or turned off (that is an answer), and 1 when looking for R is an error,
exactly when `run_r` would refuse: the R you named (`--rscript`, `GRNBENCH_RSCRIPT`) does not work,
`GRNBENCH_R_ENV` names no env, or an env that exists has no working R. Scripts can check it.
The install lines for packages outside the approved set (DoubletFinder, Slingshot) are for reference:
adding them to `grnbench-r` needs the owner's OK.
In the interface, the sidebar's **Environment** popover shows the same rows; R starts only when you
press **Check R**, and again each time you press it.

## Not wired yet (step 2c)

The knowledge base still decides availability on its own (`kb.Option.availability`, unchanged), so no
R option can be chosen by the planner yet, and planner output is the same as before. Step 2c
(after PR #133) connects `requirement_status` to the knowledge base, keeping the reason text
`missing: r:<package>`. Under the current preset, R options will be opt-in.
