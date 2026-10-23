# Install the R packages of the optional grnbench R bridge (docs/R_BRIDGE.md) into the R library of the conda
# env grnbench-r, and nowhere else. It installs only what is missing and upgrades nothing.
#   Linux:   conda env create -f requirements/r-bridge.yml
#            conda run -n grnbench-r Rscript --vanilla requirements/r-bridge.R    (SoupX from CRAN)
#   Windows: conda create -n grnbench-r -c conda-forge r-base=4.5.3 r-biocmanager r-jsonlite
#            conda run -n grnbench-r Rscript --vanilla requirements/r-bridge.R    (CRAN/Bioconductor binaries)
# The approved set: DropletUtils, scDblFinder, SoupX, celda, miQC, jsonlite (plus their dependencies).
# It prints each package's version and exits 1 if any of them does not load from the env library.
windows <- .Platform$OS.type == "windows"
lib <- .Library                         # the env's own library: never a personal R library
canon <- function(p) {
  p <- sub("/+$", "", normalizePath(p, winslash = "/", mustWork = FALSE))
  if (windows) tolower(p) else p
}
prefix <- Sys.getenv("CONDA_PREFIX")
if (!nzchar(prefix) || !startsWith(canon(lib), paste0(canon(prefix), "/"))) {
  message("r-bridge.R installs only into the R library of the conda env grnbench-r. Run it as\n",
          "  conda run -n grnbench-r Rscript --vanilla requirements/r-bridge.R\n",
          "(CONDA_PREFIX = '", prefix, "'; this R's library = '", lib, "')")
  quit(save = "no", status = 1L)
}
.libPaths(lib)
if (!nzchar(system.file(package = "BiocManager", lib.loc = lib))) {
  message("BiocManager is missing: create the env as above (it installs r-biocmanager), then run this again")
  quit(save = "no", status = 1L)
}
options(repos = BiocManager::repositories(), Ncpus = 4L)
cat(R.version.string, "\nlibrary:", lib, "\n")
absent <- function(pkgs) pkgs[!vapply(pkgs, function(p) nzchar(system.file(package = p, lib.loc = lib)), logical(1))]

if (windows) {
  cran <- absent(c("SoupX", "Matrix"))
  if (length(cran)) install.packages(cran, lib = lib, type = "binary")
  bioc <- absent(c("DropletUtils", "scDblFinder", "celda", "miQC"))
  if (length(bioc)) BiocManager::install(bioc, lib = lib, ask = FALSE, update = FALSE, type = "binary")
  # bluster (needed by scDblFinder) has no Windows binary for Bioconductor 3.22 yet, so the binary install
  # above leaves it out; the Bioconductor 3.21 build (1.17.0) works with R 4.5. Its own dependencies come
  # from the usual repositories.
  if (length(absent("bluster"))) {
    install.packages("bluster", lib = lib, type = "win.binary",
                     repos = c(bioc321 = "https://bioconductor.org/packages/3.21/bioc", getOption("repos")))
  }
} else {
  cran <- absent("SoupX")      # pure R (NeedsCompilation: no); everything else comes from r-bridge.yml
  if (length(cran)) install.packages(cran, lib = lib)
}

pkgs <- c("SoupX", "Matrix", "DropletUtils", "scDblFinder", "bluster", "celda", "miQC", "jsonlite", "BiocManager")
ok <- vapply(pkgs, function(p) {
  loads <- tryCatch({ loadNamespace(p, lib.loc = lib); TRUE },
                    error = function(e) { message(p, ": ", conditionMessage(e)); FALSE })
  cat(sprintf("%-14s %s\n", p, if (loads) as.character(utils::packageVersion(p, lib.loc = lib)) else "DOES NOT LOAD"))
  loads
}, logical(1))
cat("Bioconductor", as.character(BiocManager::version()), "\n")
if (!all(ok)) quit(save = "no", status = 1L)
