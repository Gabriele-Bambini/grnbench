# grnbench R bridge wrapper (grnbench.rbridge.run_r): Rscript --vanilla _wrapper.R, started in the work folder.
# Everything else comes from <work folder>/params.json, which is UTF-8, so paths with any characters survive
# (on Windows the command line is limited to the system code page): the method script, seed, params, inputs,
# expect and workdir. It is read with jsonlite::fromJSON(simplifyVector = FALSE): JSON arrays arrive as R
# lists, scalars as scalars.
# The wrapper fixes the random number generator and the seed, then sources the method script in a fresh
# environment that has `params`, `inputs`, `workdir` and `expect` defined, and a `set.seed` that records each
# call before calling base::set.seed (a call written base::set.seed, or made inside a package, is not seen;
# a changed generator kind still shows in rng_kind_after).
# It writes <workdir>/session.json: R and package versions, the RNG kind it
# set (rng_kind) and the one in force after the method (rng_kind_after), the method's own set.seed calls, the
# seed and the seconds. Exit status: 0 done; 1 the method failed or session.json could not be written;
# 2 bad call (arguments, no params.json, no method script, jsonlite missing).
options(warn = 1)                       # warnings go to the log where they happen, not at the end
if (length(commandArgs(trailingOnly = TRUE))) {
  message("grnbench wrapper: takes no arguments; everything is in params.json in the working folder")
  quit(save = "no", status = 2L)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  message("grnbench wrapper: the R package jsonlite is missing in ", R.version.string, " (",
          paste(.libPaths(), collapse = "; "), "). Install it in this R with install.packages(\"jsonlite\"), ",
          "or create the grnbench-r env (docs/R_BRIDGE.md).")
  quit(save = "no", status = 2L)
}
fail <- function(what, e) {
  message("grnbench wrapper: ", what, ": ", conditionMessage(e))
  call <- conditionCall(e)
  if (!is.null(call)) message("  in: ", paste(deparse(call, nlines = 3L), collapse = " "))
  quit(save = "no", status = 1L)
}
started <- proc.time()[["elapsed"]]
params_file <- file.path(getwd(), "params.json")
if (!file.exists(params_file)) {
  message("grnbench wrapper: no params.json in the working folder ", getwd())
  quit(save = "no", status = 2L)
}
run <- tryCatch(jsonlite::fromJSON(params_file, simplifyVector = FALSE),
                error = function(e) fail("cannot read params.json", e))
script <- run$script
if (!is.character(script) || length(script) != 1L || !file.exists(script)) {
  message("grnbench wrapper: the method script named in params.json does not exist: ",
          paste(format(script), collapse = " "))
  quit(save = "no", status = 2L)
}
workdir <- run$workdir
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(run$seed)
rng <- RNGkind()
seed_calls <- list()
method_env <- new.env(parent = globalenv())
method_env$params <- run$params
method_env$inputs <- run$inputs
method_env$workdir <- workdir
method_env$expect <- run$expect
method_env$set.seed <- function(seed, kind = NULL, normal.kind = NULL, sample.kind = NULL) {
  or_na <- function(x) if (is.null(x)) NA else x
  seed_calls[[length(seed_calls) + 1L]] <<- list(seed = or_na(seed), kind = or_na(kind),
                                                 normal.kind = or_na(normal.kind), sample.kind = or_na(sample.kind))
  base::set.seed(seed, kind = kind, normal.kind = normal.kind, sample.kind = sample.kind)
}
# UTF-8-BOM reads UTF-8 with or without the byte-order mark Windows editors add
tryCatch(source(script, local = method_env, encoding = "UTF-8-BOM"),
         error = function(e) fail(paste("error in", basename(script)), e))
rng_after <- RNGkind()
loaded <- sort(loadedNamespaces())
base_r <- vapply(loaded, function(p) {
  identical(suppressWarnings(utils::packageDescription(p, fields = "Priority")), "base")
}, logical(1))
versions <- lapply(loaded[!base_r], function(p) as.character(getNamespaceVersion(p)))
names(versions) <- loaded[!base_r]
if (!length(versions)) versions <- structure(list(), names = character(0))   # {} rather than []
session <- list(r_version = R.version.string, platform = R.version$platform, lib_paths = I(.libPaths()),
                rng_kind = I(rng), rng_kind_after = I(rng_after), set_seed_calls = seed_calls, seed = run$seed,
                packages = versions, seconds = round(proc.time()[["elapsed"]] - started, 3))
target <- file.path(workdir, "session.json")
written <- tryCatch({
  jsonlite::write_json(session, paste0(target, ".tmp"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
  file.rename(paste0(target, ".tmp"), target)
}, error = function(e) { message("grnbench wrapper: ", conditionMessage(e)); FALSE })
if (!isTRUE(written)) {
  message("grnbench wrapper: could not write ", target)
  quit(save = "no", status = 1L)
}
