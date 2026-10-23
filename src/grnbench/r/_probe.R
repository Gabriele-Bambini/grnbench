# grnbench R bridge probe (grnbench.rbridge.probe): Rscript --vanilla _probe.R <package> ...
# Prints ONE JSON line: R's version, its library folders, the installed version of each package named on the
# command line (null when absent) and the R version each was built for (from its DESCRIPTION "Built" field).
# Written by hand, so it works in an R without jsonlite. Reads DESCRIPTION files only: nothing is loaded.
args <- commandArgs(trailingOnly = TRUE)
quote_json <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  x <- gsub("\r", "\\r", x, fixed = TRUE)
  x <- gsub("\t", "\\t", x, fixed = TRUE)
  paste0("\"", x, "\"")
}
or_null <- function(v) if (length(v) != 1L || is.na(v)) "null" else quote_json(v)
installed_version <- function(p) {
  tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
}
built_for <- function(p) {
  b <- tryCatch(suppressWarnings(utils::packageDescription(p, fields = "Built")), error = function(e) NA_character_)
  if (length(b) != 1L || is.na(b)) return(NA_character_)
  m <- regmatches(b, regexpr("[0-9]+\\.[0-9]+(\\.[0-9]+)?", b))   # "R 4.5.1; x86_64-...; date; os" -> "4.5.1"
  if (length(m)) m else NA_character_
}
entries <- function(f) {
  paste(vapply(args, function(p) paste0(quote_json(p), ": ", or_null(f(p))), character(1)), collapse = ", ")
}
cat(paste0("{\"version\": ", quote_json(paste(R.version$major, R.version$minor, sep = ".")),
           ", \"version_string\": ", quote_json(R.version.string),
           ", \"library\": [", paste(vapply(.libPaths(), quote_json, character(1)), collapse = ", "), "]",
           ", \"packages\": {", entries(installed_version), "}",
           ", \"built\": {", entries(built_for), "}}\n"))
