# Input and output locations shared by all R scripts.
# Run scripts from the repository root, e.g. `Rscript 03_host_response/01_differential_expression.R`.
#   CMT_DATA_DIR    directory holding the input tables (default: data/)
#   CMT_OUTPUT_DIR  directory for results (default: output/)
DATA_DIR <- Sys.getenv("CMT_DATA_DIR", "data")
OUTPUT_DIR <- Sys.getenv("CMT_OUTPUT_DIR", "output")

in_file <- function(...) {
  f <- file.path(DATA_DIR, ...)
  if (!file.exists(f)) stop("Input file not found: ", f, " (see docs/input_files.md)")
  f
}

out_dir <- function(module) {
  d <- file.path(OUTPUT_DIR, module)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
