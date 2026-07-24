# Runs every tests/test_*.R in a fresh child process. From the project root:
#   Rscript tests/run_all.R
files <- sort(list.files("tests", pattern = "^test_.*\\.R$", full.names = TRUE))
failed <- character(0)
for (f in files) {
  cat("\n==", f, "==\n")
  if (system2("Rscript", f) != 0L) failed <- c(failed, f)
}
cat("\n")
if (length(failed)) {
  cat("FAILED:\n"); cat(paste0("  ", failed, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("All test files passed.\n")
