# Minimal assertion helpers shared by the test scripts.
# Each test script sources this, then calls ok() / expect_error().
ok <- function(label, condition) {
  if (!isTRUE(condition)) stop("FAIL: ", label, call. = FALSE)
  cat("  ok -", label, "\n")
}

expect_error <- function(label, expr) {
  ok(label, inherits(try(expr, silent = TRUE), "try-error"))
}
