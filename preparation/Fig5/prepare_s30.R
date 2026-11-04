# Partial recovery: compute ranks from the original supplied aggregate counts.
# No genome-wide ranking pipeline or raw association input was saved locally.
source("preparation/common.R")
config <- prep_init()
path <- prep_reference("paper", "plot_compact_percentiles_four_genes.R")
expressions <- parse(path)
is_panels <- vapply(expressions, function(expr) {
  is.call(expr) && identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("panels"))
}, logical(1))
stopifnot(sum(is_panels) == 1L)
expression <- expressions[[which(is_panels)]]
# Evaluate only the recorded tibble/rank expression, without the plot/export code.
allowed <- c("<-", "tibble", "factor", "c", "*", "/", "+", "(", "mutate",
  "if_else", "==", "-", "sprintf", "|>")
calls <- function(expr) {
  if (!is.call(expr)) return(character())
  c(as.character(expr[[1]]), unlist(lapply(as.list(expr)[-1], calls), use.names = FALSE))
}
stopifnot(all(calls(expression) %in% allowed))
environment <- new.env(parent = globalenv())
eval(expression, envir = environment)
panels <- environment$panels
expected <- read_tsv(prep_reference("paper",
  "Supplemental Figures/01_Supplemental_drafts/Figure5_trait_prioritization/trait_prioritization_data.tsv"),
  show_col_types = FALSE)
stopifnot(isTRUE(all.equal(panels %>% mutate(gene = as.character(gene)), expected,
  check.attributes = FALSE, tolerance = 1e-10)))
prep_write_tsv(panels, "paper",
  "Supplemental Figures/01_Supplemental_drafts/Figure5_trait_prioritization/trait_prioritization_data.tsv")
message("S30 recovered from the saved count/formula expression; upstream genome-wide ranks remain unavailable.")
