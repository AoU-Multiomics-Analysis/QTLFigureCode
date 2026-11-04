# Prepare the ACVRL1 locus and PromoterAI windows with the tracked pure functions.
source("preparation/common.R")
source("preparation/Fig5/functions.R")
suppressPackageStartupMessages(library(scales))
config <- prep_init()
functions <- c("normalize_variant_id", "make_variant_id", "fmt_af", "make_freq_label",
  "build_variant_freq", ".acvrl1_qtl", "prepare_promoterAI_locus",
  "build_causal_promoterAI_facets", "signed_square_trans", "plot_causal_promoterAI_facets")
found <- character()
for (expr in parse("Fig5/00_setup.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% functions &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
    eval(expr)
    found <- c(found, as.character(expr[[2]]))
  }
}
stopifnot(setequal(found, functions))
message("Read ACVRL1 eQTL and PromoterAI source rows.")
all_eqtl <- as_tibble(data.table::fread(prep_input("susie",
  "COMB.expression_SusieMerged.annotated.ancestry_skew.tsv.gz"),
  select = c("variant", "gene_name", "position", "pip", "cs_id",
    "gvs_max_af", "gvs_max_subpop", "log2_aFC"))) %>%
  mutate(rare = gvs_max_af < 0.01, chrom = str_extract(variant, "^chr[^_]+"))
eqtl <- all_eqtl %>% filter(gene_name == "ACVRL1")
qtl <- .acvrl1_qtl(eqtl)
# Preserve the original first-record-per-variant lookup before selecting the locus.
frequencies <- build_variant_freq(all_eqtl) %>%
  semi_join(qtl %>% transmute(variant_key = normalize_variant_id(variant)), by = "variant_key")
promoter <- read_tsv(prep_input("source", "Fig5/reference/acvrl1_promoterAI.tsv"),
  show_col_types = FALSE)
locus <- prepare_promoterAI_locus(promoter, qtl, "ACVRL1", "ENST00000419526.6",
  focal_position = 51912316, window_bp = 100)
levels <- unique(qtl$cs_id[!is.na(qtl$cs_id)])
palette <- setNames(rcartocolor::carto_pal(max(3, length(levels)), "Prism")[seq_along(levels)], levels)
facets <- build_causal_promoterAI_facets(locus, palette, flank_bp = 20,
  freq_lookup = frequencies)
expected_plot <- readRDS(prep_reference("source", "Fig5/objects/E.rds"))
expected <- expected_plot$data
columns <- names(expected)
ordered <- function(data) data %>% select(all_of(columns)) %>%
  mutate(across(where(is.factor), as.character)) %>% arrange(facet_id, pos, ref, alt)
stopifnot(isTRUE(all.equal(ordered(facets), ordered(expected),
  check.attributes = FALSE, tolerance = 1e-10)))
predictions <- plot_causal_promoterAI_facets(locus, palette, flank_bp = 20,
  effect_threshold = 0.5, freq_lookup = frequencies,
  motif_annotations = tibble(facet_id = "chr12_51912316_C_T", xmin = -3.5, xmax = 5.5, tf = "SP1"),
  nucleotide_cols = c(A = "#E64B35", C = "#00A087", G = "#F39B7F", T = "#4DBBD5"))
stopifnot(isTRUE(all.equal(ordered(predictions$data), ordered(facets),
  check.attributes = FALSE, tolerance = 1e-10)))
prediction_layers <- function(plot) lapply(ggplot2::ggplot_build(plot)$data,
  function(data) as_tibble(data) %>% mutate(across(where(is.factor), as.character)))
stopifnot(isTRUE(all.equal(prediction_layers(predictions), prediction_layers(expected_plot),
  check.attributes = FALSE, tolerance = 1e-10)))

# Keep the accepted cCRE, gene-track, inset, and axis layout as an explicit renderer
# dependency. Every regional and inset QTL plotting frame is replaced from qtl.
locus_template <- readRDS(prep_reference("source", "Fig5/objects/D.rds"))
locus_plot <- fig5_refresh_locus_template(locus_template, qtl)
prep_save_rds(qtl, "source", "Fig5/data/S31_qtl.rds")
prep_save_rds(frequencies, "source", "Fig5/data/S31_variant_freq.rds")
prep_save_rds(locus, "source", "Fig5/data/S31_promoter_locus.rds")
prep_save_rds(facets, "source", "Fig5/data/S31_promoter_facets.rds")
prep_save_rds(locus_plot, "source", "Fig5/objects/D.rds")
prep_save_rds(predictions, "source", "Fig5/objects/E.rds")
message("S31 D/E objects now use prepared QTL and prediction data; accepted renderer values match.")
