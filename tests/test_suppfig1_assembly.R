# Asserts the assembled Supplemental Figure 1 carries exactly six tags A-F, by
# rendering and inspecting the tag grobs. Panels A and D are 2-plot patchworks,
# so six leaves are lettered -- intended (see the design spec). Panel B's real
# builder needs a live enrichR call, so a stand-in B object is used here; the
# real Panel B is exercised on the Workbench. Run from the project root.
source("tests/fixtures.R")
source("tests/helpers.R")
fx <- make_suppfig1_fixture()
Sys.setenv(SUPPFIG1_INPUT = fx$root,
           SUPPFIG1_DATA   = file.path(fx$root, "data"),
           SUPPFIG1_OBJ    = file.path(fx$root, "objects"),
           SUPPFIG1_OUTDIR = fx$root)
source("suppfig1/00_setup.R")

# Reuse the proven extractor from tests/test_fig3_assembly.R.
extract_tag_labels <- function(patchwork_obj) {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  grid::grid.newpage(); print(patchwork_obj)
  grabbed <- grid::grid.grab()
  labels <- character(0)
  walk <- function(x) {
    if (is.null(x)) return(invisible())
    if (inherits(x, "text")) {
      lab <- paste(x$label, collapse = "")
      if (nchar(lab) == 1L && lab %in% LETTERS) labels <<- c(labels, lab)
    } else if (inherits(x, "gtable")) { for (gr in x$grobs) walk(gr) }
    else if (!is.null(x$children)) { for (ch in x$children) walk(ch) }
  }
  walk(grabbed); labels
}

# Build merged TPMs from the fixture via the pure join (no GTF/GCT read needed).
ann    <- fx$annotations
means  <- function(p) compute_mean_TPM(as_tibble(data.table::fread(p)), ann)
merged <- join_mean_tpms(means(fx$aou_tpms), means(fx$gtex_tpms))

# Panels A, C, D are the real panel plots; B is a stand-in (real B needs enrichR).
saveRDS(merged %>% filter(gene_type == "protein_coding") %>%
          plot_gene_comparison(gene_type = "coding", hex = TRUE),
        file.path(OBJ_DIR, "A.rds"))
saveRDS(ggplot() + theme_void(base_family = FONT) +
          labs(title = "GO enrichment (stand-in)"),
        file.path(OBJ_DIR, "B.rds"))
saveRDS(tf_expression_table(merged, c("TFGENE")) %>%
          select(Description, contains("TPM")) %>% pivot_longer(!Description) %>%
          ggplot(aes(name, log2(value + 1))) + geom_boxplot(),
        file.path(OBJ_DIR, "C.rds"))
saveRDS(merged %>% filter(gene_type == "lncRNA") %>%
          plot_gene_comparison(gene_type = "lncRNA", hex = FALSE),
        file.path(OBJ_DIR, "D.rds"))

source("suppfig1/assemble.R")

tags <- extract_tag_labels(SuppFig1)
ok("Supplemental Figure 1 carries exactly 6 tags", length(tags) == 6L)
ok("its tags are exactly A, B, C, D, E, F",
   identical(sort(unique(tags)), c("A", "B", "C", "D", "E", "F")))
ok("the figure was written", file.exists(file.path(fx$root, "SupplementalFig1.pdf")))

cat("test_suppfig1_assembly.R: all assertions passed\n")
