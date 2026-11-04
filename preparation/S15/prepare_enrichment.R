# Fit the annotation-on-PIP models used for the enrichment plot.
source("preparation/common.R")
source("Fig1/enrichment_glm_fit.R")

inputs <- read_tsv("Fig1/inputs.tsv", col_types = cols(.default = col_character()),
                   show_col_types = FALSE) %>% filter(key == "annot")
fits <- map2_dfr(inputs$local_subpath, inputs$Assay, function(path, assay) {
  message("Fit annotation enrichment: ", assay)
  GlmAnnotationEnrichment(prep_input("input", path)) %>% mutate(type = assay)
})
results <- fits %>% filter(term == "pip") %>% ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "fdr"),
         conf.low = estimate - 1.96 * std.error,
         conf.high = estimate + 1.96 * std.error,
         separated = std.error > 10) %>% arrange(desc(estimate))

# Reuse only the pure plotting definitions. Do not run setup or font staging.
definitions <- parse("Fig1/00_setup.R")
for (name in c("PANEL_ANNOTATION_LABELS", "PlotGlmEnrichment", "RenderGlmEnrichment")) {
  selected <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                       identical(x[[2]], as.name(name)), as.list(definitions))
  stopifnot(length(selected) == 1L)
  eval(selected[[1]])
}
definitions <- parse("Fig1/panels/panelH.R")
selected <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                     identical(x[[2]], as.name("variant_functional_groups")), as.list(definitions))
stopifnot(length(selected) == 1L)
eval(selected[[1]])
AssayColors <- c(eQTL = "#0072B2", sQTL = "#009E73", pQTL = "#D55E00")
plot <- PlotGlmEnrichment(results, variant_functional_groups)
if (prep_has_reference("source", "Fig1/objects/H.rds")) {
  expected <- readRDS(prep_reference("source", "Fig1/objects/H.rds"))$data
  fields <- c("type", "name", "estimate", "std.error", "p.value", "padj", "conf.low", "conf.high")
  before <- expected %>% select(all_of(fields)) %>% arrange(type, name)
  after <- plot$data %>% select(all_of(fields)) %>% arrange(type, name)
  stopifnot(isTRUE(all.equal(before, after, tolerance = 1e-7, check.attributes = FALSE)))
}
prep_write_tsv(results, "source", "Fig1/objects/enrichment_glm.tsv")
prep_save_rds(plot, "source", "Fig1/objects/H.rds")
message("Annotation models and enrichment plot inputs prepared.")
