# %% Read source fine-mapping rows and reduce them to credible sets
#| Each input row describes a variant in a credible set. The maximum PIP defines the lead PIP.
#| The current EAS eQTL file can be gzip CSV with a .tsv suffix. The reader checks its file signature.
assays <- c("eQTL", "sQTL", "pQTL")
groups <- fig1_groups
read_selected <- function(path, include_group = FALSE) {
  # The current EAS eQTL download is gzip CSV with a .tsv filename.
  # Inspect the magic bytes. Copy to a temporary gzip filename when required.
  connection <- file(path, "rb")
  magic <- readBin(connection, "raw", n = 2L)
  close(connection)
  input_path <- path
  if (identical(magic, as.raw(c(31, 139))) && !str_detect(path, "\\.gz$")) {
    input_path <- tempfile(fileext = ".csv.gz")
    stopifnot(file.copy(path, input_path))
    on.exit(unlink(input_path))
  }
  data.table::fread(input_path, select = c("molecular_trait_id", "cs_id", "cs_size", "variant", "pip",
    if (include_group) "group"), showProgress = FALSE, nThread = 2L) %>% as_tibble()
}
inputs <- crossing(assay = assays, population = groups) %>%
  mutate(relative = case_when(
    assay == "eQTL" ~ paste0("Fig1Data/PopulationLevelFinemapping/eQTL/AoU_", population, "_9k_eQTL_SusieMerged.tsv.gz"),
    assay == "sQTL" ~ paste0("Fig1Data/PopulationLevelFinemapping/sQTL/AoU_", population, "_9k_sQTL_SusieMerged.annotated.tsv"),
    TRUE ~ "Fig1Data/PopulationLevelFinemapping/pQTL/pQTL_fm_data_10k.tsv"))
inputs$path <- map2_chr(seq_len(nrow(inputs)), inputs$relative, function(i, relative) {
  if (inputs$assay[[i]] == "eQTL" && inputs$population[[i]] == "EAS") {
    prep_input_file("QTL_FIG1_EAS_EQTL_FILE", "input",
      "Fig1Data/PopulationLevelFinemapping/eQTL/EAS_20261102.eQTL_SusieMerged.tsv")
  } else prep_input("input", relative)
})
protein_rows <- read_selected(inputs$path[inputs$assay == "pQTL"][[1]], include_group = TRUE)
sets <- pmap_dfr(inputs, function(assay, population, relative, path) {
  message("Figure 1: Aggregate raw credible sets: ", assay, " / ", population)
  d <- if (assay == "pQTL") filter(protein_rows, group == population) else read_selected(path, assay == "sQTL")
  if (assay == "sQTL") stopifnot(all(d$group == paste0("AoU_", population, "_9k_sQTL")))
  fig1_credible_sets(d %>% select(molecular_trait_id, cs_id, cs_size, variant, pip), assay, population)
}) %>% mutate(assay = factor(assay, levels = assays), population = factor(population, levels = groups))
stopifnot(nrow(distinct(sets, assay, population)) == length(assays) * length(groups),
          !anyDuplicated(sets[c("assay", "population", "cs_id")]))

rm(protein_rows)
invisible(gc())
summary_data <- fig1_finemapping_summaries(sets)
traits <- summary_data$traits
message("Calculated ", nrow(sets), " credible sets from source variant rows.")
