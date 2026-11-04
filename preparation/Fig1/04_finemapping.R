# Prepare credible-set counts and lead PIPs from all raw fine-mapping rows.
# Use the November 2 EAS eQTL source explicitly. Do not use its older staged file.
source("preparation/common.R")
source("preparation/Fig1/functions.R")
config <- prep_init()
assays <- c("eQTL", "sQTL", "pQTL")
groups <- c("COMB", "AFR", "AMR", "EUR", "SAS", "EAS", "MID")
base <- "output/eas_fine_mapping_2026-11-02"
s11 <- "Supplemental Figures/01_Supplemental_drafts/Figure1_population_credible_set_counts/v2_with_EAS"
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
    if (include_group) "group"), showProgress = FALSE) %>% as_tibble()
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
stopifnot(nrow(distinct(sets, assay, population)) == 21L,
          !anyDuplicated(sets[c("assay", "population", "cs_id")]))
keys <- c("assay", "population", "cs_id")
cs <- sets %>% select(assay, population, cs_id, cs_size)
lead <- sets %>% select(cs_id, lead_pip, lead_variant, n_tied_lead_variants, assay, population)
saved_cs <- read_tsv(prep_reference("paper", file.path(base, "data/credible_sets.tsv.gz")), show_col_types = FALSE)
saved_lead <- read_tsv(prep_reference("paper", file.path(base, "data/lead_variant_pips.tsv.gz")), show_col_types = FALSE)
stopifnot(nrow(cs) == nrow(saved_cs),
          nrow(anti_join(cs, saved_cs, by = names(cs))) == 0L,
          nrow(anti_join(saved_cs, cs, by = names(cs))) == 0L)
lead_check <- lead %>% full_join(saved_lead, by = keys, suffix = c(".new", ".saved"), relationship = "one-to-one")
stopifnot(!anyNA(lead_check), max(abs(lead_check$lead_pip.new - lead_check$lead_pip.saved)) < 1e-12,
          all(lead_check$lead_variant.new == lead_check$lead_variant.saved),
          all(lead_check$n_tied_lead_variants.new == lead_check$n_tied_lead_variants.saved))
summary_data <- fig1_finemapping_summaries(sets)
traits <- summary_data$traits
curves <- summary_data$curves %>% arrange(assay, k)
stopifnot(max(traits$number_credible_sets) == 10L,
          identical(curves$denominator[curves$k == 1L], c(29088L, 5763L, 2529L)))
saved_curves <- read_tsv(prep_reference("paper",
  "output/figure1_flowchart_draft/science_print_size/three_assay_D/credible_set_count_curves.tsv"), show_col_types = FALSE) %>%
  mutate(assay = factor(assay, levels = assays)) %>% arrange(assay, k)
stopifnot(isTRUE(all.equal(curves, saved_curves, tolerance = 1e-12)))
histogram <- traits %>% count(assay, population, number_credible_sets, name = "number_traits") %>%
  complete(assay = factor(assays, levels = assays), population = factor(groups, levels = groups),
    number_credible_sets = 1:10, fill = list(number_traits = 0L)) %>%
  group_by(assay, population) %>% mutate(total_traits = sum(number_traits),
    percent_traits = 100 * number_traits / total_traits) %>% ungroup()
count_summary <- traits %>% group_by(assay, population) %>% summarise(number_traits = n(),
  total_credible_sets = sum(number_credible_sets), mean_credible_sets = mean(number_credible_sets),
  median_credible_sets = median(number_credible_sets), minimum = min(number_credible_sets),
  maximum = max(number_credible_sets), percent_one_set = 100 * mean(number_credible_sets == 1L), .groups = "drop")
walk(c(file.path(base, "credible_set_counts"), s11), function(folder) {
  prep_write_tsv(traits, "paper", file.path(folder, "credible_sets_per_molecular_trait.tsv.gz"))
  prep_write_tsv(histogram, "paper", file.path(folder, "histogram_counts.tsv"))
  prep_write_tsv(count_summary, "paper", file.path(folder, "population_credible_set_count_summary.tsv"))
})
prep_write_tsv(cs, "paper", file.path(base, "data/credible_sets.tsv.gz"))
prep_write_tsv(lead, "paper", file.path(base, "data/lead_variant_pips.tsv.gz"))
prep_write_tsv(sets %>% select(assay, population, cs_id, cs_size, molecular_trait_id), "paper",
               file.path(base, "data/credible_set_traits.tsv.gz"))
prep_write_tsv(inputs %>% filter(population == "EAS") %>%
  mutate(md5 = unname(tools::md5sum(path))) %>% select(assay, population, path, md5), "paper",
  file.path(base, "data/EAS_source_manifest.tsv"))
prep_write_tsv(traits %>% filter(population != "EAS"), "paper",
  "output/figure1_population_credible_set_counts/credible_sets_per_molecular_trait.tsv.gz")
prep_write_tsv(curves, "paper",
  "output/figure1_flowchart_draft/science_print_size/three_assay_D/credible_set_count_curves.tsv")

# Preserve the input contract used by the current Figure 1E notebook.
combined <- sets %>% filter(population == "COMB") %>% select(cs_id, lead_pip, assay) %>% mutate(group = "COMB")
combined_summary <- summary_data$lead_summary %>% filter(population == "COMB") %>%
  mutate(fraction = n_pip_above_09 / n_cs, label = sprintf("%.1f%%", 100 * fraction))
stopifnot(identical(combined_summary$n_cs, c(143512L, 17674L, 11922L)),
          identical(combined_summary$n_pip_above_09, c(70246L, 8301L, 6284L)))
density_plot <- ggplot(combined, aes(lead_pip, colour = assay)) +
  geom_density(bw = .04, bounds = c(0, 1), n = 2048) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0))
inset <- ggplot(combined_summary, aes(fraction, assay)) + geom_col()
prep_save_rds(list(E = density_plot, inset = inset), "paper",
  "output/figure1_flowchart_draft/bottom_row_density_options/combined_PIP_percent_inset_display.rds")
prep_write_tsv(combined_summary, "paper",
  "output/figure1_flowchart_draft/bottom_row_density_options/combined_lead_PIP_summary.tsv")
message("Figure 1 fine-mapping inputs verified: all ", nrow(cs), " credible sets and lead PIPs match the current data.")
