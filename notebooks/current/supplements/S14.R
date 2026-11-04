# %% Define source-table transformations
#| Gene types come from GENCODE v48. Credible sets are distinct trait and set-ID pairs. The GTEx ellipse is a display envelope with 2% padding; it is not a confidence interval.
# @include notebooks/current/analysis/S14_helpers.R

# %% Read the gene reference and count AoU credible sets
#| Read the original GTF and the seven source eQTL fine-mapping tables. Remove gene-version suffixes before the gene-type join. The EAS override selects the current source result.
genes <- s14_read_gene_gtf(prep_input("input", "gencode.v48.annotation.gtf.gz")) %>%
  select(gene_id, gene_name, gene_type)
groups <- c("COMB", "AFR", "AMR", "EUR", "SAS", "EAS", "MID")
aou <- map_dfr(groups, function(group) {
  path <- if (group == "EAS") prep_input_file("QTL_FIG1_EAS_EQTL_FILE", "input",
    "Fig1Data/PopulationLevelFinemapping/eQTL/EAS_20261102.eQTL_SusieMerged.tsv") else
      prep_input("input", paste0("Fig1Data/PopulationLevelFinemapping/eQTL/AoU_", group, "_9k_eQTL_SusieMerged.tsv.gz"))
  raw <- discovery_read(path, c("molecular_trait_id", "cs_id"))
  stopifnot(!anyNA(raw))
  discovery_sets(raw, "molecular_trait_id", "AoU", group, genes)
})
message("Counted credible sets in seven AoU source analyses.")

# %% Count the participants in each analysis
#| Restrict ancestry assignments to the combined eQTL sample list. Keep participant records in memory and export only group totals.
samples <- data.table::fread(prep_input("input", "reference/COMB_eQTL_9k_sample_list.tsv"),
  select = "research_id") %>% as_tibble()
ancestry <- data.table::fread(prep_input("input", "echo_v4_r2.ancestry_preds.tsv"),
  select = c("research_id", "ancestry_pred_other")) %>% as_tibble() %>% semi_join(samples, by = "research_id")
stopifnot(nrow(samples) == 8895L, nrow(ancestry) == nrow(samples), !anyDuplicated(samples$research_id))
sizes <- ancestry %>% mutate(analysis = str_to_upper(ancestry_pred_other)) %>%
  count(analysis, name = "sample_size") %>%
  bind_rows(tibble(analysis = "COMB", sample_size = nrow(samples)))
stopifnot(sizes$sample_size[sizes$analysis == "EAS"] == 1059L)
rm(samples, ancestry)
message("Counted participants by analysis; participant records were removed from memory.")

# %% Count credible sets and samples in 50 GTEx tissues
#| Each tissue count comes from its source credible-set parquet. The covariate header gives its sample count. Panel B uses all protein-coding GTEx genes; the three-study shared-gene restriction applies to the nine study points.
context <- "output/eqtl_credible_set_discovery/external_context"
gtex_dir <- file.path(prep_root("paper"), context, "reference/gtex_v10_susie")
gtex_files <- sort(list.files(gtex_dir, "[.]parquet$", full.names = FALSE))
stopifnot(length(gtex_files) == 50L)
gtex_counts <- map_dfr(gtex_files, function(name) {
  tissue <- str_remove(name, "\\.v10.*$")
  path <- prep_input("paper", file.path(context, "reference/gtex_v10_susie", name))
  raw <- arrow::read_parquet(path, col_select = c("phenotype_id", "cs_id")) %>% as_tibble()
  sets <- discovery_sets(raw, "phenotype_id", "GTEx v10", tissue, genes)
  covariates <- prep_input("paper", file.path(context, "reference/gtex_v10_covariates", paste0(tissue, ".v10.covariates.txt")))
  header <- str_split(read_lines(covariates, n_max = 1), "\t")[[1]]
  stopifnot(n_distinct(header[-1]) == length(header) - 1L, !anyNA(sets$cs_id))
  bind_rows(discovery_counts(sets) %>% mutate(panel = "A"),
    discovery_counts(filter(sets, protein_coding)) %>% mutate(panel = "B")) %>%
    rename(tissue = analysis) %>% mutate(sample_size = length(header) - 1L)
})
message("Counted GTEx credible sets and samples for 50 source tissue analyses.")

# %% Read INTERVAL and TOPMed and find shared protein-coding genes
#| INTERVAL sample size comes from its dataset metadata. TOPMed sample size is 6,454, as reported in its source publication. The shared set contains protein-coding genes with a credible set in AoU combined, INTERVAL and TOPMed.
interval_raw <- arrow::read_parquet(prep_input("paper", file.path(context, "reference/QTD001000.credible_sets.parquet")),
  col_select = c("molecular_trait_id", "cs_id")) %>% as_tibble()
interval_sets <- discovery_sets(interval_raw, "molecular_trait_id", "INTERVAL", "QTD001000", genes)
interval_meta <- read_tsv(prep_input("paper", file.path(context, "reference/INTERVAL_RNA_dataset_metadata.tsv")),
  show_col_types = FALSE) %>% filter(dataset_id == "QTD001000", quant_method == "ge")
stopifnot(nrow(interval_meta) == 1L, interval_meta$sample_size == 4729L)
topmed_file <- prep_input_file("QTL_TOPMED_CS_FILE", "input", "Fig1Data/external/Whole_blood.maf0001.cs.txt")
topmed_raw <- data.table::fread(topmed_file, select = c("phenotype_id", "cs_id")) %>% as_tibble()
topmed_sets <- discovery_sets(topmed_raw, "phenotype_id", "TOPMed", "Whole_blood.maf0001.cs.txt", genes)
stopifnot(nrow(interval_sets) == 43704L, nrow(topmed_sets) == 76545L)
publication <- read_file(prep_input("paper", file.path(context, "reference/TOPMed_full_text.xml")))
stopifnot(str_detect(publication, "6,454"), str_detect(publication, "76,545"))
sets <- bind_rows(aou, interval_sets, topmed_sets)
shared <- shared_discovery_genes(sets)
stopifnot(length(shared) == 10881L)
rm(interval_raw, topmed_raw)
message("Found 10,881 protein-coding genes shared across the three studies.")

# %% Fit the nine-point curves and calculate the GTEx envelopes
#| Panel A counts all genes. Panel B restricts the nine study points to the shared protein-coding genes. Fit log10 credible-set count against log10 sample size. The 50 GTEx tissues do not enter these descriptive fits.
study_style <- tibble(
  population = c(groups, "INTERVAL", "TOPMed"),
  label = c("AoU combined", "AFR", "AMR", "EUR", "SAS", "EAS", "MID", "INTERVAL", "TOPMed"),
  colour = c("#CF597E", "#C8D892", "#7FC385", "#30AC86", "#EBCE8D", "#ECA676", "#E47D72", "#C1C8CD", "#C1C8CD"),
  marker = c(23L, 21L, 21L, 21L, 21L, 21L, 21L, 24L, 22L))
discovery_panels <- set_names(map(c("A", "B"), function(panel) {
  selected <- if (panel == "A") sets else filter(sets, protein_coding, gene_id %in% shared)
  points <- discovery_counts(selected) %>%
    mutate(population = if_else(study == "AoU", analysis, study)) %>%
    left_join(sizes, by = "analysis", relationship = "many-to-one") %>%
    mutate(sample_size = case_when(study == "INTERVAL" ~ as.integer(interval_meta$sample_size),
      study == "TOPMed" ~ 6454L, TRUE ~ sample_size)) %>%
    left_join(study_style, by = "population", relationship = "one-to-one") %>%
    arrange(match(population, study_style$population))
  stopifnot(nrow(points) == 9L, !anyNA(points$sample_size), all(points$credible_sets > 0))
  model <- lm(log10(credible_sets) ~ log10(sample_size), data = points)
  curve <- tibble(sample_size = seq(min(points$sample_size), max(points$sample_size), length.out = 1000)) %>%
    mutate(log10_credible_sets = as.numeric(predict(model, newdata = pick(everything()))),
      credible_sets = 10^log10_credible_sets)
  tissues <- gtex_counts %>% filter(.data$panel == .env$panel) %>% select(-panel)
  envelope <- discovery_envelope(tissues)
  stopifnot(nrow(tissues) == 50L, envelope$details$all_tissues_enclosed)
  list(data = points, model = model, curve = curve, gtex = tissues, envelope = envelope)
}), c("A", "B"))
message("Calculated both nine-point curves and both GTEx display envelopes.")

# %% Draw the discovery panels
#| Diamonds show AoU combined. Circles show the six ancestry analyses. Triangles and squares show the two external blood studies. The shaded envelope contains all 50 GTEx tissues.
# @include notebooks/current/plots/S14_discovery.R
figure <- discovery_plot(discovery_panels$A, "eQTL credible sets (all genes)") +
  discovery_plot(discovery_panels$B, "eQTL credible sets (shared protein-coding genes)") +
  plot_layout(nrow = 1) + plot_annotation(tag_levels = "A")
print(figure)

# %% Save aggregate outputs
#| Save the display, study totals, fitted curves and GTEx tissue totals. No participant identifiers are exported.
save_plot(figure, "supplements/S14", 184, 85)
readr::write_tsv(imap_dfr(discovery_panels, ~mutate(.x$data, panel = .y)),
  output_file("supplements/S14_study_summaries.tsv"))
readr::write_tsv(imap_dfr(discovery_panels, ~mutate(.x$curve, panel = .y)),
  output_file("supplements/S14_regression_curves.tsv"))
readr::write_tsv(gtex_counts, output_file("supplements/S14_gtex_tissue_summaries.tsv"))
