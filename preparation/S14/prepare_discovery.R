# Prepare S14 from source credible-set tables and cohort/reference metadata.
source("preparation/common.R")
source("preparation/S14/functions.R")
suppressPackageStartupMessages({library(arrow); library(ggrepel)})

context <- "output/eqtl_credible_set_discovery/external_context"
current <- "Supplemental Figures/01_Supplemental_drafts/Figure1_eQTL_discovery/v2_with_EAS"
genes <- read_tsv(prep_input("paper", "output/qtl_circos_simplified_options/gene_reference.tsv"),
                  show_col_types = FALSE) %>% select(gene_id, gene_name, gene_type)
stopifnot(!anyDuplicated(genes$gene_id))
groups <- c("COMB", "AFR", "AMR", "EUR", "SAS", "EAS", "MID")

message("Count distinct AoU credible sets from the source fine-mapping tables.")
aou <- map_dfr(groups, function(group) {
  path <- if (group == "EAS") prep_input_file("QTL_FIG1_EAS_EQTL_FILE", "input",
    "Fig1Data/PopulationLevelFinemapping/eQTL/EAS_20261102.eQTL_SusieMerged.tsv") else
      prep_input("input", paste0("Fig1Data/PopulationLevelFinemapping/eQTL/AoU_", group, "_9k_eQTL_SusieMerged.tsv.gz"))
  raw <- discovery_read(path, c("molecular_trait_id", "cs_id"))
  stopifnot(!anyNA(raw))
  discovery_sets(raw, "molecular_trait_id", "AoU", group, genes)
})

message("Count analysis participants without exporting participant records.")
samples <- data.table::fread(prep_input("input", "reference/COMB_eQTL_9k_sample_list.tsv"), select = "research_id") %>% as_tibble()
ancestry <- data.table::fread(prep_input("input", "echo_v4_r2.ancestry_preds.tsv"),
                            select = c("research_id", "ancestry_pred_other")) %>%
  as_tibble() %>% semi_join(samples, by = "research_id")
stopifnot(nrow(samples) == 8895L, nrow(ancestry) == nrow(samples), !anyDuplicated(samples$research_id))
sizes <- ancestry %>% mutate(analysis = str_to_upper(ancestry_pred_other)) %>%
  count(analysis, name = "sample_size") %>% bind_rows(tibble(analysis = "COMB", sample_size = nrow(samples)))
stopifnot(sizes$sample_size[sizes$analysis == "EAS"] == 1059L)

message("Count GTEx credible sets and samples from the 50 source tissue analyses.")
gtex_dir <- file.path(prep_root("paper"), context, "reference/gtex_v10_susie")
gtex_files <- list.files(gtex_dir, "[.]parquet$", full.names = FALSE)
stopifnot(length(gtex_files) == 50L)
gtex <- map(gtex_files, function(name) {
  tissue <- str_remove(name, "\\.v10.*$")
  path <- prep_input("paper", file.path(context, "reference/gtex_v10_susie", name))
  raw <- arrow::read_parquet(path, col_select = c("phenotype_id", "cs_id")) %>% as_tibble()
  sets <- discovery_sets(raw, "phenotype_id", "GTEx v10", tissue, genes)
  covariates <- prep_input("paper", file.path(context, "reference/gtex_v10_covariates", paste0(tissue, ".v10.covariates.txt")))
  header <- str_split(read_lines(covariates, n_max = 1), "\t")[[1]]
  stopifnot(n_distinct(header[-1]) == length(header) - 1L, !anyNA(sets$cs_id))
  list(sets = sets, metadata = tibble(tissue, sample_size = length(header) - 1L))
})
gtex_sets <- bind_rows(map(gtex, "sets"))
gtex_sizes <- bind_rows(map(gtex, "metadata"))

message("Read the INTERVAL and TOPMed source credible sets.")
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
prep_write_tsv(genes %>% filter(gene_id %in% shared), "paper",
  file.path(context, "protein_coding/three_study_shared_genes/shared_genes.tsv"))

for (panel in c("A", "B")) {
  selected <- if (panel == "A") sets else sets %>% filter(protein_coding, gene_id %in% shared)
  counts <- discovery_counts(selected)
  points <- counts %>% mutate(population = if_else(study == "AoU", analysis, study)) %>%
    left_join(sizes %>% select(analysis, sample_size), by = "analysis") %>%
    mutate(sample_size = case_when(study == "INTERVAL" ~ as.integer(interval_meta$sample_size),
                                  study == "TOPMed" ~ 6454L, TRUE ~ sample_size)) %>%
    arrange(match(population, c(groups, "INTERVAL", "TOPMed")))
  stopifnot(nrow(points) == 9L, !anyNA(points$sample_size), all(points$credible_sets > 0))
  target <- readRDS(prep_reference("paper", file.path(current, paste0("source_", panel), "plot.rds")))
  compared <- target$data %>% select(population, sample_size, credible_sets, genes_with_credible_sets) %>%
    left_join(points %>% select(population, sample_size, credible_sets, genes_with_credible_sets),
              by = "population", suffix = c("_saved", "_new"), relationship = "one-to-one")
  for (field in c("sample_size", "credible_sets", "genes_with_credible_sets")) {
    stopifnot(identical(as.numeric(compared[[paste0(field, "_saved")]]), as.numeric(compared[[paste0(field, "_new")]])))
  }
  # Preserve the accepted style columns while replacing every calculated value.
  display <- target$data %>% select(-sample_size, -credible_sets, -genes_with_credible_sets) %>%
    left_join(points %>% select(population, sample_size, credible_sets, genes_with_credible_sets),
              by = "population", relationship = "one-to-one")
  model <- lm(log10(credible_sets) ~ log10(sample_size), data = display)
  curve <- tibble(sample_size = seq(min(display$sample_size), max(display$sample_size), length.out = 1000)) %>%
    mutate(log10_credible_sets = as.numeric(predict(model, newdata = pick(everything()))), credible_sets = 10^log10_credible_sets)
  stopifnot(isTRUE(all.equal(curve, target$curve, tolerance = 1e-10, check.attributes = FALSE)))
  gtex_selected <- if (panel == "A") gtex_sets else gtex_sets %>% filter(protein_coding)
  tissue_counts <- discovery_counts(gtex_selected) %>% rename(tissue = analysis) %>%
    left_join(gtex_sizes, by = "tissue", relationship = "one-to-one")
  check <- target$gtex %>% select(tissue, credible_sets, sample_size) %>%
    inner_join(tissue_counts %>% select(tissue, credible_sets, sample_size), by = "tissue", suffix = c("_saved", "_new"))
  stopifnot(nrow(check) == 50L, all(check$credible_sets_saved == check$credible_sets_new),
            all(check$sample_size_saved == check$sample_size_new))
  gtex_display <- target$gtex %>% select(-credible_sets, -genes_with_credible_sets, -sample_size) %>%
    left_join(tissue_counts %>% select(tissue, credible_sets, genes_with_credible_sets, sample_size), by = "tissue")
  envelope <- discovery_envelope(gtex_display)
  stopifnot(isTRUE(all.equal(envelope$boundary, target$envelope$boundary, tolerance = 1e-10)))
  plot <- target$plot
  plot$data <- display
  plot$layers[[1]]$data <- envelope$boundary
  plot$layers[[2]]$data <- curve
  prepared <- list(plot = plot, data = display, curve = curve, model = model, gtex = gtex_display, envelope = envelope)
  prep_save_rds(prepared, "paper", file.path(current, paste0("source_", panel), "plot.rds"))
  prep_write_tsv(display, "paper", file.path(current, paste0("panel_", panel, "_plotted_points.tsv")))
  prep_write_tsv(curve, "paper", file.path(current, paste0("panel_", panel, "_regression_curve.tsv")))
  prep_write_tsv(gtex_display, "paper", file.path(current, paste0("panel_", panel, "_gtex_v10_tissue_summary.tsv")))
}
message("S14 source counts, shared genes, nine-point curves, and GTEx envelopes match the accepted figure.")
