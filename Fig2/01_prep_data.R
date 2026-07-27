# =============================================================================
# Figure 2 — data prep. Reads raw inputs ONCE, derives every per-panel data
# frame, and writes them to fig2/data/<panel>.rds for the panel scripts.
# Run from project root:  Rscript fig2/01_prep_data.R
# =============================================================================
source("fig2/00_setup.R")

message("Loading QTL inputs ...")
eQTL_COMB_susie <- add_derived(read_qtl("susie_eqtl"))
sQTL_COMB_susie <- add_derived(read_qtl("susie_sqtl"))
# Proteomics is optional: input_path() errors if the file is absent (no gs_uri to
# stage it), so treat that as "omit pQTL" — mirrors fig5's optional `nohom` input.
pqtl_path  <- tryCatch(input_path("susie_pqtl"), error = function(e) NA_character_)
pQTL_susie <- if (!is.na(pqtl_path)) add_derived(fread(pqtl_path, nThread = 4)) else NULL
if (is.null(pQTL_susie)) message("Note: proteomics (pQTL) input not found — omitted.")

saveit <- function(x, name) {
  saveRDS(x, file.path(DATA_DIR, paste0(name, ".rds")))
  message("  wrote ", file.path(DATA_DIR, paste0(name, ".rds")))
}

## ---- Panel A: MAF CDF of causal variants ----------------------------------
assay_tbl <- function(df, assay) df %>% transmute(variant, pip, gvs_max_af, posterior_mean, Assay = assay)
all_causal_variants <- bind_rows(
    if (!is.null(pQTL_susie)) assay_tbl(pQTL_susie, "Proteomics"),
    assay_tbl(eQTL_COMB_susie, "Expression"),
    assay_tbl(sQTL_COMB_susie, "Splicing")
  ) %>%
  filter(pip > 0.9, !is.na(gvs_max_af), gvs_max_af > 0)
xmin_box <- min(log10(all_causal_variants$gvs_max_af), na.rm = TRUE)
x_center <- mean(c(xmin_box, log10(0.01)))
saveit(list(df = all_causal_variants, x_center = x_center), "A")

## ---- Panel B: ancestry-skew, ALL fine-mapped variants (PIP > 0.9) ---------
# Ancestry-skew is variant-level (identical across assays) -> one row per variant.
skew_cols <- c("variant", "pip", "gvs_max_subpop", "gvs_max_af", "odds_ratio", "p_value", "rare")
skew_df <- bind_rows(
    eQTL_COMB_susie %>% select(all_of(skew_cols)),
    sQTL_COMB_susie %>% select(all_of(skew_cols)),
    if (!is.null(pQTL_susie)) pQTL_susie %>% select(all_of(skew_cols))
  ) %>%
  filter(pip > 0.9, !is.na(gvs_max_subpop)) %>%
  distinct(variant, .keep_all = TRUE) %>%
  mutate(padj = p.adjust(p_value, method = "fdr"),
         enriched = odds_ratio > 5 & padj < 0.05,
         pop = pretty_pop(gvs_max_subpop),
         rare_lab = if_else(rare, "Rare", "Common"))

skew_counts <- skew_df %>% filter(enriched) %>% distinct(variant, pop) %>% count(pop)
or_min      <- min(skew_df$odds_ratio[skew_df$odds_ratio > 0], na.rm = TRUE)
pair_line_df <- skew_df %>% distinct(pop) %>%
  mutate(x = as.numeric(factor(pop)) - 0.32, xend = as.numeric(factor(pop)) + 0.32,
         y = or_min * 0.55, yend = or_min * 0.55)
admix_x  <- match("ADMIX", sort(unique(skew_df$pop)))
label_df <- tibble::tibble(rare_lab = c("Common", "Rare"),
                           x = admix_x + c(-0.2, 0.2), y = or_min * 0.32)
saveit(list(skew_df = skew_df, skew_counts = skew_counts, or_min = or_min,
            pair_line_df = pair_line_df, label_df = label_df), "B")

## ---- Panel C: |log2 aFC| by variant consequence (expression) --------------
remove_list <- c(
  "stop_lost", "inframe_deletion", "stop_retained_variant", "inframe_insertion",
  "CTCF.bound", "mature_mirna_variant", "upstream_gene_variant", "intron_variant",
  "intergenic_variant", "FANTOM5", "CTCF.only", "downstream_gene_variant",
  "synonymous_variant", "start_lost", "missense_variant",
  "non_coding_transcript_variant", "non_coding_transcript_exon_variant",
  "nmd_transcript_variant", "coding_sequence_variant", "NMD_transcript_variant")
consequence_order <- c(
  "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
  "splice_acceptor_variant", "splice_donor_variant", "splice_donor_5th_base_variant",
  "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
  "splice_region_variant", "PLS", "pELS", "dELS", "inframe_insertion",
  "inframe_deletion", "missense_variant", "coding_sequence_variant",
  "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
  "upstream_gene_variant", "downstream_gene_variant", "intergenic_variant",
  "CTCF.bound", "CTCF.only", "DNase.H3K4me3", "FANTOM5")

worst_csq <- eQTL_COMB_susie %>%
  filter(pip > 0.9, consequence != ".") %>%
  select(variant, log2_aFC, pip, rare, posterior_mean, consequence, molecular_trait_id,
         FANTOM5, pELS, dELS, PLS, `CTCF.bound`, `DNase.H3K4me3`, `CTCF.only`) %>%
  separate_rows(consequence, sep = ",") %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = consequence, values_from = value, values_fill = FALSE) %>%
  pivot_longer(!c(variant, log2_aFC, pip, rare, posterior_mean, molecular_trait_id)) %>%
  mutate(name = trimws(name)) %>%
  filter(!name %in% remove_list, !str_detect(name, "DNase"), value == TRUE) %>%
  mutate(consequence_rank = match(name, consequence_order)) %>%
  filter(!is.na(consequence_rank)) %>%
  group_by(variant) %>% slice_min(consequence_rank, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(name = if_else(str_detect(name, "splice"), "Splice region", name))

afc_plot_df <- worst_csq %>%
  mutate(abs_log2_aFC = pmin(abs(log2_aFC), 3),
         direction = if_else(log2_aFC >= 0, "Up", "Down")) %>%
  mutate(name = str_replace(name, "_prime_", "' "),
         name = str_replace(name, "UTR_variant", " UTR"),
         name = case_when(name == "stop_gained" ~ "Stop Gain",
                          name == "frameshift_variant" ~ "Frameshift", TRUE ~ name)) %>%
  mutate(name = fct_reorder(name, abs_log2_aFC, .fun = median, .desc = TRUE)) %>%
  mutate(ExpectedHet = (1 + 2^log2_aFC) / 2,
         high_effect = ExpectedHet < 0.6 | ExpectedHet > 1.4,
         color_group = paste0(name, "_", rare))

log2_afc_regression <- afc_plot_df %>%
  group_by(name) %>% do(broom::tidy(lm(abs(log2_aFC) ~ rare, data = .))) %>%
  filter(!str_detect(term, "Int")) %>% ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "fdr"),
         stars = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""),
         label = paste0("+", round(estimate, 2), stars), x = 3.1)

base_names  <- levels(afc_plot_df$name)
base_cols   <- as.character(paletteer::paletteer_d("futurevisions::pegasi")[seq_along(base_names)])
names(base_cols) <- base_names
afc_palette <- c(setNames(colorspace::darken(base_cols, 0.3),  paste0(base_names, "_TRUE")),
                 setNames(colorspace::lighten(base_cols, 0.3), paste0(base_names, "_FALSE")))
saveit(list(afc_plot_df = afc_plot_df, log2_afc_regression = log2_afc_regression,
            afc_palette = afc_palette), "C")

## ---- Panel D: cumulative expected het expression (expression) -------------
cumulative_effects <- eQTL_COMB_susie %>%
  filter(!is.na(log2_aFC)) %>%
  make_cumulative_afc_curve_both(gene_col = "molecular_trait_id")
saveit(cumulative_effects, "D")

## ---- Panel E: QTL-burden outlier enrichment -------------------------------
saveit(read_qtl("burden_enrichment"), "E")

## ---- Panel F: median genes per individual ---------------------------------
median_pc <- read_qtl("median_ko") %>%
  filter(gene_type == "protein_coding") %>%
  mutate(gene_type = "Protein Coding",
         lower = as.numeric(str_extract(as.character(pct_bin), "-?[0-9]+\\.?[0-9]*")))
saveit(median_pc, "F")

## ---- Panel G: validated dosage-like events (dumbbell) ---------------------
outliers_df <- readRDS(input_path("burden_count")) %>%
  mutate(label = str_replace(label, "\\(", "\n"),
         label = str_replace(label, ", ", "\n"),
         label = str_remove(label, "\\)"),
         label = if_else(str_detect(label, "Drug"),
                         str_replace(label, "DrugInteraction", "Drug Target"), label),
         label = fct_reorder(label, total))
saveit(outliers_df, "G")

message("Data prep complete.")
