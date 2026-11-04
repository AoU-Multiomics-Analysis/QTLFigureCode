# Current panel and supplemental data transformations. Preserve the original
# analysis boundary: saved colocalization, TWAS, annotation, and SuSiE results.

# The local-only resolver supplies the first parquet file selected by the
# original directory reader. No ClinGen file is downloaded or modified.
gene_sets <- function(refresh = FALSE) {
  build_or_cache("gene_sets", refresh = refresh, {
    clingen <- arrow::read_parquet(input_path("clingen")) %>%
      select(targetFromSourceId, allelicRequirements, confidence, score, studyId, targetId, diseaseFromSource)
    derive_gene_sets(clingen, read_any(input_path("drug_interactions")))
  })
}

fig4_r2_data <- function(data) {
  data %>% mutate(H2_upper = SE * 1.96 + H2, adj_R2 = r.squared / H2) %>%
    filter(Pval < 1e-5) %>% select(phenotype_id, AF_threshold, adj_R2, H2, Pval) %>%
    pivot_wider(names_from = AF_threshold, values_from = adj_R2) %>% mutate(difference = `0` - `0.01`)
}

# These are the exact current labelled-panel filters, after gene counting and
# the lenient union-rate join. Keep their source order within each category.
fig4_current_trait_data <- function(data) {
  data %>% filter(trait_category != "neuropsychiatric", trait_category != "renal_urological", count > 5) %>%
    group_by(trait_category) %>% arrange(desc(count)) %>% filter(row_number() < 11) %>% ungroup()
}

fig4_enrichment_data <- function(gene_data) {
  full <- tidyr::expand_grid(outcome = names(ENRICH_OUTCOMES), model = names(ENRICH_PREDICTORS)) %>%
    mutate(res = map2(outcome, model, ~fit_enrich(.x, ENRICH_PREDICTORS[[.y]], gene_data))) %>%
    unnest(res) %>% mutate(type = ENRICH_OUTCOMES[outcome]) %>% relocate(type, model)
  plot <- full %>% filter(term %in% c("max_rare", "max_common")) %>%
    select(type, model, term, OR = estimate, conf.low, conf.high, p.value) %>%
    mutate(term = str_remove_all(term, "max_|ever_|TRUE"), term = recode(term, rare = "Rare", common = "Common"),
           term = factor(term, c("Rare", "Common")),
           type = factor(type, c("Mendelian Gene\n(Any)", "Mendelian Gene\n(Immune/Hematological)", "Druggable gene", "Approved Drug Target")))
  list(all_fits = full, plot_data = plot)
}

fig4_upset_counts <- function(coloc) {
  qtl_sets <- c("eQTL", "pQTL", "sQTL")
  selected <- coloc %>% filter(layer %in% qtl_sets) %>% distinct(trait, gene, layer)
  membership <- selected %>% mutate(value = 1L) %>%
    pivot_wider(names_from = layer, values_from = value, values_fill = 0L)
  for (assay in setdiff(qtl_sets, names(membership))) membership[[assay]] <- 0L
  membership <- membership %>% mutate(qtl_combination = pmap_chr(pick(all_of(qtl_sets)), function(...) paste(qtl_sets[as.logical(c(...))], collapse = "+")))
  intersections <- membership %>% count(qtl_combination, name = "gene_trait_pairs") %>% arrange(desc(gene_trait_pairs), qtl_combination)
  list(membership = membership, intersections = intersections,
       multilayer_n = sum(rowSums(membership[qtl_sets]) >= 2),
       all_three_n = sum(rowSums(membership[qtl_sets]) == 3))
}

fig4_cd40_data <- function(signal, source_variants, gene_id = "ENSG00000101017") {
  columns <- c("study", "trait", "layer", "gwas_cs", "qtl_sig", "gene", "RCP", "CLPP", "RCP_pass_FDR")
  signal <- signal %>% filter(gene == gene_id, layer %in% c("eQTL", "pQTL", "sQTL"), !is.na(CLPP)) %>%
    select(all_of(columns)) %>% mutate(qtl_cs = sub("^.*:", "", qtl_sig), locus_key = paste(study, gwas_cs, sep = "|"))
  stopifnot(nrow(signal) > 0)
  source_variants <- source_variants %>% semi_join(signal, by = c("layer", "qtl_cs"))
  # Base merge order is retained: the original factor order and shape order use
  # this table order. The many-to-many expansion across trait/locus is intended.
  data <- merge(source_variants, signal, by = c("layer", "qtl_cs"), all = FALSE)
  trait_order <- c("systemic lupus erythematosus", "multiple sclerosis", "rheumatoid arthritis",
                   "inflammatory bowel disease", "crohn's disease", "ldl cholesterol", "monocyte")
  trait_order <- trait_order[trait_order %in% unique(data$trait)]
  data$trait <- factor(data$trait, trait_order)
  loci <- data %>% distinct(trait, locus_key) %>% arrange(trait, locus_key) %>% group_by(trait) %>%
    mutate(locus_rank = row_number(), n_loci = n(), row_offset = (locus_rank - (n_loci + 1) / 2) * 0.13) %>% ungroup()
  data <- merge(data, loci, by = c("trait", "locus_key"), all.x = TRUE)
  data %>% mutate(trait_y = length(trait_order) - as.integer(trait) + 1 + row_offset,
                  assay = factor(layer, c("eQTL", "pQTL", "sQTL")),
                  credible_set = factor(paste(layer, sub(".*_", "", qtl_cs)), levels = unique(paste(layer, sub(".*_", "", qtl_cs)))),
                  position_mb = position / 1e6)
}
