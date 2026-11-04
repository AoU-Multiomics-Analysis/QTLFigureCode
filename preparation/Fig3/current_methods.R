# Portable entry points and panel summaries. The statistical methods are in
# methods.R and symmetric_methods.R, retained from the current local sources.

# Keep canonical-CDS selection and interval overlap unchanged. Read compressed
# GTF directly rather than invoking macOS-only gzcat through a shell.
.canonical_cds <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    path <- input_path("gene_gtf")
    g <- data.table::fread(path, sep = "\t", header = FALSE, skip = "chr",
                          select = c(1, 3, 4, 5, 9),
                          col.names = c("chrom", "feature", "start", "end", "annotation"), showProgress = FALSE) %>%
      as_tibble() %>% filter(feature == "CDS") %>%
      mutate(gene = str_match(annotation, 'gene_id "([^".]+)')[, 2],
             canonical = str_detect(annotation, fixed('tag "Ensembl_canonical"'))) %>%
      filter(canonical) %>% select(chrom, start, end, gene)
    stopifnot(nrow(g) > 0, all(str_detect(g$gene, "^ENSG[0-9]+$")))
    cache <<- data.table::as.data.table(g)
    data.table::setkey(cache, chrom, start, end)
    cache
  }
})

fig3_gene_counts <- function(gene_results, genes_all_three) {
  gene_results %>% filter(Gene %in% genes_all_three) %>% group_by(Gene) %>%
    summarise(NumColocs = sum(Coloc), .groups = "drop") %>% count(NumColocs) %>% mutate(prop = n / sum(n))
}

fig3_signal_counts <- function(components, genes_all_three) {
  signal_counts <- components %>% filter(!is.na(coloc_component)) %>% distinct(Gene, coloc_component) %>%
    count(Gene) %>% count(n, name = "nn")
  coloc_rate <- components %>% filter(Gene %in% genes_all_three) %>% group_by(type) %>%
    summarise(n_cs = n_distinct(cs_id), n_coloc = n_distinct(cs_id[has_coloc]), coloc_rate = n_coloc / n_cs, .groups = "drop")
  list(signal_counts = signal_counts, coloc_rate = coloc_rate)
}

fig3_interaction_data <- function(model_input) {
  formulas <- list(
    `eQTL-sQTL` = log1p(abs(posterior_mean.y)) ~ log1p(abs(posterior_mean.x)) * log_maf,
    `eQTL-pQTL` = log1p(posterior_mean.y) ~ log1p(posterior_mean.x) * log_maf,
    `sQTL-pQTL` = log1p(abs(posterior_mean.y)) ~ log1p(abs(posterior_mean.x)) * log_maf)
  imap_dfr(model_input, function(d, pair) {
    d <- d %>% filter(!is.na(posterior_mean.x), !is.na(posterior_mean.y))
    broom::tidy(lm(formulas[[pair]], data = d), conf.int = TRUE) %>%
      filter(!str_detect(term, "Int")) %>% mutate(type = pair)
  })
}

fig3_effect_data <- function(model_input) {
  binned <- model_input[["eQTL-pQTL"]] %>% filter(!is.na(max_af)) %>%
    mutate(maf_bin = case_when(max_af >= 0.05 ~ "Common", max_af >= 0.01 ~ "Low-Freq", TRUE ~ "Rare"),
           posterior_mean.x = log1p(posterior_mean.x), posterior_mean.y = log1p(posterior_mean.y))
  slopes <- binned %>% group_by(maf_bin) %>% summarise(n = n(),
    slope = coef(lm(posterior_mean.y ~ posterior_mean.x))[["posterior_mean.x"]], .groups = "drop")
  list(data = binned, slopes = slopes,
       limit = quantile(abs(c(binned$posterior_mean.x, binned$posterior_mean.y)), 0.97, na.rm = TRUE))
}

fig3_posterior_thresholds <- function(gene, cs, variant) {
  inputs <- list(GLCP = gene, RCP = cs, SCP = variant)
  imap_dfr(inputs, function(d, level) map_dfr(c(0.01, 0.1, 0.5, 0.9), function(threshold) {
    d %>% group_by(type) %>% summarise(n_tested = n(), n_coloc = sum(.data[[level]] > threshold, na.rm = TRUE), .groups = "drop") %>%
      mutate(threshold = threshold, posterior = level)
  }))
}
