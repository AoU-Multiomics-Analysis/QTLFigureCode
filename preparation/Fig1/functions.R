# Pure transformations used by the current Figure 1 preparation steps.
# All identifiers are stable gene IDs. Source q values are not recalculated.
fig1_assert_local_input <- function(path) {
  if (grepl("^(gs|s3|https?)://", path)) stop("Input localization error: unresolved cloud URI.", call. = FALSE)
  if (!file.exists(path) || file.access(path, 4) != 0) stop("Missing or unreadable local input: ", path, call. = FALSE)
  invisible(path)
}

fig1_read_gene_gtf <- function(path) {
  fig1_assert_local_input(path)
  connection <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  chunks <- list()
  repeat {
    lines <- readLines(connection, n = 50000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[str_detect(lines, "^[^#][^\\t]*\\t[^\\t]*\\tgene\\t")]
    if (length(lines)) chunks[[length(chunks) + 1L]] <- read_tsv(I(paste(lines, collapse = "\n")),
      col_names = c("chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"),
      col_types = cols(.default = col_character()), quote = "", progress = FALSE)
  }
  genes <- bind_rows(chunks) %>% transmute(
    gene_id = str_extract(str_match(attributes, 'gene_id "([^"]+)"')[, 2], "ENSG[0-9]+"),
    gene_name = str_match(attributes, 'gene_name "([^"]+)"')[, 2],
    gene_type = str_match(attributes, 'gene_type "([^"]+)"')[, 2],
    chromosome, start = as.integer(start), end = as.integer(end), strand)
  stopifnot(nrow(genes) > 0, !anyNA(genes), !anyDuplicated(genes$gene_id))
  genes
}

fig1_ukbb_target_map <- function(gene_ref, hgnc, panel) {
  hgnc <- hgnc %>% filter(status == "Approved", str_detect(ensembl_gene_id, "^ENSG[0-9]+$"))
  approved <- hgnc %>% transmute(gene_name = symbol, gene_id = ensembl_gene_id,
    approved_symbol = symbol, mapping_source = "HGNC approved symbol", priority = 1L)
  previous <- hgnc %>% separate_longer_delim(prev_symbol, delim = "|") %>%
    filter(!is.na(prev_symbol), prev_symbol != "") %>%
    transmute(gene_name = prev_symbol, gene_id = ensembl_gene_id, approved_symbol = symbol,
      mapping_source = "HGNC previous symbol", priority = 2L)
  aliases <- hgnc %>% separate_longer_delim(alias_symbol, delim = "|") %>%
    filter(!is.na(alias_symbol), alias_symbol != "") %>%
    transmute(gene_name = alias_symbol, gene_id = ensembl_gene_id, approved_symbol = symbol,
      mapping_source = "HGNC alias", priority = 4L)
  gencode <- gene_ref %>% filter(gene_type == "protein_coding") %>%
    transmute(gene_name, gene_id, approved_symbol = gene_name,
      mapping_source = "GENCODE v48 exact protein symbol", priority = 0L)
  best <- bind_rows(approved, previous, gencode, aliases) %>% distinct() %>%
    group_by(gene_name) %>% filter(priority == min(priority)) %>% ungroup()
  ambiguous <- best %>% group_by(gene_name) %>%
    filter(n_distinct(gene_id) > 1L, first(priority) != 0L) %>% ungroup()
  mapping <- best %>% anti_join(ambiguous %>% distinct(gene_name), by = "gene_name") %>%
    arrange(gene_name, priority) %>% distinct(gene_name, gene_id, .keep_all = TRUE)
  protein_names <- hgnc %>% select(symbol, ensembl_gene_id, uniprot_ids) %>%
    separate_longer_delim(uniprot_ids, delim = "|") %>% filter(!is.na(uniprot_ids)) %>%
    group_by(uniprot_ids) %>% filter(n_distinct(symbol) == 1L) %>% ungroup()
  protein_identity <- panel %>% inner_join(protein_names,
    by = c("uniprot_id" = "uniprot_ids"), relationship = "many-to-many") %>%
    group_by(gene_name) %>% filter(n_distinct(symbol) == 1L) %>% ungroup() %>%
    left_join(gencode %>% select(current_symbol = gene_name, current_gene_id = gene_id),
      by = c("symbol" = "current_symbol"), relationship = "many-to-many") %>%
    transmute(gene_name, gene_id = coalesce(current_gene_id, ensembl_gene_id),
      approved_symbol = symbol,
      mapping_source = "UKBB ST3 UniProt identity; GENCODE protein ID where available", priority = -1L) %>% distinct()
  mapping <- bind_rows(mapping %>% anti_join(protein_identity %>% distinct(gene_name), by = "gene_name"),
                        protein_identity)
  list(map = mapping, ambiguous = ambiguous, protein_identity = protein_identity)
}

fig1_ukbb_cis_targets <- function(st9, mapping) {
  targets <- st9 %>% filter(`cis/trans` == "cis") %>%
    transmute(assay_target = `Assay Target`, reference_symbol = `cis gene`,
      variant_annotation_gene_id = `Ensembl gene ID`) %>%
    separate_longer_delim(reference_symbol, delim = ",") %>%
    mutate(reference_symbol = str_trim(reference_symbol)) %>%
    left_join(mapping %>% distinct(gene_name, gene_id),
      by = c("reference_symbol" = "gene_name"), relationship = "many-to-many") %>% distinct()
  list(mapped = filter(targets, !is.na(gene_id)), unresolved = filter(targets, is.na(gene_id)), audit = targets)
}

fig1_classify_novelty <- function(raw, genes, known, measured_ids = character()) {
  stopifnot(!anyDuplicated(raw[c("assay", "phenotype_id")]), !anyDuplicated(genes$gene_id))
  x <- raw %>% mutate(gene_id = str_extract(phenotype_id, "ENSG[0-9]+")) %>%
    left_join(genes, by = "gene_id", relationship = "many-to-one") %>%
    left_join(known %>% group_by(assay, gene_id) %>%
      summarise(prior_studies = paste(sort(unique(study)), collapse = ";"), .groups = "drop"),
      by = c("assay", "gene_id"), relationship = "many-to-one") %>%
    mutate(significant = as.double(source_qval) < .05,
      paper_biotype_retained = !is.na(gene_type) & !str_detect(gene_type, "pseudo|IG|TEC|TR"),
      pQTL_lncRNA_excluded_by_paper = assay == "pQTL" & gene_type == "lncRNA",
      previously_reported = !is.na(prior_studies),
      novel_in_paper = significant & !previously_reported & paper_biotype_retained &
        !pQTL_lncRNA_excluded_by_paper,
      novel_for_figure = novel_in_paper & gene_type != "artifact",
      group = case_when(gene_type == "protein_coding" ~ "Protein-Coding",
                        gene_type == "lncRNA" ~ "lncRNA", TRUE ~ "non-coding RNA"),
      type = recode(assay, eQTL = "eGene", sQTL = "sGene", pQTL = "pGene"),
      ukbb_measured = if_else(assay == "pQTL", gene_id %in% measured_ids, NA),
      novelty_context = case_when(!novel_for_figure ~ "Not selected as novel",
        assay != "pQTL" ~ "Absent from the specified significant blood-reference gene sets",
        ukbb_measured ~ "Measured by UKBB; no significant cis association in ST9",
        TRUE ~ "Not represented by a mapped UKBB ST3 protein target"))
  stopifnot(!anyNA(x$gene_id), !anyNA(x$significant), !anyNA(x$novel_for_figure))
  x
}

fig1_credible_sets <- function(variants, assay_name, population_name) {
  stopifnot(nrow(variants) > 0, !anyNA(variants), all(is.finite(variants$pip)),
    all(between(variants$pip, 0, 1)), all(variants$cs_size >= 1),
    all(variants$cs_size == as.integer(variants$cs_size)),
    all(str_remove(variants$cs_id, "_L[0-9]+$") == variants$molecular_trait_id))
  consistency <- variants %>% group_by(cs_id) %>%
    summarise(n_traits = n_distinct(molecular_trait_id), n_sizes = n_distinct(cs_size), .groups = "drop")
  stopifnot(all(consistency$n_traits == 1L), all(consistency$n_sizes == 1L))
  variants %>% group_by(cs_id) %>%
    summarise(molecular_trait_id = first(molecular_trait_id), cs_size = first(cs_size),
      lead_pip = max(pip), lead_variant = first(variant[pip == max(pip)]),
      n_tied_lead_variants = n_distinct(variant[pip == max(pip)]), .groups = "drop") %>%
    mutate(assay = assay_name, population = population_name)
}

fig1_finemapping_summaries <- function(sets) {
  traits <- sets %>% distinct(assay, population, cs_id, molecular_trait_id) %>%
    count(assay, population, molecular_trait_id, name = "number_credible_sets")
  lead_summary <- sets %>% group_by(assay, population) %>% summarise(n_cs = n(),
    median_lead_pip = median(lead_pip), mean_lead_pip = mean(lead_pip),
    n_pip_above_09 = sum(lead_pip > .9), pct_pip_above_09 = 100 * mean(lead_pip > .9), .groups = "drop")
  curves <- map_dfr(1:5, function(k) traits %>% filter(population == "COMB") %>% group_by(assay) %>%
    summarise(k = k, numerator = sum(number_credible_sets >= k), denominator = n(),
      proportion = mean(number_credible_sets >= k), .groups = "drop"))
  list(traits = traits, lead_summary = lead_summary, curves = curves)
}

fig1_category_counts <- function(membership, categories) {
  long <- membership %>% select(assay, gene_id, gene_name, all_of(categories)) %>%
    distinct() %>% pivot_longer(all_of(categories), names_to = "category", values_to = "member")
  counts <- long %>% group_by(assay, category) %>%
    summarise(n_genes = n_distinct(gene_id[member %in% TRUE]), unknown = sum(is.na(member)), .groups = "drop")
  lists <- long %>% filter(member %in% TRUE) %>% distinct(assay, category, gene_id, gene_name)
  list(counts = counts, lists = lists,
       distinct = lists %>% group_by(category) %>% summarise(n_genes = n_distinct(gene_id), .groups = "drop"))
}

fig1_go_descendants <- function(nodes, edges, root) {
  edges <- edges %>% filter(child %in% nodes$go_id, parent %in% nodes$go_id,
    predicate %in% c("is_a", "http://purl.obolibrary.org/obo/BFO_0000050"))
  seen <- root
  repeat {
    expanded <- union(seen, edges$child[edges$parent %in% seen])
    if (length(expanded) == length(seen)) break
    seen <- expanded
  }
  seen
}
