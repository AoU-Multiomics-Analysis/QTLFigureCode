# Symmetric credible-set background and MAF models from Fig3/symmetric.R.
# Sourced after preparation/Fig3/methods.R. Rendering is kept in notebooks.
# The portable GTF reader in current_methods.R replaces the gzcat command.

SYM_PAIRS <- list(
  `eQTL-pQTL` = c(ann = "eQTL", oth = "pQTL"),
  `sQTL-pQTL` = c(ann = "sQTL", oth = "pQTL"),
  `eQTL-sQTL` = c(ann = "eQTL", oth = "sQTL"))

SYM_SUSIE_KEY <- c(eQTL = "susie_eqtl", sQTL = "susie_sqtl", pQTL = "susie_pqtl")

SPLICEAI_COLS <- c("splice_ai_acceptor_gain_score", "splice_ai_acceptor_loss_score",
                   "splice_ai_donor_gain_score",    "splice_ai_donor_loss_score")

CLEAN_CODING  <- c("missense_variant", "synonymous_variant",
                   "stop_gained", "frameshift_variant")

SPLICEAI_MAX  <- 0.2

SYM_SUSIE_COLS <- function()
  c("molecular_trait_id", "variant", "cs_id", "pip", "cs_size", "gene_id",
    "gene_type", "distTSS", "consequence", "phylop",
    FEATURE_CCRE_COLS, AF_POP_COLS, SPLICEAI_COLS)

.SYM_CS_SPLIT <- "_(?=ENSG[0-9]+:)"

.strip_gene_prefix <- function(x) sub("^ENSG[0-9]+:", "", x)

.partner_cs_ids <- function(right) {
  right <- right[nzchar(right)]
  if (!length(right)) return(character(0))
  unique(.strip_gene_prefix(unlist(
    strsplit(right, .SYM_CS_SPLIT, perl = TRUE), use.names = FALSE)))
}

.bare_gene <- function(x) stringr::str_extract(x, "ENSG[0-9]+")

.canonical_cds <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    p <- input_path("gene_gtf")
    message("  reading canonical CDS from ", basename(p), " ...")
    g <- data.table::fread(
      cmd = sprintf("gzcat %s | awk -F'\\t' '$3==\"CDS\"'", shQuote(p)),
      sep = "\t", header = FALSE, showProgress = FALSE,
      select = c(1, 4, 5, 9), col.names = c("chrom", "start", "end", "attr"))
    g[, gene  := sub('.*gene_id "([^".]+).*', "\\1", attr)]
    g[, canon := grepl('tag "Ensembl_canonical"', attr, fixed = FALSE)]
    g[, attr := NULL]
    # A backreference written "\\1" inside an Rscript -e shell string reaches R
    # as "\1" (a control character) and blanks every id, which then reads as
    # "no variant is in any CDS". Assert rather than trust.
    if (!all(grepl("^ENSG[0-9]+$", utils::head(g$gene, 200))))
      stop("GENCODE gene_id parse failed -- got e.g. '", g$gene[1], "'", call. = FALSE)
    g <- g[canon == TRUE, .(chrom, start, end, gene)]
    if (!nrow(g)) stop("no Ensembl_canonical CDS records found in the GTF", call. = FALSE)
    data.table::setkey(g, chrom, start, end)
    cache <<- g
    g
  }
})

.in_canonical_cds <- function(variant, egene) {
  cds <- .canonical_cds()
  q <- data.table::data.table(
    chrom = sub("_.*", "", variant),
    pos   = suppressWarnings(as.integer(data.table::tstrsplit(variant, "_")[[2]])))
  qu <- unique(q[!is.na(pos), .(chrom, start = pos, end = pos)])
  # foverlaps prefixes the QUERY columns with "i.": i.start is the variant
  # position, `start` is the CDS start. Grouping on `start` matches ~nothing.
  h  <- data.table::foverlaps(qu, cds, by.x = c("chrom", "start", "end"), nomatch = 0L)
  bygene <- h[, .(genes = list(unique(gene))), by = .(chrom, pos = i.start)]
  idx <- bygene[data.table::data.table(chrom = q$chrom, pos = q$pos),
                on = c("chrom", "pos"), which = TRUE]
  mapply(function(i, e) !is.na(i) && e %in% bygene$genes[[i]],
         idx, egene, USE.NAMES = FALSE)
}

.clean_coding_features <- function(d) {
  d <- data.table::as.data.table(d)
  sai <- intersect(SPLICEAI_COLS, names(d))
  for (s in sai) data.table::set(d, j = s, value = suppressWarnings(as.numeric(d[[s]])))
  d[, sai_max := if (length(sai))
      do.call(pmax, c(.SD, list(na.rm = TRUE))) else NA_real_, .SDcols = sai]
  d[, egene := sub("\\..*", "", .bare_gene(gene_id))]
  d[, in_egene_canon := .in_canonical_cds(variant, egene)]
  cc <- intersect(FEATURE_CCRE_COLS, names(d))
  d[, ccre := Reduce(`|`, .SD), .SDcols = cc]
  d[, is_clean := in_egene_canon & (is.na(sai_max) | sai_max <= SPLICEAI_MAX) & !ccre]
  for (f in intersect(CLEAN_CODING, names(d))) {
    raw <- as.logical(d[[f]])
    data.table::set(d, j = paste0(f, "_raw"), value = raw)
    data.table::set(d, j = f, value = raw & d$is_clean)
  }
  d[]
}

.read_susie_symmetric <- local({
  cache <- new.env(parent = emptyenv())
  function(layer) {
    if (!is.null(cache[[layer]])) return(cache[[layer]])
    d <- data.table::as.data.table(
      data.table::fread(input_path(SYM_SUSIE_KEY[[layer]]),
                        select = SYM_SUSIE_COLS(), showProgress = FALSE))
    d[, gene := .bare_gene(gene_id)]
    d[, variant := sub("^chrchr", "chr", variant)]
    cache[[layer]] <- d
    message("    ", layer, ": ", format(nrow(d), big.mark = ","), " SuSiE rows")
    d
  }
})

symmetric_universe <- function(refresh = FALSE)
  build_or_cache("symmetric_universe", refresh = refresh, {
    sigf <- input_files("enloc_cs")
    out <- lapply(names(SYM_PAIRS), function(pname) {
      p <- SYM_PAIRS[[pname]]; ann <- p[["ann"]]; oth <- p[["oth"]]
      row <- sigf[sigf$type == pname, ]
      if (nrow(row) != 1L)
        stop("expected exactly one enloc_cs file for pair '", pname, "', got ",
             nrow(row), call. = FALSE)

      sig <- data.table::as.data.table(
        LoadEnlocTable(row$path[1], ENLOC_CS_COLS))
      sig <- data.table::as.data.table(
        add_fdr_threshold_label(sig, "RCP", 0.05, "Coloc", quiet = TRUE))
      sig[, `:=`(left  = sub("\\(@\\).*$", "", Signal),
                 right = sub("^[^(]*\\(@\\)", "", Signal))]
      sig[, `:=`(gene = .bare_gene(left), left_cs = .strip_gene_prefix(left))]
      genes <- unique(stats::na.omit(sig$gene))

      ann_cs <- sig[, .(coloc = any(Coloc)), by = .(cs_id = left_cs)]
      oth_ids <- .partner_cs_ids(sig[Coloc == TRUE]$right)

      su_a <- .read_susie_symmetric(ann)
      su_o <- .read_susie_symmetric(oth)
      oth_cs <- unique(su_o[gene %in% genes, .(cs_id)])
      oth_cs[, coloc := cs_id %in% oth_ids]

      orphan <- setdiff(oth_ids, oth_cs$cs_id)
      if (length(orphan))
        message(sprintf(paste0("    %s: %d of %d colocalizing %s credible sets ",
                               "are absent from the SuSiE table and cannot be ",
                               "recovered"),
                        pname, length(orphan), length(oth_ids), oth))

      va <- su_a[ann_cs, on = "cs_id", nomatch = 0][, layer := ann]
      vo <- su_o[gene %in% genes][oth_cs, on = "cs_id", nomatch = 0][, layer := oth]
      u  <- data.table::rbindlist(list(va, vo), use.names = TRUE, fill = TRUE)
      u[, pair := pname]
      message(sprintf("    %s: %s rows (%s %s / %s %s)", pname,
              format(nrow(u), big.mark = ","),
              format(nrow(va), big.mark = ","), ann,
              format(nrow(vo), big.mark = ","), oth))
      u
    })
    u <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)

    # max_af mirrors coloc_model_input(): fold each per-ancestry frequency to a
    # MINOR allele frequency, then take the maximum across ancestries.
    folded <- u[, AF_POP_COLS, with = FALSE]
    folded[] <- lapply(folded, function(x) ifelse(!is.na(x) & x > .5, 1 - x, x))
    u[, max_af := do.call(pmax, c(as.list(folded), na.rm = FALSE))]
    u <- u[!is.na(max_af) & max_af > 0 & !is.na(pip) & !is.na(cs_size)]
    u[, log_maf := log10(max_af)]

    u <- .clean_coding_features(build_feature_matrix(as.data.frame(u)))
    u[]
  })

.EMPTY_MAF_FIT <- data.table::data.table(
  feature = character(0), fold = numeric(0), lo = numeric(0), hi = numeric(0),
  p = numeric(0), n = integer(0))

fit_maf_enrichment <- function(d, features, min_variants = 5) {
  d <- data.table::as.data.table(d)
  # `layers` is dropped when it has a single level: lm() errors outright on a
  # one-level factor ("contrasts can be applied only to factors with 2 or more
  # levels"), which tryCatch would swallow into an empty result. In production
  # it always has 2-3 levels; on a small fixture it may not.
  adj <- if (length(unique(d$layers)) > 1L) "x + layers + pip" else "x + pip"
  fo  <- stats::as.formula(paste("log_maf ~", adj))
  data.table::rbindlist(lapply(features, function(f) {
    x <- as.logical(d[[f]]); x[is.na(x)] <- FALSE
    if (min(table(x)) < min_variants) return(NULL)
    m <- tryCatch(suppressWarnings(stats::lm(
           fo, data = cbind(d, x = x), weights = d$pip)),
         error = function(e) NULL)
    if (is.null(m)) return(NULL)
    co <- stats::summary.lm(m)$coefficients
    if (!"xTRUE" %in% rownames(co)) return(NULL)
    # NEGATE so >1 means RARER; 10^ NOT exp() -- the outcome is log10 MAF, and
    # exponentiating base e understates every effect (see NOTES.md).
    b <- -co["xTRUE", 1]; se <- co["xTRUE", 2]
    data.table::data.table(
      feature = f, fold = 10^b, lo = 10^(b - 1.959964 * se),
      hi = 10^(b + 1.959964 * se), p = co["xTRUE", 4], n = sum(x))
  }), fill = TRUE) -> r
  if (!nrow(r)) return(data.table::copy(.EMPTY_MAF_FIT))
  r
}

collapse_to_variants <- function(u, features) {
  u <- data.table::as.data.table(u)
  gmean <- function(x) exp(mean(log(pmax(x, 1e-12))))
  u[, c(.(pip = gmean(pip), log_maf = log_maf[1],
          layers = paste(sort(unique(layer)), collapse = "+")),
        lapply(.SD, function(z) any(as.logical(z)))),
    by = variant, .SDcols = features]
}

coloc_maf_enrichment <- function(refresh = FALSE)
  build_or_cache("coloc_maf_enrichment", refresh = refresh, {
    u <- symmetric_universe()
    u <- u[pair == "eQTL-pQTL" & coloc == TRUE]
    feats <- setdiff(feature_names(as.data.frame(u)),
                     c("rare", "is_clean", "ccre", "in_egene_canon",
                       paste0(CLEAN_CODING, "_raw")))
    v <- collapse_to_variants(u, feats)
    v[, layers := factor(layers)]
    r <- fit_maf_enrichment(v, feats)
    r[, padj := if (nrow(r)) stats::p.adjust(p, "BH") else numeric(0)]
    message(sprintf("  %s colocalizing variants | %d annotations fitted",
                    format(nrow(v), big.mark = ","), nrow(r)))
    r[]
  })

PANEL_E_LABELS <- c(
  stop_gained = "Stop gained", frameshift_variant = "Frameshift",
  splice_acceptor_variant = "Splice acceptor", splice_donor_variant = "Splice donor",
  splice_region_variant = "Splice region", conserved = "Conserved",
  PLS = "Promoter (PLS)", `3_prime_UTR_variant` = "3' UTR",
  non_coding_transcript_exon_variant = "Non-coding exon",
  non_coding_transcript_variant = "Non-coding transcript",
  NMD_transcript_variant = "NMD transcript",
  `5_prime_UTR_variant` = "5' UTR", tss_proximal = "TSS-proximal",
  CTCF.bound = "CTCF-bound", pELS = "Proximal enhancer",
  upstream_gene_variant = "Upstream", intron_variant = "Intron",
  dELS = "Distal enhancer", missense_variant = "Missense",
  synonymous_variant = "Synonymous")

PANEL_E_MIN_N <- 5
