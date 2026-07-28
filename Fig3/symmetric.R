# =============================================================================
# Figure 3 — symmetric colocalization background and the panel E rare-variant
# enrichment. Sourced by fig3/00_setup.R; not run directly.
#
# WHY THIS EXISTS
# ---------------
# FastENLOC enumerates credible sets from ONE layer only. In `Signal`, the half
# before "(@)" is the QTL/annotation side; the half after it is the partner
# layer, and it is populated ONLY when the two colocalize. So every
# non-colocalizing row of .snp.out/.sig.out is an annotation-side-specific
# credible set, and the partner layer's own unshared signals are absent from the
# output entirely.
#
# For eQTL-pQTL that meant the background was 11,065 eQTL-specific credible sets
# and ZERO pQTL-specific ones -- colocalizing variants were being compared only
# against eQTL-only signals. This file reconstructs the missing side: the
# partner layer's credible sets at the genes the pair was tested on, minus the
# ones enloc reports as colocalizing.
#
# See Fig3/NOTES.md for what this changed.
# =============================================================================

# Pair -> (annotation layer = enloc's "qtl" side, partner layer = "gwas" side).
# The protein layer is always the partner; it is never the annotation side.
SYM_PAIRS <- list(
  `eQTL-pQTL` = c(ann = "eQTL", oth = "pQTL"),
  `sQTL-pQTL` = c(ann = "sQTL", oth = "pQTL"),
  `eQTL-sQTL` = c(ann = "eQTL", oth = "sQTL"))

SYM_SUSIE_KEY <- c(eQTL = "susie_eqtl", sQTL = "susie_sqtl", pQTL = "susie_pqtl")

SPLICEAI_COLS <- c("splice_ai_acceptor_gain_score", "splice_ai_acceptor_loss_score",
                   "splice_ai_donor_gain_score",    "splice_ai_donor_loss_score")

# Consequence classes redefined as "clean canonical" calls (see
# .clean_coding_features). splice_region_variant is deliberately NOT in this
# list: the SpliceAI criterion would discard exactly the variants that make a
# splice-region call meaningful, and it belongs with splice_donor/acceptor,
# which are not cleaned either.
CLEAN_CODING  <- c("missense_variant", "synonymous_variant",
                   "stop_gained", "frameshift_variant")
SPLICEAI_MAX  <- 0.2      # above this a coding variant may act through splicing

# Columns pulled from the SuSiE tables for the symmetric universe. Deliberately
# explicit: a silently-missing column here would shrink the universe rather than
# fail.
SYM_SUSIE_COLS <- function()
  c("molecular_trait_id", "variant", "cs_id", "pip", "cs_size", "gene_id",
    "gene_type", "distTSS", "consequence", "phylop",
    FEATURE_CCRE_COLS, AF_POP_COLS, SPLICEAI_COLS)

# -----------------------------------------------------------------------------
# enloc Signal parsing
# -----------------------------------------------------------------------------
# Two traps, both of which fail SILENTLY by returning an empty join rather than
# an error, and both of which cost real debugging time:
#
#  1. BOTH halves of `Signal` carry an "ENSG<digits>:" gene prefix that the
#     SuSiE cs_id does not (expression "ENSG...L2", splicing
#     "chr:...:clu:ENSG...L2", proteomics "UNIPROT_ENSG...L2"). Leaving it on
#     matches nothing.
#  2. A signal may name SEVERAL partner credible sets, separated by "_" followed
#     by a NEW gene-prefixed id. Splitting on a plain "_ENSG" instead shreds the
#     UniProt-prefixed pQTL ids (Q6PCB0_ENSG00000179403.12_L2) into fragments.
.SYM_CS_SPLIT <- "_(?=ENSG[0-9]+:)"

.strip_gene_prefix <- function(x) sub("^ENSG[0-9]+:", "", x)

# Partner-side credible-set ids named by a signal, as bare SuSiE cs_ids.
.partner_cs_ids <- function(right) {
  right <- right[nzchar(right)]
  if (!length(right)) return(character(0))
  unique(.strip_gene_prefix(unlist(
    strsplit(right, .SYM_CS_SPLIT, perl = TRUE), use.names = FALSE)))
}

# gene_id is strand-prefixed in the splicing table ("-:ENSG..."), so pull the
# accession out rather than assuming the field starts with it.
.bare_gene <- function(x) stringr::str_extract(x, "ENSG[0-9]+")

# -----------------------------------------------------------------------------
# GENCODE canonical CDS
# -----------------------------------------------------------------------------
# Positions falling in the CDS of an Ensembl-canonical transcript, keyed by
# gene. Used to demand that a "missense" call refers to the eGENE's own dominant
# protein rather than a neighbouring gene or a minor isoform.
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

# TRUE where the variant lies in a canonical CDS of `egene`.
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

# -----------------------------------------------------------------------------
# Clean canonical coding calls
# -----------------------------------------------------------------------------
# `consequence` is a VARIANT-level VEP call on whichever transcript VEP picked.
# A coding class is kept only where the amino-acid change is plausibly the
# mechanism AND refers to the eGene's dominant protein:
#   (a) in a canonical CDS of the eGene itself,
#   (b) SpliceAI max delta <= SPLICEAI_MAX (not acting through splicing),
#   (c) no ENCODE cCRE overlap (not plausibly regulatory).
# The pre-cleaning column is retained as <feature>_raw so nothing is lost.
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

# -----------------------------------------------------------------------------
# The symmetric universe
# -----------------------------------------------------------------------------
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

# One row per (variant, credible set, layer) for BOTH layers of each pair, each
# row labelled with whether its credible set colocalizes.
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

# -----------------------------------------------------------------------------
# Panel E: MAF enrichment among colocalizing variants
# -----------------------------------------------------------------------------
# Pure seam: takes the collapsed per-variant frame, returns the tidy fit. Kept
# free of I/O so it is testable without the staged inputs.
# Zero fitted annotations is a legitimate outcome (a small fixture, or a pair
# where nothing clears the floor) -- not an error. Returning an UNtyped empty
# table makes `padj := p.adjust(p)` coerce NAs and every downstream filter fail
# on a missing column, so the empty shape is declared once, here.
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

# Collapse the universe to one row per variant. A colocalizing variant sits in
# BOTH layers' credible sets for ~31% of keys; without this it enters twice as
# though independent. The weight is the GEOMETRIC mean of its PIPs, which is
# high only when the variant is confidently fine-mapped in BOTH layers -- what
# colocalization evidence actually requires.
collapse_to_variants <- function(u, features) {
  u <- data.table::as.data.table(u)
  gmean <- function(x) exp(mean(log(pmax(x, 1e-12))))
  u[, c(.(pip = gmean(pip), log_maf = log_maf[1],
          layers = paste(sort(unique(layer)), collapse = "+")),
        lapply(.SD, function(z) any(as.logical(z)))),
    by = variant, .SDcols = features]
}

# Panel E's estimand: among COLOCALIZING eQTL-pQTL variants, does each
# annotation mark rarer or commoner variants?
#
# Restricted to eQTL-pQTL because it is the only pair whose partner side is
# fully recoverable (0 orphan credible sets; eQTL-sQTL loses 24.5%).
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

# -----------------------------------------------------------------------------
# Panel E
# -----------------------------------------------------------------------------
# Styling deliberately mirrors panel C: grey40 confidence lines at the default
# weight and a solid size-4 point, so the two read as siblings. Colour is the
# one departure -- it marks direction, which C does not need.
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

PANEL_E_MIN_N <- 5    # matches fit_feature_enrichment()'s floor; see NOTES.md

build_panelE_maf <- function(enr, min_n = PANEL_E_MIN_N) {
  d <- dplyr::filter(as.data.frame(enr), padj < 0.05, n >= min_n)

  # An annotation missing from PANEL_E_LABELS would vanish from the panel
  # without a word -- exactly how stop_gained/frameshift disappeared during
  # development. Name them instead.
  unlabelled <- setdiff(d$feature, names(PANEL_E_LABELS))
  if (length(unlabelled))
    message("  note: significant but unlabelled, so not plotted: ",
            paste(unlabelled, collapse = ", "))
  d <- dplyr::filter(d, feature %in% names(PANEL_E_LABELS))
  if (!nrow(d)) return(blank_plot)

  d <- d %>%
    dplyr::mutate(dir = ifelse(fold > 1, "rare", "common")) %>%
    dplyr::arrange(fold) %>%
    dplyr::mutate(label = factor(PANEL_E_LABELS[feature],
                                 levels = PANEL_E_LABELS[feature]))

  dir_col <- c(rare = unname(SIGNAL_FILL["5"]), common = unname(LAYER_FILL["3"]))
  # Limits are derived from the data, not hard-coded: a CI end outside the scale
  # is dropped by geom_linerange with only a warning.
  x_min <- min(d$lo) * 0.85
  x_max <- max(d$hi) * 1.25
  breaks <- c(0.5, 0.7, 1, 2, 5, 10, 20)
  breaks <- breaks[breaks > x_min & breaks < x_max]
  lx <- function(f) 10^(f * (log10(x_max) - log10(x_min)))
  y_ann <- nrow(d) + 0.85

  ggplot(d, aes(x = fold, y = label)) +
    annotate("rect", xmin = x_min, xmax = 1, ymin = -Inf, ymax = Inf,
             fill = "grey95", color = NA) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_linerange(aes(xmin = lo, xmax = hi), color = "grey40") +
    geom_point(aes(color = dir), size = 4) +
    annotate("segment", x = lx(-0.004), xend = lx(-0.055), y = y_ann, yend = y_ann,
             arrow = arrow(length = unit(5, "pt"), type = "closed"),
             color = dir_col[["common"]], linewidth = 0.4) +
    annotate("segment", x = lx(0.004), xend = lx(0.055), y = y_ann, yend = y_ann,
             arrow = arrow(length = unit(5, "pt"), type = "closed"),
             color = dir_col[["rare"]], linewidth = 0.4) +
    annotate("text", x = lx(-0.065), y = y_ann, label = "common", hjust = 1,
             size = SZ_ANNOT, color = dir_col[["common"]]) +
    annotate("text", x = lx(0.065), y = y_ann, label = "rare", hjust = 0,
             size = SZ_ANNOT, color = dir_col[["rare"]]) +
    scale_color_manual(values = dir_col, guide = "none") +
    scale_x_log10(breaks = breaks, labels = breaks,
                  limits = c(x_min, x_max), expand = expansion(mult = 0)) +
    scale_y_discrete(expand = expansion(add = c(0.6, 1.5))) +
    labs(x = "MAF fold change vs. other colocalizing variants", y = NULL) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(l = 12, r = 24, t = 12, b = 6, unit = "pt"))
}
