# =============================================================================
# Figure 3 — shared setup (sourced by every panel / supplemental / assembly).
# Run all scripts from the PROJECT ROOT. Each script starts with:
#   source("fig3/00_setup.R")
# Re-sourcing is cheap and idempotent.
#
# Figure 3 covers cross-QTL colocalization (FastENLOC) between each pair of
# molecular assays — eQTL-pQTL, sQTL-pQTL, eQTL-sQTL — at gene, credible-set
# (signal) and variant level.
#
# See docs/superpowers/specs/2026-07-23-fig3-cleanup-design.md
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(igraph)
  library(broom)
  library(scales)
  library(showtext)
  library(sysfonts)
  library(paletteer)
})

# --- Paths (override via env vars if needed) --------------------------------
# INPUT is the staging root; it defaults to the current working directory, so
# gsutil-staged files land under the cwd in their manifest subpaths (enloc/,
# finemapping/, ...). Set FIG3_INPUT to keep inputs in a separate folder.
INPUT    <- Sys.getenv("FIG3_INPUT",    ".")            # staging root (cwd)
DATA_DIR <- Sys.getenv("FIG3_DATA",     "fig3/data")    # memoized derived tables
OBJ_DIR  <- Sys.getenv("FIG3_OBJ",      "fig3/objects") # per-panel ggplot objects
OUTDIR   <- Sys.getenv("FIG3_OUTDIR", "fig3")         # final Figure3.pdf/png
MANIFEST <- Sys.getenv("FIG3_MANIFEST", "fig3/inputs.tsv")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OBJ_DIR,  showWarnings = FALSE, recursive = TRUE)

dir.create(OUTDIR,   showWarnings = FALSE, recursive = TRUE)

source("common/manifest.R")
source("common/vep.R")

# --- Font: Inter via Google if reachable, else default sans -----------------
FONT <- "sans"
if (!"inter" %in% sysfonts::font_families()) {
  if (isTRUE(tryCatch({ font_add_google("Inter", "inter"); TRUE }, error = function(e) FALSE)))
    FONT <- "inter"
} else FONT <- "inter"

DPI <- 400                     # showtext DPI; must match every ggsave() DPI
showtext_opts(dpi = DPI)
showtext_auto()

# --- Design system ----------------------------------------------------------
# element_text() sizes are in pt; geom_text()/geom_label()/annotate() are in mm.
# Figure 3 is an 18x12" canvas, so these are much smaller than Figure 1's.
SZ_AXIS_TITLE <- 11; SZ_AXIS_TEXT <- 9;  SZ_STRIP  <- 12
SZ_LEG_TEXT   <- 9;  SZ_LEG_TITLE <- 10; SZ_TAG    <- 14
SZ_TITLE      <- 11
SZ_ANNOT      <- 2.7           # mm — geom_label()/annotate() only

theme_set(theme_classic())
theme_update(
  text             = element_text(family = FONT),
  axis.title       = element_text(size = SZ_AXIS_TITLE),
  axis.text        = element_text(size = SZ_AXIS_TEXT),
  strip.text       = element_text(size = SZ_STRIP),
  strip.background = element_blank(),
  legend.text      = element_text(size = SZ_LEG_TEXT),
  legend.title     = element_text(size = SZ_LEG_TITLE)
)

blank_plot <- ggplot() + theme_void(base_family = FONT)

# Panel A fill — number of colocalized molecular layers.
LAYER_FILL <- c("0" = "#D9D9D9", "1" = "#A9C1D9", "2" = "#7B9AB8", "3" = "#577590")

# Panel B fill — number of colocalized signals.
SIGNAL_FILL <- c("0" = "#E9DCDC", "1" = "#D8BBBB", "2" = "#C89A9A", "3" = "#B37979",
                 "4" = "#9E5E5E", "5" = "#8A4D4D", "6" = "#763F3F", "7" = "#643232")

# Per-ancestry allele-frequency columns, used to derive a per-variant maximum
# MAF. They ride along on the annotated SuSiE tables, which is why Fig3 no
# longer reads the standalone AoU annotation table (see NOTES.md). Named
# explicitly so a schema change fails loudly rather than silently changing
# which columns are folded.
AF_POP_COLS <- c("gvs_afr_af", "gvs_amr_af", "gvs_eas_af", "gvs_eur_af",
                 "gvs_mid_af", "gvs_sas_af", "gvs_oth_af")

# cCRE flags carried through for the panel E enrichment. Already logical in the
# SuSiE tables.
FEATURE_CCRE_COLS <- c("PLS", "pELS", "dELS", "CTCF.bound", "CTCF.only",
                       "DNase.H3K4me3", "FANTOM5")
# Extra columns the enrichment needs beyond the base fine-mapping fields.
FEATURE_EXTRA_COLS <- c("consequence", "phylop", "distTSS", "gene_type", "cs_size")

# =============================================================================
# FastENLOC output readers
# -----------------------------------------------------------------------------
# The .sig.out and .snp.out files are whitespace-aligned and their header line
# does not line up with the data columns, so they are parsed positionally: drop
# the header, split each row on runs of whitespace, and apply known names.
# (.gene.out is a normal TSV and is read with fread.)
# =============================================================================
ENLOC_CS_COLS      <- c("Signal", "Num_SNP", "CPIP_qtl", "CPIP_gwas_marginal",
                        "CPIP_gwas_qtl_prior", "RCP", "LCP")
ENLOC_VARIANT_COLS <- c("Signal", "SNP", "PIP_qtl", "PIP_gwas_marginal",
                        "PIP_gwas_qtl_prior", "SCP")

LoadEnlocTable <- function(path, cols) {
  x <- readLines(path)
  # Line numbers in error messages below are 1-based into this original file,
  # so track them before dropping the header line.
  line_no <- seq_along(x)[-1]
  x <- x[-1]                                   # drop the header line
  keep <- nzchar(trimws(x))                    # drop blank lines
  x <- x[keep]
  line_no <- line_no[keep]

  # Column type policy, declared once and applied identically on both the
  # empty and non-empty paths below: "Signal"/"SNP" are identifiers and stay
  # character; every other FastENLOC column (SNP counts, posterior/prior
  # probabilities) is a numeric measurement. This makes the two paths
  # structurally identical by construction rather than by coincidence of what
  # type.convert() infers from a particular file's values -- previously the
  # empty-file path returned all-character columns, which broke bind_rows()
  # against a populated file's integer/numeric columns.
  is_id <- cols %in% c("Signal", "SNP")

  if (length(x) == 0L) {
    # Header-only file: a gene pair with no signals is a legitimate FastENLOC
    # output, not an error. Return a 0-row tibble with the right columns and
    # types so downstream bind_rows()/split() code keeps working.
    out <- lapply(is_id, function(id) if (id) character(0) else numeric(0))
    names(out) <- cols
    return(tibble::as_tibble(out))
  }

  rows <- strsplit(trimws(x), "\\s+")
  # do.call(rbind, ...) SILENTLY RECYCLES when row lengths differ, so a ragged
  # row (wrong token count) would land its values in the wrong columns without
  # ever tripping the length check below. Check lengths before rbind instead.
  n <- lengths(rows)
  bad <- which(n != length(cols))
  if (length(bad))
    stop("Expected ", length(cols), " columns in ", path, ", found ",
         paste(unique(n[bad]), collapse = "/"), " on line(s) ",
         paste(head(line_no[bad], 5), collapse = ", "),
         if (length(bad) > 5) ", ..." else "")

  mat <- do.call(rbind, rows)
  out <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(out) <- cols
  # strsplit() yields character; apply the same identifier/measurement policy
  # declared above instead of inferring types via type.convert(). A
  # measurement column that holds something non-numeric surfaces as NA (with
  # R's standard "NAs introduced by coercion" warning) rather than silently
  # staying character and breaking a later numeric comparison.
  out[!is_id] <- lapply(out[!is_id], as.numeric)
  tibble::as_tibble(out)
}

# Label rows as colocalized by controlling the Bayesian FDR at `control_level`:
# sort local FDRs (1 - posterior) ascending, take the largest prefix whose
# running mean stays at or below the target, and threshold at its posterior.
#
# `prob_col` and `label_col` are plain COLUMN-NAME STRINGS. The original used
# rlang::ensym() here, but every call is internal and passes a name held in a
# variable, which is exactly the case non-standard evaluation cannot handle.
#
# The original stashed the threshold in attr(out, "threshold"), which the
# following dplyr verbs silently dropped; it is reported via message() instead.
add_fdr_threshold_label <- function(df, prob_col, control_level = 0.05,
                                    label_col = "Coloc", quiet = FALSE) {
  stopifnot(is.character(prob_col), length(prob_col) == 1L)
  if (!prob_col %in% names(df))
    stop("Column '", prob_col, "' not found; have: ", paste(names(df), collapse = ", "))

  probs <- as.numeric(df[[prob_col]])

  # Rank by probability (best first) and take the local FDR from that order,
  # rather than sorting lfdr directly. This lets the eventual cutoff be read
  # back as an actual element of `probs` (see below) instead of reconstructed
  # via `1 - lfdr`.
  ord  <- order(probs, decreasing = TRUE, na.last = NA)  # NA probs excluded
  lfdr <- 1 - probs[ord]
  FDR  <- cumsum(lfdr) / seq_along(lfdr)
  # Tolerance guards an exact-boundary case (e.g. probs = 0.95, control_level =
  # 0.05) where `1 - 0.95` rounds to a double a hair above 0.05, which would
  # otherwise make a should-pass row fail `FDR <= control_level`. Double-
  # precision error is ~1e-16, so 1e-8 is generous without risking false passes.
  valid <- which(FDR <= control_level + 1e-8)
  k <- if (length(valid)) max(valid) else 0L

  # cutoff must be an observed value from `probs` itself (never `1 - lfdr`):
  # that round-trip is inexact for p < 0.5 in IEEE-754, which could silently
  # exclude the very row that defined the boundary from `probs >= cutoff`.
  # Using an actual element also guarantees every row tied with the boundary
  # is included, so significance never depends on arbitrary sort order.
  cutoff <- if (k > 0L) probs[ord[k]] else Inf

  # A row with no posterior is not a colocalization; label it FALSE rather
  # than letting NA propagate into the label column and poison downstream
  # sum()/mean() calls that count colocalizations.
  label <- !is.na(probs) & probs >= cutoff

  if (!quiet)
    message(sprintf("    %s >= %.4g at FDR %.0f%%  (%d / %d rows labelled)",
                    prob_col, cutoff, 100 * control_level,
                    sum(label), length(probs)))

  df[[label_col]] <- label
  df
}

# =============================================================================
# Derived-data cache
# -----------------------------------------------------------------------------
# Fig3's panels share a heavy upstream: the three SuSiE tables, which are also
# the source of the per-variant allele frequencies (see .variant_feature_map()
# below) now that the standalone AoU annotation table is gone from the manifest. Under
# strict per-panel loading those tables would be read, collapsed and joined once
# per panel. Each builder below memoizes its result to DATA_DIR/<name>.rds, so a
# panel still calls only the builders it needs while the expensive work happens
# once.
#
# `expr` is a lazily-evaluated promise: on a cache hit it is never evaluated, so
# a cached panel never reads its inputs at all.
#
# THE CACHE IS KEYED BY NAME ONLY, not by input content. Clear DATA_DIR (or pass
# refresh = TRUE) after changing an input file or a builder.
# =============================================================================
build_or_cache <- function(name, expr, refresh = FALSE) {
  f <- file.path(DATA_DIR, paste0(name, ".rds"))
  if (!refresh && file.exists(f)) return(readRDS(f))
  message("building ", name, " ...")
  val <- force(expr)
  saveRDS(val, f)
  val
}

# =============================================================================
# SuSiE fine-mapping (COMB)
# =============================================================================
.add_cs_stats <- function(d) {
  d %>%
    group_by(molecular_trait_id, cs_index) %>%
    mutate(cpip = sum(pip), number_variants = dplyr::n()) %>%
    ungroup()
}

# Columns every SuSiE table must supply. AF_POP_COLS are carried through so
# .variant_feature_map() can source allele frequencies here instead of from a
# separate annotation table.
SUSIE_BASE_COLS <- c("molecular_trait_id", "variant", "pip", "cs_index",
                     "posterior_mean", "cs_id", "posterior_sd")

# Read one SuSiE table for `key`, validating its schema first.
#
# Fig3 reads the SAME tables as Fig2. These are combined-ancestry (COMB)
# fine-mapping exports, one per assay, so there is nothing to select on: the old
# pQTL read's filter(group == "COMB") is gone. It would in fact have matched
# NOTHING here -- the `group` column holds "COMB.proteomics", not "COMB" -- and
# handed the pipeline a silently empty table.
#
# The schema is checked up front and reports EVERY missing column: a bare
# select() error names only the first, and these files come from an upstream
# pipeline this repo does not control, so a partial message sends you round the
# loop once per column.
.read_susie_table <- function(key, extra_cols = character()) {
  path <- input_path(key)
  d    <- fread(path)
  need <- c(SUSIE_BASE_COLS, extra_cols, AF_POP_COLS)
  absent <- setdiff(need, names(d))
  if (length(absent))
    stop("SuSiE table for key '", key, "' (", basename(path), ") is missing ",
         "required column(s): ", paste(absent, collapse = ", "),
         ".\n  Fig3 reads the same SuSiE tables as Fig2; a differently ",
         "post-processed export may not carry the credible-set columns ",
         "(cs_id / cs_index / posterior_sd) that Fig3's colocalization needs.",
         call. = FALSE)
  d
}

susie_finemapping <- function(refresh = FALSE)
  build_or_cache("susie_finemapping", refresh = refresh, {
    eQTL <- .read_susie_table("susie_eqtl", "gene_name") %>%
      select(all_of(SUSIE_BASE_COLS), gene_name, all_of(AF_POP_COLS),
             any_of(c(FEATURE_CCRE_COLS, FEATURE_EXTRA_COLS))) %>%
      mutate(molecular_trait_id = str_remove(molecular_trait_id, "\\..*"),
             variant = str_replace(variant, "chrchr", "chr")) %>%
      .add_cs_stats()

    sQTL <- .read_susie_table("susie_sqtl") %>%
      mutate(molecular_trait_id = str_extract(molecular_trait_id, "ENSG[0-9]+")) %>%
      select(all_of(SUSIE_BASE_COLS), all_of(AF_POP_COLS),
             any_of(c(FEATURE_CCRE_COLS, FEATURE_EXTRA_COLS))) %>%
      .add_cs_stats()

    pQTL <- .read_susie_table("susie_pqtl") %>%
      mutate(molecular_trait_id = str_extract(molecular_trait_id, "ENSG[0-9]+")) %>%
      select(all_of(SUSIE_BASE_COLS), all_of(AF_POP_COLS),
             any_of(c(FEATURE_CCRE_COLS, FEATURE_EXTRA_COLS))) %>%
      .add_cs_stats()

    cs_list <- bind_rows(
      eQTL %>% distinct(molecular_trait_id, cs_id) %>% mutate(type = "eQTL"),
      pQTL %>% distinct(molecular_trait_id, cs_id) %>% mutate(type = "pQTL"),
      sQTL %>% distinct(molecular_trait_id, cs_id) %>% mutate(type = "sQTL")
    )

    # NOTE: this is the INTERSECTION of the three assays. The original called it
    # QTL_union_gene_list, which was a misnomer.
    genes_all_three_assays <- eQTL %>%
      filter(molecular_trait_id %in% pQTL$molecular_trait_id,
             molecular_trait_id %in% sQTL$molecular_trait_id) %>%
      distinct(molecular_trait_id) %>%
      pull(molecular_trait_id)

    protein_map <- pQTL %>%
      mutate(UniProt = str_remove(cs_id, "_.*")) %>%
      distinct(UniProt, molecular_trait_id)

    gene_map <- eQTL %>%
      distinct(molecular_trait_id, gene_name) %>%
      left_join(protein_map, by = "molecular_trait_id")

    list(eQTL = eQTL, sQTL = sQTL, pQTL = pQTL, cs_list = cs_list,
         gene_map = gene_map, genes_all_three_assays = genes_all_three_assays)
  })

# =============================================================================
# FastENLOC colocalization
# -----------------------------------------------------------------------------
# The FDR threshold is computed PER ASSAY PAIR, matching the original (which
# labelled each file before binding). `type` arrives attached from the manifest,
# so the original's per-frame mutate(type = '...') calls are gone.
# =============================================================================
.label_by_pair <- function(d, prob_col, level_label) {
  d %>%
    split(.$type) %>%
    purrr::imap_dfr(function(x, ty) {
      message("  ", ty, " (", level_label, "):")
      add_fdr_threshold_label(x, prob_col)
    })
}

coloc_gene_level <- function(refresh = FALSE)
  build_or_cache("coloc_gene_level", refresh = refresh, {
    load_input("enloc_gene") %>%
      .label_by_pair("GLCP", "gene level")
  })

coloc_cs_level <- function(refresh = FALSE)
  build_or_cache("coloc_cs_level", refresh = refresh, {
    d <- load_input("enloc_cs", reader = LoadEnlocTable, cols = ENLOC_CS_COLS) %>%
      .label_by_pair("RCP", "signal level")

    # The documented FastENLOC Signal format is exactly two sides,
    # "gene:cs1(@)gene:cs2". A second side listing MULTIPLE credible sets
    # separates them with "_", not another "(@)" (a later task handles that
    # case), so more than one "(@)" here really does mean a malformed row
    # rather than a legitimate shape. Validate explicitly and stop loudly --
    # a warning is not enough because build_or_cache would cache the bad
    # result -- rather than let extra = "merge" silently fold any surplus
    # parts into cs_2.
    n_at <- str_count(d$Signal, fixed("(@)"))
    bad  <- which(n_at != 1L)
    if (length(bad))
      stop(length(bad), " row(s) of Signal do not contain exactly one '(@)' ",
           "separator (found ", paste(unique(n_at[bad]), collapse = "/"),
           "); first offending value: '", d$Signal[bad[1]], "'")

    # Exactly one "(@)" per row is now guaranteed, so the split below can
    # safely assume two parts.
    d %>%
      separate(Signal, into = c("cs_1", "cs_2"), sep = "\\(@\\)",
               remove = FALSE, fill = "right", extra = "merge") %>%
      mutate(cs_1 = str_remove(cs_1, "^ENSG[0-9]+:"),
             cs_2 = str_remove(cs_2, "^ENSG[0-9]+:"))
  })

# One row per variant, carrying allele frequencies AND the annotation features
# panel E's enrichment needs.
#
# cs_size is a property of a CREDIBLE SET and one variant can sit in several
# across genes and assays, so the map carries the median over that variant's
# rows. Keeping exactly one row per variant is what the many-to-one join guard
# in coloc_variant_level() relies on.
.variant_feature_map <- function(fm) {
  keep <- c("variant", AF_POP_COLS, FEATURE_CCRE_COLS, FEATURE_EXTRA_COLS)
  raw <- bind_rows(lapply(list(fm$eQTL, fm$sQTL, fm$pQTL),
                          function(d) dplyr::select(d, any_of(keep)))) %>%
    # The eQTL reader already repairs "chrchr"; normalise all three here so a
    # stray prefix cannot split one variant into two rows.
    mutate(variant = str_replace(variant, "chrchr", "chr"),
           across(all_of(AF_POP_COLS), as.numeric))

  size <- raw %>%
    group_by(variant) %>%
    summarise(cs_size_med = stats::median(as.numeric(cs_size), na.rm = TRUE),
              .groups = "drop")

  af <- raw %>% select(-any_of("cs_size"))
  af <- af[rowSums(!is.na(af[AF_POP_COLS])) > 0, , drop = FALSE] %>% distinct()

  # A duplicate here means two SuSiE tables disagree about the same variant's
  # variant-level annotations, which is a genuine upstream inconsistency rather
  # than something to silently pick a winner for.
  dup <- unique(af$variant[duplicated(af$variant)])
  if (length(dup))
    stop(length(dup), " variant(s) carry CONFLICTING variant-level values ",
         "across the SuSiE fine-mapping tables -- e.g. '", dup[1],
         "'. These are properties of the variant alone; resolve the ",
         "disagreement upstream.", call. = FALSE)
  left_join(af, size, by = "variant")
}

# Cut-points for the binarized covariates. Conventional defaults, NOT tuned --
# if a result turns on the exact value that is a finding about fragility and
# belongs in NOTES.md.
PHYLOP_CONSERVED <- 2
TSS_PROXIMAL_BP  <- 2000

# Columns that are logical but are NOT features (outcomes or bookkeeping).
NON_FEATURE_LOGICAL <- c("Coloc", "coloc", "flagged")

# Turn the raw annotation columns into the binary feature set: four binarized
# covariates plus one logical column per VEP consequence term.
build_feature_matrix <- function(d) {
  d %>%
    mutate(
      rare           = !is.na(max_af)    & max_af   <  0.01,
      conserved      = !is.na(phylop)    & phylop   >  PHYLOP_CONSERVED,
      # Sign-agnostic: distTSS is signed by strand, proximity is not.
      tss_proximal   = !is.na(distTSS)   & abs(distTSS) < TSS_PROXIMAL_BP,
      protein_coding = !is.na(gene_type) & gene_type == "protein_coding") %>%
    expand_vep_consequences()
}

# Logical columns eligible as enrichment features.
feature_names <- function(d)
  setdiff(names(d)[vapply(d, is.logical, logical(1))], NON_FEATURE_LOGICAL)

# Adjusted enrichment of each feature on a binary outcome.
#
# The PIP + log10(cs_size) adjustment is the point of the design: colocalization
# requires good fine-mapping in BOTH layers, so an unadjusted test conflates
# biology with detectability.
#
# A feature is fitted only with at least `min_carriers` observations in the
# SMALLER of its two groups. Below that, glm returns separated or wildly
# unstable estimates -- the infinite odds ratios that silently vanished from
# Fig1 panel H. Dropped features are named in a message rather than disappearing.
fit_feature_enrichment <- function(d, outcome, features, min_carriers = 10) {
  empty <- tibble::tibble(feature = character(), odds_ratio = numeric(),
                          conf.low = numeric(), conf.high = numeric(),
                          p.value = numeric(), n_carriers = integer())
  y <- as.logical(d[[outcome]])
  keep <- character(0); dropped <- character(0)
  for (f in features) {
    x <- as.logical(d[[f]])
    tab <- table(x, y)
    if (nrow(tab) == 2L && ncol(tab) == 2L && min(tab) >= min_carriers)
      keep <- c(keep, f) else dropped <- c(dropped, f)
  }
  if (length(dropped))
    message("  ", length(dropped), " feature(s) below the ", min_carriers,
            "-carrier floor, not tested: ", paste(dropped, collapse = ", "))
  if (!length(keep)) return(empty)

  # A feature clearing the carrier floor can still be a hard failure: e.g.
  # complete/quasi-separation makes stats::confint() error ("need at least
  # two non-NA values to interpolate") rather than return an unstable-but-
  # finite CI. Fitting each surviving feature is therefore wrapped in its own
  # tryCatch so ONE such feature cannot abort every other feature in `keep`
  # (an unguarded error inside map_dfr would). On failure the feature is
  # dropped and named in a message -- worded differently from the below-floor
  # case above, so the two reasons a feature is missing from the result are
  # distinguishable in the log -- rather than silently returning a row of NAs,
  # which would render cleanly while being wrong.
  purrr::map_dfr(keep, function(f) {
    tryCatch({
      form <- stats::as.formula(
        sprintf("%s ~ `%s` + PIP_qtl + log10(cs_size_med)", outcome, f))
      fit <- stats::glm(form, data = d, family = stats::binomial())
      ci  <- suppressMessages(stats::confint(fit))
      co  <- summary(fit)$coefficients
      # R's model.matrix names a logical predictor's coefficient
      # `<variable name>TRUE`, with any backticks used to make an
      # irregular name syntactic in the formula stripped by the parser --
      # never `` `<name>`TRUE ``. Match that exact name; a substring
      # grep() fallback previously "worked" only by the accident of the
      # feature's row preceding the covariates' rows in this formula's
      # term order, which is untested and not a real invariant.
      row <- paste0(f, "TRUE")
      if (!row %in% rownames(co))
        stop("coefficient row '", row, "' not found among: ",
             paste(rownames(co), collapse = ", "))
      tibble::tibble(
        feature    = f,
        odds_ratio = unname(exp(stats::coef(fit)[row])),
        conf.low   = unname(exp(ci[row, 1])),
        conf.high  = unname(exp(ci[row, 2])),
        p.value    = unname(co[row, 4]),
        n_carriers = sum(as.logical(d[[f]]), na.rm = TRUE))
    }, error = function(e) {
      message("  feature '", f, "' failed to fit/profile, not tested: ",
              conditionMessage(e))
      empty
    })
  })
}

coloc_variant_level <- function(refresh = FALSE)
  build_or_cache("coloc_variant_level", refresh = refresh, {
    af <- .variant_feature_map(susie_finemapping())

    variants <- load_input("enloc_variant", reader = LoadEnlocTable, cols = ENLOC_VARIANT_COLS) %>%
      .label_by_pair("SCP", "variant level") %>%
      mutate(Gene = str_remove(Signal, ":.*"))

    # relationship = "many-to-one" asserts SNP is unique in `af`. .variant_feature_map()
    # already enforces that, so this is a belt-and-braces guard: without it a
    # duplicated key would silently multiply rows (dplyr's default only warns)
    # and build_or_cache would cache the inflated result on this very run.
    joined <- tryCatch(
      left_join(variants, af, by = c("SNP" = "variant"), relationship = "many-to-one"),
      error = function(e)
        stop("SuSiE-derived allele-frequency table has duplicate variant keys, ",
             "violating the expected many-to-one variant-level join: ",
             conditionMessage(e), call. = FALSE)
    )

    # Coverage is no longer guaranteed the way it was with the genome-wide
    # annotation table: an enloc variant absent from every SuSiE table now gets
    # no AF and is dropped downstream by the max_af NA filters in panels C/D.
    # Report it rather than let the sample size shrink silently.
    n_missing <- sum(is.na(joined[[AF_POP_COLS[1]]]))
    if (n_missing)
      message("  note: ", n_missing, " of ", nrow(joined), " variant-level rows ",
              "have no allele frequency in the SuSiE tables (dropped downstream).")

    joined
  })

# =============================================================================
# Collapsing colocalized credible sets into components
# -----------------------------------------------------------------------------
# For each assay pair, keep colocalized signals whose two sides belong to the
# SAME gene, and treat each as an undirected edge between two credible sets.
# Connected components are then signals shared across molecular layers.
#
# A Signal's second side can list several credible sets joined by "_", so it is
# split before matching. eQTL-sQTL sQTL ids embed a chromosome, hence its
# different split pattern.
# =============================================================================
.clean_pair_edges <- function(cs_colocs, pair, split_pattern) {
  cs_colocs %>%
    filter(Coloc == TRUE, type == pair) %>%
    separate(Signal, into = c("sig_1", "sig_2"), sep = "\\(@\\)",
             fill = "right", extra = "merge") %>%
    mutate(sig_2_list = str_split(sig_2, split_pattern)) %>%
    unnest(sig_2_list) %>%
    mutate(gene_1 = str_extract(sig_1,      "^ENSG[0-9]+"),
           gene_2 = str_extract(sig_2_list, "^ENSG[0-9]+")) %>%
    filter(gene_1 == gene_2) %>%
    select(cs_1, sig_2_list, RCP, LCP, type, gene_1)
}

# Pure function of its argument so it can be tested without a manifest.
coloc_cs_edges <- function(cs_colocs) {
  bind_rows(
    .clean_pair_edges(cs_colocs, "eQTL-sQTL", "_(?=ENSG[0-9]+:chr)"),
    .clean_pair_edges(cs_colocs, "eQTL-pQTL", "_(?=ENSG[0-9]+:)"),
    .clean_pair_edges(cs_colocs, "sQTL-pQTL", "_(?=ENSG[0-9]+:)")
  ) %>%
    mutate(cs_2 = str_remove(sig_2_list, "^ENSG.*:")) %>%
    select(gene = gene_1, RCP, LCP, type, cs_1, cs_2)
}

coloc_cs_components <- function(refresh = FALSE)
  build_or_cache("coloc_cs_components", refresh = refresh, {
    edge_list <- coloc_cs_edges(coloc_cs_level()) %>%
      transmute(from = cs_1, to = cs_2) %>%
      filter(!is.na(from), !is.na(to)) %>%
      distinct(from, to)

    g    <- igraph::graph_from_data_frame(edge_list, directed = FALSE)
    memb <- igraph::components(g)$membership

    components_df <- tibble(cs_id           = names(memb),
                            coloc_component = as.integer(memb))

    susie_finemapping()$cs_list %>%
      dplyr::rename(Gene = molecular_trait_id) %>%
      left_join(components_df, by = "cs_id") %>%
      mutate(has_coloc = !is.na(coloc_component))
  })

# =============================================================================
# Variant-level model input
# -----------------------------------------------------------------------------
# Annotation allele frequencies are reported against the ALT allele, so they are
# folded to minor-allele frequency before the per-variant maximum across
# ancestry groups is taken. log_maf is standardised across ALL pairs pooled (as
# in the original) and the per-pair models are fit on that common scale.
# =============================================================================
coloc_model_input <- function(refresh = FALSE)
  build_or_cache("coloc_model_input", refresh = refresh, {
    v  <- coloc_variant_level()
    fm <- susie_finemapping()

    missing_af <- setdiff(AF_POP_COLS, names(v))
    if (length(missing_af))
      stop("SuSiE fine-mapping tables are missing allele-frequency columns: ",
           paste(missing_af, collapse = ", "))

    folded <- v %>%
      select(all_of(AF_POP_COLS)) %>%
      mutate(across(everything(), ~ ifelse(. > .5, 1 - ., .)))

    cleaned <- v %>%
      mutate(max_af  = do.call(pmax, as.list(folded)),
             rare    = max_af < .01,
             log_maf = as.numeric(scale(-log10(max_af)))) %>%
      filter(Coloc == TRUE)

    # One row per (variant, gene) for the posterior join. Splicing
    # molecular_trait_ids are INTRON-level and collapsed to a gene above, so a
    # variant sitting in credible sets for two introns of the same gene yields
    # two rows -- 13,950 such pairs in the sQTL table, plus a handful in eQTL and
    # pQTL. Panels C/D regress cross-assay EFFECT SIZES, so the strongest signal
    # is the faithful single representative of what that variant does to that
    # gene; keeping both would double-count it and inflate the regression's N,
    # which is precisely what the many-to-one guard below refuses to allow.
    #
    # Collapsing happens HERE and not in susie_finemapping(): the full tables
    # still feed cs_list / gene_map / the credible-set graph behind panels A and
    # B, and dropping rows there would delete whole credible sets from them.
    #
    # arrange() + distinct() keeps the first row per key, which after sorting on
    # descending |posterior_mean| is the strongest; it is far cheaper than a
    # grouped slice_max() over ~840k groups and is deterministic.
    posteriors <- function(d)
      d %>%
        select(variant, molecular_trait_id, posterior_mean, posterior_sd) %>%
        arrange(desc(abs(posterior_mean))) %>%
        distinct(variant, molecular_trait_id, .keep_all = TRUE)

    # relationship = "many-to-one" asserts (variant, gene) is unique on the
    # fine-mapping side -- nothing upstream guarantees this, and a variant
    # sitting in two credible sets for the same gene (a real SuSiE outcome) is
    # exactly the case that would otherwise fan the join out and silently
    # inflate the regression's sample size, double-counting that variant.
    # Wrap the bare dplyr error so it names which fine-mapping table (eQTL,
    # sQTL or pQTL) holds the duplicate and shows an offending (variant, gene)
    # pair, mirroring the annotation join above in coloc_variant_level().
    join_posteriors <- function(d, fm_table, table_name) {
      p <- posteriors(fm_table)
      tryCatch(
        left_join(d, p, by = c("SNP" = "variant", "Gene" = "molecular_trait_id"),
                  relationship = "many-to-one"),
        error = function(e) {
          key <- p[c("variant", "molecular_trait_id")]
          dup <- p[duplicated(key) | duplicated(key, fromLast = TRUE), , drop = FALSE]
          offender <- if (nrow(dup))
            sprintf("(variant = '%s', gene = '%s')", dup$variant[1], dup$molecular_trait_id[1])
          else
            "<offending pair not found>"
          stop(table_name, " fine-mapping table has duplicate (variant, gene) keys, ",
               "violating the expected many-to-one posterior join -- e.g. ", offender,
               ": ", conditionMessage(e), call. = FALSE)
        }
      )
    }

    pair_input <- function(pair, left, left_name, right, right_name) {
      # Each join_posteriors() call must be forced to a concrete data frame
      # (via an explicit <- binding) before the next one starts. Chaining the
      # two calls through %>% instead would make the outer call's `d` argument
      # a lazy promise wrapping the inner join_posteriors() call; forcing it
      # inside the outer tryCatch's monitored left_join() would let the outer
      # handler catch and mis-attribute the inner join's error to the
      # right-hand table's name.
      d0 <- cleaned %>% filter(type == pair)
      d1 <- join_posteriors(d0, left, left_name)
      join_posteriors(d1, right, right_name)
    }

    list(
      `eQTL-sQTL` = pair_input("eQTL-sQTL", fm$eQTL, "eQTL", fm$sQTL, "sQTL"),
      `eQTL-pQTL` = pair_input("eQTL-pQTL", fm$eQTL, "eQTL", fm$pQTL, "pQTL"),
      `sQTL-pQTL` = pair_input("sQTL-pQTL", fm$sQTL, "sQTL", fm$pQTL, "pQTL")
    )
  })

# =============================================================================
# Panel E — colocalization feature enrichment
# =============================================================================
# A pair is flagged underpowered for an analysis when fewer than this fraction
# of its candidate features clear the carrier floor.
UNDERPOWERED_FRAC <- 0.5

# Pure orchestration over an already-built feature matrix `v` (the output of
# build_feature_matrix()): cross both analyses with every observed pair, fit
# each combination, correct p-values WITHIN (analysis, pair), and flag pairs
# where fewer than UNDERPOWERED_FRAC of candidate features survive.
#
#   A  what makes a variant shared    coloc vs non-coloc, within pair
#   B  what makes a shared variant rare   rare vs common, within Coloc == TRUE
#
# Analysis B is attempted for every pair; where a pair has too few rare
# colocalizing variants the carrier floor empties it and the underpowered flag
# fires on its own. The criterion decides -- nothing is excluded a priori.
#
# `min_carriers` is a parameter HERE, not on coloc_feature_enrichment(), so
# tests can lower the floor against a small synthetic matrix without ever
# threading a non-default floor through the memoized builder:
# build_or_cache() keys its cache file on NAME ONLY, so a non-default
# min_carriers passed to the public function could silently read back a
# stale cache built under a different floor.
.enrich_pairs <- function(v, min_carriers = 10) {
  feats <- feature_names(v)
  message("  ", length(feats), " candidate features")

  specs <- tidyr::expand_grid(analysis = c("A", "B"), pair = unique(v$type))

  # The documented shape, kept even when NO (analysis, pair) combination
  # survives the carrier floor (a real outcome on a small fixture/cohort,
  # not an error): bind_rows() against this schema is what keeps the
  # result a well-formed zero-row tibble instead of the columnless 0x0
  # frame purrr::pmap_dfr() returns when every iteration answers NULL.
  empty <- tibble::tibble(
    analysis = character(), pair = character(), feature = character(),
    odds_ratio = numeric(), conf.low = numeric(), conf.high = numeric(),
    p.value = numeric(), padj = numeric(), n_carriers = integer(),
    flagged = logical())

  out <- purrr::pmap_dfr(specs, function(analysis, pair) {
    d <- v %>% filter(type == pair)
    if (analysis == "A") {
      d$.outcome <- as.logical(d$Coloc)
    } else {
      d <- d %>% filter(Coloc == TRUE)
      d$.outcome <- as.logical(d$rare)
    }
    # `rare` is a tested FEATURE in A but the OUTCOME in B. Regressing a
    # variable on itself is meaningless, so it is dropped from B only.
    use <- if (analysis == "A") feats else setdiff(feats, "rare")
    message("  analysis ", analysis, " / ", pair, " (n = ", nrow(d), ")")
    res <- fit_feature_enrichment(d, ".outcome", use, min_carriers = min_carriers)
    # Coverage assertion: a shrinking feature set must not look like a clean run.
    n_expected <- length(use)
    if (nrow(res) < n_expected)
      message("    ", nrow(res), "/", n_expected, " features tested")
    if (!nrow(res)) return(NULL)
    res %>%
      mutate(analysis = analysis, pair = pair,
             padj = stats::p.adjust(p.value, method = "fdr"),
             flagged = nrow(res) < UNDERPOWERED_FRAC * n_expected)
  })
  dplyr::bind_rows(empty, out)
}

# Memoized wrapper: build the annotated feature matrix from the fine-mapping /
# colocalization tables and hand it to .enrich_pairs() at the default floor.
coloc_feature_enrichment <- function(refresh = FALSE)
  build_or_cache("coloc_feature_enrichment", refresh = refresh, {
    vl <- coloc_variant_level()
    # Fold ALT-allele frequencies to MAF, then take the per-variant maximum
    # across ancestries. Two steps, matching coloc_model_input() exactly --
    # inlining across() into do.call(pmax, ...) does not evaluate the way it
    # reads.
    folded <- vl %>%
      select(all_of(AF_POP_COLS)) %>%
      mutate(across(everything(), ~ ifelse(. > .5, 1 - ., .)))
    v <- vl %>%
      mutate(max_af = do.call(pmax, as.list(folded))) %>%
      build_feature_matrix()

    .enrich_pairs(v)
  })

# Pair-contrast scatter: adjusted log-OR in eQTL-pQTL against eQTL-sQTL.
# Points on the diagonal are generic shared-signal markers; off-diagonal ones
# distinguish protein- from splicing-colocalization. sQTL-pQTL is structurally
# absent -- this plot has two axes -- and appears in the supplement instead.
build_panelE <- function(enr) {
  pairs_needed <- c("eQTL-pQTL", "eQTL-sQTL")

  wide <- enr %>%
    filter(analysis == "A", pair %in% pairs_needed) %>%
    mutate(log_or = log2(odds_ratio), sig = padj < 0.05) %>%
    select(analysis, feature, pair, log_or, sig) %>%
    tidyr::pivot_wider(names_from = pair,
                       values_from = c(log_or, sig))

  # pivot_wider() only creates a `log_or_<pair>`/`sig_<pair>` column for a
  # pair value actually present in the (already analysis-A-only) input.
  # coloc_feature_enrichment() can legitimately return zero rows (nothing
  # clears the carrier floor -- exactly what the small test fixture does),
  # or rows for only one of the two required pairs, in which case one or
  # both expected columns never appear. Add any missing one as all-NA so the
  # filter/aes() below have something to reference instead of erroring on a
  # column that doesn't exist -- an all-NA column drops every row from the
  # !is.na() filter next, which is the correct behaviour: no feature can be
  # plotted without both pairs.
  needed_cols <- paste0(c("log_or_", "sig_"), rep(pairs_needed, each = 2))
  for (col in setdiff(needed_cols, names(wide))) wide[[col]] <- NA

  wide <- wide %>%
    filter(!is.na(`log_or_eQTL-pQTL`), !is.na(`log_or_eQTL-sQTL`))

  ggplot(wide, aes(`log_or_eQTL-pQTL`, `log_or_eQTL-sQTL`)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey85") +
    geom_point(aes(fill = `sig_eQTL-pQTL` & `sig_eQTL-sQTL`),
               shape = 21, size = 3.4, stroke = 0.4, colour = "grey25") +
    ggrepel::geom_text_repel(aes(label = feature), size = 3, max.overlaps = 12,
                             segment.colour = "grey70") +
    scale_fill_manual(values = c(`TRUE` = "#364B63", `FALSE` = "grey85"),
                      name = "FDR < 0.05 in both") +
    labs(x = "Adjusted log2 OR — eQTL-pQTL",
         y = "Adjusted log2 OR — eQTL-sQTL") +
    theme(legend.position = "bottom")
}

# One forest-plot row for the FigS_ColocEnrichment supplement (one analysis
# half: A or B, already filtered to a single `analysis` value by the caller).
#
# `d` can legitimately have zero rows -- not only when coloc_feature_enrichment()
# returns fully empty (e.g. the 9-row test fixture, too small for anything to
# clear the carrier floor: see .enrich_pairs() above), but also when the OTHER
# analysis half has real rows and this one does not (a pair can clear the
# floor for "what makes a variant shared" and not for "what makes a shared
# variant rare", or vice versa). `facet_wrap(~pair)` on a 0-row frame errors
# at render time ("Faceting variables must have at least one value") because
# the `pair` factor carries no levels to facet on -- there is nothing to draw.
#
# A placeholder panel that names the analysis is rendered instead of silently
# dropping the row: a missing half of a supplement reads as a rendering bug,
# not as "this analysis produced no results" -- the two must be
# distinguishable at a glance. Extracted here (out of the standalone
# Fig3/supp_coloc_enrichment.R) specifically so it is unit-testable; a helper
# buried in a script that CI never runs is how this went uncaught.
enrichment_forest_panel <- function(d, title) {
  if (nrow(d) == 0L) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0,
                 label = paste0("No feature cleared the carrier floor\nfor analysis ", title),
                 size = 3.2, colour = "grey40", family = FONT) +
        theme_void(base_family = FONT) +
        labs(title = title) +
        theme(plot.title = element_text(size = SZ_STRIP, face = "bold"))
    )
  }
  ggplot(d, aes(x = odds_ratio, y = reorder(feature, odds_ratio),
                colour = flagged)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.35) +
    facet_wrap(~pair, nrow = 1) +
    scale_x_continuous(trans = "log10") +
    # Grey = the pair was flagged underpowered for this analysis.
    scale_colour_manual(values = c(`FALSE` = "#364B63", `TRUE` = "grey70"),
                        guide = "none") +
    labs(x = "Adjusted odds ratio", y = NULL, title = title) +
    theme(plot.title = element_text(size = SZ_STRIP, face = "bold"))
}
