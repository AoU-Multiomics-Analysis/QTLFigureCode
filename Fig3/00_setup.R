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
OUTDIR   <- Sys.getenv("FIG3_OUTDIR",   ".")            # final Figure3.pdf/png
MANIFEST <- Sys.getenv("FIG3_MANIFEST", "fig3/inputs.tsv")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OBJ_DIR,  showWarnings = FALSE, recursive = TRUE)

source("common/manifest.R")

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

blank_plot <- ggplot() + theme_void()

# The three assay pairs, in the order used throughout the figure.
COLOC_PAIRS <- c("eQTL-pQTL", "sQTL-pQTL", "eQTL-sQTL")

# Panel A fill — number of colocalized molecular layers.
LAYER_FILL <- c("0" = "#D9D9D9", "1" = "#A9C1D9", "2" = "#7B9AB8", "3" = "#577590")

# Panel B fill — number of colocalized signals.
SIGNAL_FILL <- c("0" = "#E9DCDC", "1" = "#D8BBBB", "2" = "#C89A9A", "3" = "#B37979",
                 "4" = "#9E5E5E", "5" = "#8A4D4D", "6" = "#763F3F", "7" = "#643232")

# Per-ancestry allele-frequency columns of the AoU annotation table, used to
# derive a per-variant maximum MAF. Named explicitly so a schema change fails
# loudly rather than silently changing which columns are folded.
AF_POP_COLS <- c("gvs_afr_af", "gvs_amr_af", "gvs_eas_af", "gvs_eur_af",
                 "gvs_mid_af", "gvs_sas_af", "gvs_oth_af")

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
  # ever tripping the ncol() guard below. Check lengths before rbind instead.
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
# Fig3's panels share a heavy upstream: the AoU annotation table (the largest
# input in the repo) and three SuSiE tables. Under strict per-panel loading the
# annotation join would run three times per build. Each builder below memoizes
# its result to DATA_DIR/<name>.rds, so a panel still calls only the builders it
# needs while the expensive work happens once.
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

susie_finemapping <- function(refresh = FALSE)
  build_or_cache("susie_finemapping", refresh = refresh, {
    eQTL <- fread(input_path("susie_eqtl")) %>%
      select(molecular_trait_id, variant, pip, cs_index, posterior_mean, cs_id,
             posterior_sd, gene_name) %>%
      mutate(molecular_trait_id = str_remove(molecular_trait_id, "\\..*"),
             variant = str_replace(variant, "chrchr", "chr")) %>%
      .add_cs_stats()

    sQTL <- fread(input_path("susie_sqtl")) %>%
      mutate(molecular_trait_id = str_extract(molecular_trait_id, "ENSG[0-9]+")) %>%
      select(molecular_trait_id, variant, pip, cs_index, posterior_mean, cs_id,
             posterior_sd) %>%
      .add_cs_stats()

    pQTL <- fread(input_path("susie_pqtl")) %>%
      filter(group == "COMB") %>%
      mutate(molecular_trait_id = str_extract(molecular_trait_id, "ENSG[0-9]+")) %>%
      select(molecular_trait_id, variant, pip, cs_index, posterior_mean, cs_id,
             posterior_sd) %>%
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
