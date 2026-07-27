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

# Columns every SuSiE table must supply. AF_POP_COLS are carried through so
# .variant_af_map() can source allele frequencies here instead of from a
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
      select(all_of(SUSIE_BASE_COLS), gene_name, all_of(AF_POP_COLS)) %>%
      mutate(molecular_trait_id = str_remove(molecular_trait_id, "\\..*"),
             variant = str_replace(variant, "chrchr", "chr")) %>%
      .add_cs_stats()

    sQTL <- .read_susie_table("susie_sqtl") %>%
      mutate(molecular_trait_id = str_extract(molecular_trait_id, "ENSG[0-9]+")) %>%
      select(all_of(SUSIE_BASE_COLS), all_of(AF_POP_COLS)) %>%
      .add_cs_stats()

    pQTL <- .read_susie_table("susie_pqtl") %>%
      mutate(molecular_trait_id = str_extract(molecular_trait_id, "ENSG[0-9]+")) %>%
      select(all_of(SUSIE_BASE_COLS), all_of(AF_POP_COLS)) %>%
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

# Allele frequency is a property of the VARIANT, not of an assay, so the three
# SuSiE tables are pooled and reduced to one row per variant. Rows whose AF
# columns are entirely NA carry no information and are dropped before the
# uniqueness check, so a table that simply lacks a variant cannot manufacture a
# spurious conflict with one that has it.
.variant_af_map <- function(fm) {
  af <- bind_rows(
      fm$eQTL %>% select(variant, all_of(AF_POP_COLS)),
      fm$sQTL %>% select(variant, all_of(AF_POP_COLS)),
      fm$pQTL %>% select(variant, all_of(AF_POP_COLS))
    ) %>%
    mutate(across(all_of(AF_POP_COLS), as.numeric),
           # The eQTL reader already repairs "chrchr"; normalise all three here
           # so a stray prefix cannot split one variant into two AF rows.
           variant = str_replace(variant, "chrchr", "chr"))
  af <- af[rowSums(!is.na(af[AF_POP_COLS])) > 0, , drop = FALSE] %>% distinct()

  # One AF row per variant, or the variant-level join below fans out. A
  # duplicate here means two SuSiE tables disagree about the same variant's
  # allele frequencies, which is a genuine upstream inconsistency rather than
  # something to silently pick a winner for.
  dup <- unique(af$variant[duplicated(af$variant)])
  if (length(dup))
    stop(length(dup), " variant(s) carry CONFLICTING per-ancestry allele ",
         "frequencies across the SuSiE fine-mapping tables -- e.g. '", dup[1],
         "'. Allele frequency is a property of the variant alone; resolve the ",
         "disagreement upstream.", call. = FALSE)
  af
}

coloc_variant_level <- function(refresh = FALSE)
  build_or_cache("coloc_variant_level", refresh = refresh, {
    af <- .variant_af_map(susie_finemapping())

    variants <- load_input("enloc_variant", reader = LoadEnlocTable, cols = ENLOC_VARIANT_COLS) %>%
      .label_by_pair("SCP", "variant level") %>%
      mutate(Gene = str_remove(Signal, ":.*"))

    # relationship = "many-to-one" asserts SNP is unique in `af`. .variant_af_map()
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

    posteriors <- function(d)
      d %>% select(variant, molecular_trait_id, posterior_mean, posterior_sd)

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
