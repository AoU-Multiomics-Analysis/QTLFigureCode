# Figure 4 Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the flat 774-line `Fig4/Fig4.R` (GWAS-trait colocalization) with a modular, manifest-driven, Fig3-style pipeline, fixing the recorded defects — chiefly the twice-computed `ColocAdjustedConditionalTWAS` (whose first copy references an undefined variable) and the dead `GeneBayes` read.

**Architecture:** `Fig4/00_setup.R` sources the shared `common/manifest.R`, then adds Fig4's paths, theme, palettes, folded helper functions, and six memoized derived-data builders that cache to `fig4/data/`. Each panel script sources setup, calls only the builders it needs, and caches one ggplot object to `fig4/objects/`. `assemble.R` composes them; `run_all.R` runs everything.

**Tech Stack:** R 4.3, tidyverse, data.table, patchwork, ggh4x, MatchIt, arrow, broom, showtext/sysfonts, paletteer. Tests are plain `Rscript` files using the existing `tests/helpers.R` (`ok`/`expect_error`).

**Spec:** `docs/superpowers/specs/2026-07-24-fig4-cleanup-design.md`

## Global Constraints

- **Branch:** `fig4-cleanup` (off `main`; spec committed as `4bff2ae`).
- **Run every script from the project root**, referring to it as lowercase `fig4/...` (macOS is case-insensitive; matches Fig1/Fig2/Fig3).
- **Never commit inputs or figure outputs** — *All of Us* controlled-access data.
- **`showtext` DPI must equal every `ggsave()` DPI.** Fig4 uses `DPI <- 400`.
- **`element_text()` sizes are pt; `geom_*`/`annotate()` sizes are mm.** Never reuse a constant across the two.
- **Tagging is `assemble.R`'s job.** No panel script may call `labs(tag = ...)`.
- **Reserved manifest columns are exactly `key`, `local_subpath`, `gs_uri`.**
- **The derived-data cache is keyed by name only**, not by input content.
- Commit after every task. Never `git add -A` (controlled-access files may be untracked in the tree).

---

## File Structure

| File | Responsibility |
|---|---|
| `Fig4/00_setup.R` | **Create.** Libraries, `FIG4_*` paths, theme, `SZ_*`, palettes, folded helpers, `build_or_cache`, six builders. |
| `Fig4/inputs.tsv` | **Create.** 6 keys / 6 rows. |
| `Fig4/panels/panelA.R`–`panelD.R` | **Create.** One ggplot object each. |
| `Fig4/assemble.R` | **Create.** Compose → `Figure4.pdf` / `Figure4.png`. |
| `Fig4/run_all.R` | **Create.** Panels → assemble. |
| `Fig4/README.md`, `Fig4/NOTES.md`, `Fig4/.gitignore` | **Create.** |
| `Fig4/Fig4.R` | **Commit unchanged (Task 1), then `git rm` (Task 9).** |
| `.gitignore` | **Modify.** Add Fig4 input dirs + re-apply `/data/` `/susie_local/`. |
| `tests/fixtures.R` | **Modify.** Add `make_fig4_fixture()` / `use_fig4_fixture()`. |
| `tests/test_fig4_helpers.R` | **Create.** Pure-transform tests (Tasks 2–5). |
| `tests/test_fig4_assembly.R` | **Create.** Rendered A/B/C/D tag sequence (Task 8). |

---

## Task 1: Preserve original, setup skeleton, manifest, gitignore

**Files:**
- Commit then later remove: `Fig4/Fig4.R`
- Create: `Fig4/00_setup.R` (paths/libs/theme/palettes/`build_or_cache` only; helpers added in later tasks)
- Create: `Fig4/inputs.tsv`, `Fig4/.gitignore`
- Modify: `.gitignore`

**Interfaces:**
- Consumes: `common/manifest.R` (`input_path`, `input_files`, `load_input`).
- Produces: globals `INPUT`, `DATA_DIR`, `OBJ_DIR`, `OUTDIR`, `MANIFEST`, `FONT`, `DPI`, `SZ_*`, `GROUP_LABELS`, `GROUP_COLORS`, `RATE_COLUMNS`, `BLOOD_TRAITS`; function `build_or_cache(name, expr, refresh = FALSE)`.

- [ ] **Step 1: Preserve the original in history**

`Fig4/Fig4.R` is currently untracked. Commit it unchanged so it stays recoverable:

```bash
git add Fig4/Fig4.R
git commit -m "Add the original Fig4.R before cleanup

Fig4/Fig4.R was never committed. Recording it unchanged so the original analysis
stays recoverable from history once the cleanup replaces it, and so the
restructuring is reviewable as a real diff."
```

- [ ] **Step 2: Create the manifest**

Create `Fig4/inputs.tsv` — **tab-separated**, 1 header + 6 rows:

```
key	local_subpath	gs_uri
coloc_gene	high_level_coloc_outputs/harmonized/harmonized_coloc.gene.tsv.gz	-
coloc_rate	high_level_coloc_outputs/final_summaries/coloc_rate_by_trait.tsv	-
r2_twas	R2_TWAS_combined.tsv	-
conditional_twas	ConditionalTWAS.cleaned.tsv.gz	-
drug_interactions	reference/interactions.tsv	-
clingen	reference/evidence_clingen	-
```

(The `clingen` row points at the parquet *directory*; the builder globs `*.parquet` inside it.)

- [ ] **Step 3: Create the setup skeleton**

Create `Fig4/00_setup.R`:

```r
# =============================================================================
# Figure 4 — shared setup (sourced by every panel / assembly script).
# Run all scripts from the PROJECT ROOT. Each script starts with:
#   source("fig4/00_setup.R")
#
# Figure 4 covers colocalization of molecular QTLs with GWAS traits and the
# downstream rare-variant / conditional-TWAS analyses.
#
# See docs/superpowers/specs/2026-07-24-fig4-cleanup-design.md
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggh4x)
  library(MatchIt)
  library(arrow)
  library(broom)
  library(scales)
  library(showtext)
  library(sysfonts)
  library(paletteer)
})

# --- Paths (override via env vars if needed) --------------------------------
INPUT    <- Sys.getenv("FIG4_INPUT",    ".")            # staging root (cwd)
DATA_DIR <- Sys.getenv("FIG4_DATA",     "fig4/data")    # memoized derived tables
OBJ_DIR  <- Sys.getenv("FIG4_OBJ",      "fig4/objects") # per-panel ggplot objects
OUTDIR   <- Sys.getenv("FIG4_OUTDIR",   ".")            # final Figure4.pdf/png
MANIFEST <- Sys.getenv("FIG4_MANIFEST", "fig4/inputs.tsv")
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
SZ_AXIS_TITLE <- 14; SZ_AXIS_TEXT <- 12; SZ_STRIP <- 12
SZ_LEG_TEXT   <- 10; SZ_LEG_TITLE <- 11; SZ_TAG   <- 14
SZ_TITLE      <- 12

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

# --- Trait-category palette / labels (used by Panel A) ----------------------
GROUP_LABELS <- c(
  autoimmune = "AutoImm", inflammatory_biomarker = "Inflam",
  allergic_respiratory = "Allerg", dermatological = "Derm",
  cardiometabolic = "CardioM", lipid_metabolic = "Lipid", lipid = "Lipid",
  metabolic = "Metabolic", clinical_chemistry = "ClinChem",
  hematological = "Heme", anthropometric = "Anthro", renal = "Renal",
  renal_urological = "Renal/Urol", gastrointestinal = "GI",
  reproductive_endocrine = "ReproEnd", musculoskeletal = "MSK",
  cancer = "Cancer", neurodegenerative = "Neurodeg", neurological = "Neurol",
  neuropsychiatric = "NeuroPsy", ophthalmological = "Ophth"
)
GROUP_COLORS <- c(
  autoimmune = "#D62828", inflammatory_biomarker = "#F77F00",
  allergic_respiratory = "#2A9D8F", dermatological = "#E76F51",
  cardiometabolic = "#277DA1", lipid_metabolic = "#6A994E", lipid = "#6A994E",
  metabolic = "#80B918", clinical_chemistry = "#F4A261",
  hematological = "#A23B72", anthropometric = "#E9C46A", renal = "#457B9D",
  renal_urological = "#8D99AE", gastrointestinal = "#BC6C25",
  reproductive_endocrine = "#E5989B", musculoskeletal = "#7F5539",
  cancer = "#9B5DE5", neurodegenerative = "#4361EE", neurological = "#3A0CA3",
  neuropsychiatric = "#B5179E", ophthalmological = "#48CAE4"
)
RATE_COLUMNS <- c(stringent = "pct_stringent_union", lenient = "pct_lenient_union")

# --- Blood traits used by Panel C -------------------------------------------
BLOOD_TRAITS <- c("eosinophil", "hematocrit", "lymphocyte", "MCH", "MCHC",
                  "monocyte", "neutrophil", "platelet_count", "RDW")

# =============================================================================
# Derived-data cache (identical contract to Fig3): lazy `expr`, cache to
# DATA_DIR/<name>.rds, keyed by NAME only. Clear DATA_DIR after changing an input.
# =============================================================================
build_or_cache <- function(name, expr, refresh = FALSE) {
  f <- file.path(DATA_DIR, paste0(name, ".rds"))
  if (!refresh && file.exists(f)) return(readRDS(f))
  message("building ", name, " ...")
  val <- force(expr)
  saveRDS(val, f)
  val
}
```

- [ ] **Step 4: Create Fig4/.gitignore**

```
# Per-panel ggplot object cache (regenerated by the pipeline)
objects/

# Memoized derived tables (regenerated by the builders in 00_setup.R)
data/

# macOS
.DS_Store
```

- [ ] **Step 5: Extend the root .gitignore**

Add to the controlled-access inputs block of `.gitignore`:

```
/high_level_coloc_outputs/
/R2_TWAS_combined.tsv
/ConditionalTWAS.cleaned.tsv.gz
/data/
/susie_local/
```

Add to the "Figure outputs" block:

```
/Figure4.pdf
/Figure4.png
```

(`/reference/` is already present.)

- [ ] **Step 6: Verify setup sources and the manifest resolves**

Run:
```bash
Rscript -e '
  source("fig4/00_setup.R")
  stopifnot(nrow(input_manifest()) == 6L)
  stopifnot(length(unique(input_manifest()$key)) == 6L)
  stopifnot(is.function(build_or_cache), DPI == 400)
  cat("Fig4 setup OK: 6 rows, 6 keys\n")'
```
Expected: `Fig4 setup OK: 6 rows, 6 keys`

- [ ] **Step 7: Verify staged paths and outputs are git-ignored**

Run:
```bash
git check-ignore -v high_level_coloc_outputs/x.tsv R2_TWAS_combined.tsv \
  ConditionalTWAS.cleaned.tsv.gz reference/interactions.tsv Figure4.pdf \
  Fig4/objects/A.rds Fig4/data/gene_sets.rds data/x susie_local/x
git check-ignore Fig2/objects/A.rds && echo "REGRESSION" || echo "Fig2 objects still tracked - OK"
```
Expected: a match line per path in the first command; `Fig2 objects still tracked - OK`

- [ ] **Step 8: Commit**

```bash
git add Fig4/00_setup.R Fig4/inputs.tsv Fig4/.gitignore .gitignore
git commit -m "Add Fig4 setup skeleton and input manifest

Paths, font/DPI, theme, trait-category palette/labels and the derived-data cache;
the manifest declares 6 keys. Inputs are recorded local-only (only the coloc
tar.gz has a known gs:// origin, documented in the README as a manual staging
step). Also re-applies the /data/ and /susie_local/ ignores as a safety measure."
```

---

## Task 2: Coloc-plot helpers and the `coloc_by_trait` builder (Panel A upstream)

**Files:**
- Modify: `Fig4/00_setup.R` (append)
- Create: `tests/fixtures.R` additions, `tests/test_fig4_helpers.R`

**Interfaces:**
- Consumes: `build_or_cache`, `input_path`, `GROUP_COLORS`, `GROUP_LABELS`, `RATE_COLUMNS`.
- Produces: `read_any(path, ...)`, `require_columns`, `as_logical_flag`, `complete_category_colors`, `blend_with_white`, `make_fill`, `reorder_within`, `scale_y_reordered`, `build_plot_data(gene, coloc_rate, min_coloc_genes = 1L, rate_type = "stringent")`, `build_empty_plot()`, `build_plot(combined)`, `coloc_by_trait(refresh = FALSE)`.

- [ ] **Step 1: Add the Fig4 fixture**

Append to `tests/fixtures.R`:

```r
# --- Fig4 mock project -------------------------------------------------------
# Small TSV stand-ins for Fig4's tabular inputs, enough to exercise the pure
# transforms. The ClinGen parquet and gene-set derivation are tested from an
# in-memory frame (see tests/test_fig4_helpers.R), so no parquet is written here.
make_fig4_fixture <- function(root = file.path(tempdir(), "fig4_fixture")) {
  unlink(root, recursive = TRUE)
  for (d in c("high_level_coloc_outputs/harmonized",
              "high_level_coloc_outputs/final_summaries", "reference"))
    dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
  wtsv <- function(df, p) write.table(df, p, sep = "\t", row.names = FALSE, quote = FALSE)

  # coloc gene table: gene / trait / trait_category / any_coloc
  genes <- expand.grid(g = 1:8, trait = c("eosinophil", "LDL"),
                       stringsAsFactors = FALSE)
  gene_tbl <- data.frame(
    gene = paste0("ENSG", sprintf("%08d", genes$g)),
    trait = genes$trait,
    trait_category = ifelse(genes$trait == "eosinophil", "hematological", "lipid"),
    any_coloc = TRUE)
  # fwrite gzips by the .gz extension, so read_any()'s `gzip -dc` path works.
  data.table::fwrite(
    gene_tbl,
    file.path(root, "high_level_coloc_outputs/harmonized/harmonized_coloc.gene.tsv.gz"),
    sep = "\t")

  # coloc rate: one union row per trait
  rate_tbl <- data.frame(
    trait = c("eosinophil", "LDL"), layer = "union",
    pct_stringent_union = c(30, 20), pct_lenient_union = c(45, 33))
  wtsv(rate_tbl, file.path(root, "high_level_coloc_outputs/final_summaries/coloc_rate_by_trait.tsv"))

  list(root = root)
}

use_fig4_fixture <- function(fx) {
  Sys.setenv(FIG4_INPUT    = fx$root,
             FIG4_DATA     = file.path(fx$root, "data"),
             FIG4_OBJ      = file.path(fx$root, "objects"),
             FIG4_OUTDIR   = fx$root)
}
```

- [ ] **Step 2: Write the failing test**

Create `tests/test_fig4_helpers.R`:

```r
# Tests for Fig4's pure helper transforms. Run from the project root:
#   Rscript tests/test_fig4_helpers.R
source("tests/fixtures.R")
fx <- make_fig4_fixture(); use_fig4_fixture(fx)
source("fig4/00_setup.R")
source("tests/helpers.R")

# --- build_plot_data -------------------------------------------------------
gene <- as_tibble(data.table::fread(input_path("coloc_gene")))
rate <- as_tibble(data.table::fread(input_path("coloc_rate")))

combined <- build_plot_data(gene, rate, min_coloc_genes = 1, rate_type = "lenient")
ok("build_plot_data returns one row per trait", nrow(combined) == 2L)
ok("carries the selected union rate",
   all(c("pct_coloc_union", "count") %in% names(combined)))
ok("count is the number of colocalizing genes per trait",
   all(combined$count == 8L))
ok("lenient rate is selected", combined$pct_coloc_union[combined$trait == "eosinophil"] == 45)

expect_error("an unknown rate_type errors",
             build_plot_data(gene, rate, rate_type = "nope"))
expect_error("a non-positive min_coloc_genes errors",
             build_plot_data(gene, rate, min_coloc_genes = 0))

# min_coloc_genes trims traits below the threshold
trimmed <- build_plot_data(gene, rate, min_coloc_genes = 9, rate_type = "lenient")
ok("min_coloc_genes trims traits with too few genes", nrow(trimmed) == 0L)

# as_logical_flag coerces common encodings
ok("as_logical_flag maps yes/1/true",
   identical(as_logical_flag(c("yes", "1", "true", "no")), c(TRUE, TRUE, TRUE, FALSE)))

cat("test_fig4_helpers.R (Task 2): all assertions passed\n")
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `Rscript tests/test_fig4_helpers.R`
Expected: FAIL — `could not find function "build_plot_data"`

- [ ] **Step 4: Write the implementation**

Append to `Fig4/00_setup.R`:

```r
# =============================================================================
# Panel A helpers — coloc-by-trait-category plot
# =============================================================================
read_any <- function(path, ...) {
  args <- list(showProgress = FALSE, ...)
  if (grepl("\\.gz$", path)) args$cmd <- paste("gzip -dc", shQuote(path))
  else                         args$input <- path
  as_tibble(do.call(data.table::fread, args))
}

require_columns <- function(tbl, cols, label) {
  missing <- setdiff(cols, names(tbl))
  if (length(missing) > 0)
    stop(label, " is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
}

as_logical_flag <- function(x) {
  if (is.logical(x)) return(replace_na(x, FALSE))
  if (is.numeric(x)) return(replace_na(x != 0, FALSE))
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

complete_category_colors <- function(categories) {
  categories <- sort(unique(as.character(categories)))
  missing <- setdiff(categories, names(GROUP_COLORS))
  fallback <- if (length(missing) > 0)
    setNames(scales::hue_pal(l = 55, c = 90)(length(missing)), missing) else character()
  c(GROUP_COLORS, fallback)
}

blend_with_white <- function(hex, amount = 0.65) {
  x <- grDevices::col2rgb(hex) / 255
  y <- x * (1 - amount) + amount
  grDevices::rgb(y[1, ], y[2, ], y[3, ])
}

make_fill <- function(category, intensity, colors) {
  base  <- unname(colors[as.character(category)])
  light <- blend_with_white(base, 0.72)
  vapply(seq_along(base), function(i) {
    palette <- grDevices::colorRampPalette(c(light[[i]], base[[i]]))(100)
    palette[pmax(1, pmin(100, round(intensity[[i]] * 99) + 1))]
  }, character(1))
}

reorder_within <- function(x, by, within, fun = mean, sep = "___")
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun)

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  scale_y_discrete(..., labels = function(x) gsub(reg, "", x))
}

# Build the per-trait count + union-coloc-rate table. `rate_type` selects which
# union rate column drives pct_coloc_union. Error messages name the argument
# plainly (the original used optparse-style "--gene" labels).
build_plot_data <- function(gene, coloc_rate, min_coloc_genes = 1L,
                            rate_type = "stringent") {
  rate_type <- tolower(trimws(rate_type))
  if (length(rate_type) != 1 || is.na(rate_type) || !rate_type %in% names(RATE_COLUMNS))
    stop("rate_type must be one of: ", paste(names(RATE_COLUMNS), collapse = ", "),
         call. = FALSE)
  rate_column <- unname(RATE_COLUMNS[[rate_type]])
  require_columns(gene, c("gene", "trait", "trait_category", "any_coloc"), "gene table")
  require_columns(coloc_rate, c("trait", "layer", unname(RATE_COLUMNS)), "coloc_rate table")
  if (length(min_coloc_genes) != 1 || is.na(min_coloc_genes) ||
      min_coloc_genes < 1 || min_coloc_genes != as.integer(min_coloc_genes))
    stop("min_coloc_genes must be a positive integer.", call. = FALSE)

  gene$any_coloc <- as_logical_flag(gene$any_coloc)
  coloc_counts <- gene %>%
    filter(any_coloc, !is.na(gene), gene != "", !is.na(trait_category), trait_category != "") %>%
    distinct(gene, trait, trait_category) %>%
    count(trait, trait_category, name = "count")

  pct_coloc <- coloc_rate %>%
    filter(layer == "union") %>%
    transmute(trait, pct_stringent_union, pct_lenient_union,
              pct_coloc_union = .data[[rate_column]])
  dup <- pct_coloc %>% count(trait) %>% filter(n > 1) %>% pull(trait)
  if (length(dup) > 0)
    stop("coloc_rate has multiple union rows for trait(s): ",
         paste(dup, collapse = ", "), call. = FALSE)

  combined <- coloc_counts %>% left_join(pct_coloc, by = "trait")
  missing_rate <- combined %>% filter(is.na(pct_coloc_union)) %>% pull(trait)
  if (length(missing_rate) > 0)
    message("Dropped (no union ", rate_type, " coloc rate): ",
            paste(unique(missing_rate), collapse = ", "))
  combined <- combined %>%
    filter(!is.na(pct_coloc_union)) %>% mutate(rate_type = .env$rate_type)
  trimmed <- combined %>% filter(count < min_coloc_genes) %>% pull(trait)
  if (length(trimmed) > 0)
    message("Trimmed ", length(trimmed), " trait(s) with fewer than ",
            min_coloc_genes, " colocalizing genes.")
  combined %>% filter(count >= min_coloc_genes)
}

build_empty_plot <- function() {
  ggplot() +
    annotate("text", x = 0, y = 0,
             label = "No traits have both colocalizing genes and the selected union coloc rate",
             color = "grey35", size = 4) +
    xlim(-1, 1) + ylim(-1, 1) + theme_void()
}

build_plot <- function(combined) {
  colors <- complete_category_colors(combined$trait_category)
  plotdat <- combined %>%
    group_by(trait_category) %>%
    mutate(
      intensity = if (max(count) == min(count)) rep(0.75, n())
                  else (count - min(count)) / (max(count) - min(count)),
      intensity = 0.25 + 0.75 * intensity,
      fill_col  = make_fill(trait_category, intensity, colors)) %>%
    ungroup() %>%
    mutate(
      trait_category = forcats::fct_drop(forcats::fct_reorder(trait_category, count, .fun = sum)),
      trait          = reorder_within(trait, count, trait_category),
      base_col       = unname(colors[as.character(trait_category)]))

  long <- bind_rows(
    plotdat %>% transmute(trait_category, trait, base_col, fill_col, metric = "a_count", value = count),
    plotdat %>% transmute(trait_category, trait, base_col, fill_col, metric = "b_rate",  value = pct_coloc_union))
  xline <- -max(plotdat$count) * 0.03
  seg  <- plotdat %>% count(trait_category, name = "nbar") %>% mutate(metric = "a_count")
  vref <- tibble(metric = "b_rate", xi = median(plotdat$pct_coloc_union, na.rm = TRUE))
  metric_labels <- c(a_count = "Colocalizing Genes", b_rate = "CS coloc rate ")
  category_labeller <- c(GROUP_LABELS, setNames(
    setdiff(as.character(unique(plotdat$trait_category)), names(GROUP_LABELS)),
    setdiff(as.character(unique(plotdat$trait_category)), names(GROUP_LABELS))))

  ggplot(long, aes(y = trait)) +
    geom_col(data = ~ subset(.x, metric == "a_count"),
             aes(x = value, fill = fill_col), width = 0.85, show.legend = FALSE) +
    geom_segment(data = seg, inherit.aes = FALSE,
                 aes(y = 0.5, yend = nbar + 0.5, x = xline, xend = xline),
                 color = "grey60", linewidth = 0.5) +
    geom_vline(data = vref, aes(xintercept = xi),
               linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_segment(data = ~ subset(.x, metric == "b_rate"),
                 aes(x = 0, xend = value, yend = trait, color = base_col), linewidth = 0.4) +
    geom_point(data = ~ subset(.x, metric == "b_rate"),
               aes(x = value, color = base_col), size = 1.6) +
    scale_fill_identity() + scale_color_identity() + scale_y_reordered() +
    facet_grid(trait_category ~ metric, scales = "free", space = "free_y", switch = "y",
               labeller = labeller(
                 trait_category = as_labeller(category_labeller, default = label_value),
                 metric = as_labeller(metric_labels))) +
    ggh4x::facetted_pos_scales(x = list(
      scale_x_continuous(),
      scale_x_continuous(labels = function(x) paste0(x, "%"),
                         expand = expansion(mult = c(0, 0.08))))) +
    ggh4x::force_panelsizes(cols = c(3, 1.3)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.spacing.x = grid::unit(8, "pt"), panel.spacing.y = grid::unit(2, "pt"),
          strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 12),
          strip.text.x = element_text(size = 12))
}

# Builder: the per-trait coloc table (lenient union rate, matching the figure).
coloc_by_trait <- function(refresh = FALSE)
  build_or_cache("coloc_by_trait", refresh = refresh, {
    build_plot_data(read_any(input_path("coloc_gene")),
                    read_any(input_path("coloc_rate")),
                    min_coloc_genes = 1, rate_type = "lenient")
  })
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `Rscript tests/test_fig4_helpers.R`
Expected: PASS — 9 `ok - ...` lines then `test_fig4_helpers.R (Task 2): all assertions passed`

- [ ] **Step 6: Verify the builder runs against the fixture**

Run:
```bash
Rscript -e '
  source("tests/fixtures.R"); fx <- make_fig4_fixture(); use_fig4_fixture(fx)
  source("fig4/00_setup.R")
  cb <- coloc_by_trait()
  stopifnot(nrow(cb) == 2L, all(cb$count == 8L))
  cat("coloc_by_trait OK\n")'
```
Expected: `coloc_by_trait OK`

- [ ] **Step 7: Commit**

```bash
git add Fig4/00_setup.R tests/fixtures.R tests/test_fig4_helpers.R
git commit -m "Add Fig4 coloc-by-trait helpers and builder

Folds Panel A's plotting/data helpers into 00_setup.R and adds the coloc_by_trait
builder. build_plot_data's error messages are reworded from the original optparse
'--gene' style to plain argument names. tests/fixtures.R gains a mock Fig4 project."
```

---

## Task 3: `read_conditional_twas` and the `gene_sets` builder (Panel D upstream)

**Files:**
- Modify: `Fig4/00_setup.R` (append)
- Modify: `tests/test_fig4_helpers.R` (append)

**Interfaces:**
- Consumes: `build_or_cache`, `input_path`, `input_files`, `read_any`.
- Produces: `read_conditional_twas(refresh = FALSE)`; `derive_gene_sets(clingen_df, drug_df) -> list(druggable, approved, mendelian, heme_immune, clingen_cleaned)`; `gene_sets(refresh = FALSE)` (same list, read from the manifest).

- [ ] **Step 1: Write the failing test**

Append to `tests/test_fig4_helpers.R`, before the final `cat(...)`:

```r
# --- derive_gene_sets ------------------------------------------------------
clingen_mock <- tibble::tibble(
  targetId          = c("ENSG1", "ENSG1", "ENSG2", "ENSG3"),
  targetFromSourceId = c("S1", "S1", "S2", "S3"),
  studyId           = c("st1", "st1", "st2", "st3"),
  confidence        = c("Definitive", "Definitive", "Strong", "Limited"),
  allelicRequirements = list("AR", "AR", "AD", "AD"),
  diseaseFromSource = c("hereditary anemia", "hereditary anemia",
                        "dilated cardiomyopathy", "some limited disease"))
drug_mock <- tibble::tibble(
  gene_name = c("geneA", "geneB", "geneB", ""),
  approved  = c(TRUE, FALSE, FALSE, TRUE))

gs <- derive_gene_sets(clingen_mock, drug_mock)
ok("mendelian keeps only Definitive/Strong genes",
   setequal(gs$mendelian, c("ENSG1", "ENSG2")))
ok("Limited-confidence genes are excluded", !"ENSG3" %in% gs$mendelian)
ok("heme-immune picks the anemia gene", "ENSG1" %in% gs$heme_immune)
ok("heme-immune excludes the cardiomyopathy gene", !"ENSG2" %in% gs$heme_immune)
ok("druggable upper-cases symbols and drops blanks",
   setequal(gs$druggable, c("GENEA", "GENEB")))
ok("approved keeps only approved==TRUE non-blank symbols",
   setequal(gs$approved, "GENEA"))
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript tests/test_fig4_helpers.R`
Expected: FAIL — `could not find function "derive_gene_sets"`

- [ ] **Step 3: Write the implementation**

Append to `Fig4/00_setup.R`:

```r
# =============================================================================
# Panel C/D upstream — conditional TWAS read, and gene-set annotations
# =============================================================================
read_conditional_twas <- function(refresh = FALSE)
  build_or_cache("conditional_twas", refresh = refresh, {
    read_any(input_path("conditional_twas"))
  })

# Heme / immune disease keyword pattern (hand-curated; see NOTES.md).
HEME_IMMUNE_KW <- str_c(c(
  "an[a]?emia", "thrombocyt", "neutropen", "leuk[a]?emia", "lymphoma",
  "h[a]?emoglobin", "platelet", "h[a]?emolytic", "spherocytosis", "sideroblastic",
  "myelodysplas", "granulocyt", "erythrocyt", "reticulocyt", "coagulation",
  "thromb", "bleeding", "h[a]?emophilia", "neutrophil", "eosinophil", "basophil",
  "monocyt", "myeloid", "h[a]?ematopoie", "pancytopenia", "h[a]?emochromatosis",
  "polycythemia", "thalass", "macrocyt", "blood cell", "bone marrow", "hemostasis",
  "porphyria", "fanconi", "immunodefic", "immune", "autoimmun", "autoinflammat",
  "immunodysreg", "agammaglobulin", "lymphoproliferat", "complement",
  "immunoglobulin", "hypogammaglobulin", "granulomatous", "histiocyt",
  "hemophagocytic", "interferonopathy"), collapse = "|")

# Pure gene-set derivation (testable without arrow / the parquet read).
derive_gene_sets <- function(clingen_df, drug_df) {
  clingen_cleaned <- clingen_df %>%
    distinct(targetId, targetFromSourceId, studyId, confidence,
             allelicRequirements, diseaseFromSource) %>%
    filter(confidence == "Definitive" | confidence == "Strong") %>%
    unnest(allelicRequirements)
  heme_immune_genes <- clingen_cleaned %>%
    filter(str_detect(tolower(diseaseFromSource), HEME_IMMUNE_KW)) %>%
    pull(targetId) %>% unique()

  dg <- as.data.table(drug_df)
  dg[, sym := toupper(gene_name)]
  druggable <- unique(dg[sym != "", sym])
  approved  <- unique(dg[approved == TRUE & sym != "", sym])

  list(druggable = druggable, approved = approved,
       mendelian = unique(clingen_cleaned$targetId),
       heme_immune = heme_immune_genes, clingen_cleaned = clingen_cleaned)
}

gene_sets <- function(refresh = FALSE)
  build_or_cache("gene_sets", refresh = refresh, {
    clingen_dir <- input_path("clingen")
    pq <- list.files(clingen_dir, pattern = "\\.parquet$", full.names = TRUE)
    if (length(pq) == 0) stop("No .parquet files under ", clingen_dir)
    clingen_df <- arrow::read_parquet(pq[1]) %>%
      select(targetFromSourceId, allelicRequirements, confidence, score,
             studyId, targetId, diseaseFromSource)
    drug_df <- read_any(input_path("drug_interactions"))
    derive_gene_sets(clingen_df, drug_df)
  })
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/test_fig4_helpers.R`
Expected: PASS — 15 `ok - ...` lines then the Task 2 completion line

- [ ] **Step 5: Commit**

```bash
git add Fig4/00_setup.R tests/test_fig4_helpers.R
git commit -m "Add Fig4 conditional-TWAS read and gene-set builders

read_conditional_twas caches the raw table for Panels C and D. gene-set
derivation (druggable/approved/Mendelian/heme-immune) is factored into a pure
derive_gene_sets() tested on in-memory frames; the gene_sets() builder wraps it
around the ClinGen parquet + drug-interaction reads."
```

---

## Task 4: QQ helpers and the `qq_match_source` builder (Panel C upstream)

**Files:**
- Modify: `Fig4/00_setup.R` (append)
- Modify: `tests/test_fig4_helpers.R` (append)

**Interfaces:**
- Consumes: `build_or_cache`, `read_conditional_twas`.
- Produces: `prepare_qq_match_source(data)`, `calculate_qq_data(matched_data, traits, common_cap = 40, rare_cap = 12.5)`, `match_colocalized_pairs(data, seed = 20260713, ratio = 1, caliper = 0.2, replace = FALSE) -> list(fit, data)`, `plot_qq(qq_data)`, `qq_match_source(refresh = FALSE)`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_fig4_helpers.R`, before the final `cat(...)`:

```r
# --- calculate_qq_data (pure; no MatchIt needed) ---------------------------
matched_mock <- tibble::tibble(
  group = rep(c("Colocalized genes", "Matched non-colocalized genes"), each = 4),
  common_pvalue = c(1e-8, 1e-4, 1e-2, 0.5, 0.9, 0.5, 0.2, 0.05),
  p_rare_cond   = c(1e-20, 1e-3, 0.02, 0.7, 0.8, 0.6, 0.3, 0.1))
qq <- calculate_qq_data(matched_mock, traits = character(0)) # no trait filter col
```

Wait — `calculate_qq_data` filters on `coloc_trait`; add that column:

```r
matched_mock$coloc_trait <- "eosinophil"
qq <- calculate_qq_data(matched_mock, traits = "eosinophil")
ok("calculate_qq_data returns both tests", setequal(unique(qq$test), c("common", "rare")))
ok("expected/observed columns present",
   all(c("expected", "observed", "observed_display", "capped") %in% names(qq)))
ok("common p-values are capped at 40",  all(qq$observed_display[qq$test == "common"] <= 40))
ok("rare p-values are capped at 12.5",  all(qq$observed_display[qq$test == "rare"]   <= 12.5))
ok("the 1e-20 rare hit is flagged capped",
   any(qq$capped[qq$test == "rare"]))

# --- prepare_qq_match_source (pure filtering/slicing) ----------------------
cond_mock <- tibble::tibble(
  status = "ok",
  gene = c("g1", "g1", "g2"),
  coloc_trait = "eosinophil",
  molecular_trait_id = c("m1", "m2", "m3"),
  any_coloc = c(TRUE, TRUE, FALSE),
  z_rare_cond = c(2, 2, 1), rare_zscore = c(2, 2, 1), Vcor = c(0.1, 0.1, 0.1),
  GLCP = c(0.9, 0.8, NA), p_rare_cond = c(0.01, 0.2, 0.5),
  mean_TPM = c(10, 10, 5), shet = c(0.01, 0.01, 0.02), n_rare = c(3, 3, 1))
src <- prepare_qq_match_source(cond_mock)
ok("prepare_qq_match_source keeps one row per (gene, coloc_trait)",
   nrow(src) == 2L)
ok("it adds the modelling covariates",
   all(c("log_tpm", "log_shet", "log_n_rare", "abs_vcor", "coloc") %in% names(src)))

# --- match_colocalized_pairs (only if MatchIt is installed) ----------------
if (requireNamespace("MatchIt", quietly = TRUE)) {
  big <- dplyr::bind_rows(lapply(1:40, function(i) tibble::tibble(
    status = "ok", gene = paste0("g", i), coloc_trait = "eosinophil",
    molecular_trait_id = paste0("m", i),
    any_coloc = i %% 2 == 0,
    z_rare_cond = rnorm(1, 2), rare_zscore = rnorm(1, 2), Vcor = 0.1,
    GLCP = 0.9, p_rare_cond = runif(1),
    mean_TPM = runif(1, 1, 50), shet = runif(1, 0.001, 0.1), n_rare = sample(1:10, 1))))
  set.seed(1)
  mr <- match_colocalized_pairs(prepare_qq_match_source(big))
  ok("match returns a fit and a matched data frame",
     all(c("fit", "data") %in% names(mr)) && nrow(mr$data) > 0)
  ok("matched data has the two comparison groups",
     setequal(levels(mr$data$group),
              c("Matched non-colocalized genes", "Colocalized genes")))
} else {
  cat("  (skipped match_colocalized_pairs: MatchIt not installed)\n")
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript tests/test_fig4_helpers.R`
Expected: FAIL — `could not find function "calculate_qq_data"`

- [ ] **Step 3: Write the implementation**

Append to `Fig4/00_setup.R`:

```r
# =============================================================================
# Panel C helpers — conditional-TWAS QQ with propensity matching
# =============================================================================
prepare_qq_match_source <- function(data) {
  data %>%
    filter(status == "ok") %>%
    distinct(gene, coloc_trait, molecular_trait_id, .keep_all = TRUE) %>%
    mutate(any_coloc = coalesce(any_coloc, FALSE),
           inflation = abs(z_rare_cond) / pmax(abs(rare_zscore), 0.5),
           artifact  = coalesce(inflation > 3 & abs(Vcor) > 0.3, FALSE)) %>%
    filter(!artifact) %>%
    group_by(gene, coloc_trait) %>%
    arrange(desc(coalesce(GLCP, -Inf)), p_rare_cond, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(coloc = as.integer(any_coloc), log_tpm = log1p(mean_TPM),
           log_shet = log10(shet), log_n_rare = log1p(n_rare), abs_vcor = abs(Vcor)) %>%
    filter(is.finite(log_tpm), is.finite(log_shet),
           is.finite(log_n_rare), is.finite(abs_vcor))
}

calculate_qq_data <- function(matched_data, traits, common_cap = 40, rare_cap = 12.5) {
  matched_data %>%
    filter(as.character(coloc_trait) %in% traits) %>%
    select(group, common = common_pvalue, rare = p_rare_cond) %>%
    pivot_longer(c(common, rare), names_to = "test", values_to = "p_value") %>%
    filter(is.finite(p_value), between(p_value, 0, 1)) %>%
    group_by(test, group) %>%
    arrange(p_value, .by_group = TRUE) %>%
    mutate(rank = row_number(), n = n(),
           expected = -log10(rank / (n + 1)),
           observed = -log10(pmax(p_value, .Machine$double.xmin)),
           y_cap = case_when(test == "common" ~ common_cap,
                             test == "rare" ~ rare_cap, .default = Inf),
           observed_display = pmin(observed, y_cap),
           capped = observed > y_cap) %>%
    ungroup()
}

match_colocalized_pairs <- function(data, seed = 20260713, ratio = 1,
                                    caliper = 0.2, replace = FALSE) {
  set.seed(seed)
  match_fit <- MatchIt::matchit(
    coloc ~ log_tpm + log_shet + log_n_rare + abs_vcor, data = data,
    method = "nearest", distance = "glm", exact = ~ coloc_trait,
    ratio = ratio, replace = replace, caliper = caliper, std.caliper = TRUE)
  matched_data <- MatchIt::match.data(match_fit, drop.unmatched = TRUE) %>%
    mutate(group = if_else(coloc == 1L, "Colocalized genes", "Matched non-colocalized genes"),
           group = factor(group, levels = c("Matched non-colocalized genes", "Colocalized genes")))
  list(fit = match_fit, data = matched_data)
}

plot_qq <- function(qq_data) {
  ggplot(qq_data, aes(x = expected, y = observed_display, colour = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.55, colour = "#666666") +
    geom_point(size = 1, alpha = 0.8) +
    facet_wrap(~test, scales = "free_y",
               labeller = as_labeller(c(common = "Common", rare = "Rare"))) +
    scale_colour_manual(values = c("Matched non-colocalized genes" = "#2B2B2B",
                                   "Colocalized genes" = "#275D8E"), name = NULL) +
    guides(colour = guide_legend(override.aes = list(shape = 16, size = 2.5, alpha = 1))) +
    labs(x = expression(Expected ~ -log[10](P)), y = expression(Observed ~ -log[10](P))) +
    theme_classic(base_size = 12) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 12),
          panel.spacing = grid::unit(2, "lines"), legend.position = "top",
          legend.text = element_text(size = 10), legend.key.width = grid::unit(1.2, "lines"))
}

qq_match_source <- function(refresh = FALSE)
  build_or_cache("qq_match_source", refresh = refresh, {
    match_colocalized_pairs(prepare_qq_match_source(read_conditional_twas()))$data
  })
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/test_fig4_helpers.R`
Expected: PASS — the `calculate_qq_data` and `prepare_qq_match_source` assertions pass; the MatchIt block either passes or prints the skip line.

- [ ] **Step 5: Commit**

```bash
git add Fig4/00_setup.R tests/test_fig4_helpers.R
git commit -m "Add Fig4 conditional-TWAS QQ helpers and matching builder

prepare_qq_match_source / calculate_qq_data / plot_qq and the MatchIt
propensity-matching wrapper, plus the qq_match_source builder. calculate_qq_data
and prepare_qq_match_source are tested directly; the MatchIt step is exercised
when the package is installed and skipped otherwise."
```

---

## Task 5: The `gene_twas_data` builder and `fit_enrich` (Panel D upstream)

**Files:**
- Modify: `Fig4/00_setup.R` (append)
- Modify: `tests/test_fig4_helpers.R` (append)

**Interfaces:**
- Consumes: `build_or_cache`, `read_conditional_twas`, `gene_sets`.
- Produces: `gene_twas_data(refresh = FALSE)` (per-gene summary with `max_rare`, `max_common`, `log_shet`, `log_tpm` = `log10(TPM+1)`, `ever_rare`, `ever_common`, `druggable`, `approved`, `MendelianGene`, `HemeImmuneGene`, `TPM`); `fit_enrich(outcome, predictors, data)`; `ENRICH_OUTCOMES`, `ENRICH_PREDICTORS`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_fig4_helpers.R`, before the final `cat(...)`:

```r
# --- fit_enrich ------------------------------------------------------------
enr_mock <- tibble::tibble(
  MendelianGene = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
  max_rare = c(3, 0.2, 2.5, 0.1, 2.8, 0.3, 3.1, 0.05),
  max_common = c(1, 1, 1, 1, 1, 1, 1, 1),
  log_shet = rnorm(8), TPM = runif(8, 1, 50))
res <- fit_enrich("MendelianGene", c("max_rare", "max_common", "log_shet", "log2(TPM)"),
                  data = enr_mock)
ok("fit_enrich returns a tidy OR table with a max_rare term",
   "max_rare" %in% res$term && all(c("estimate", "conf.low", "conf.high") %in% names(res)))
ok("estimates are exponentiated (odds ratios, > 0)", all(res$estimate > 0))
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript tests/test_fig4_helpers.R`
Expected: FAIL — `could not find function "fit_enrich"`

- [ ] **Step 3: Write the implementation**

Append to `Fig4/00_setup.R`:

```r
# =============================================================================
# Panel D upstream — colocalized conditional-TWAS enrichment
# -----------------------------------------------------------------------------
# NOTE: the original computed ColocAdjustedConditionalTWAS TWICE. The first copy
# referenced an undefined `fdr_rare_sig` and fed GeneTWASData; the second copy
# recomputed it cleanly (BH adjustment) but went unused. We keep ONLY the clean
# version and build GeneTWASData from it (see NOTES.md).
# =============================================================================
gene_twas_data <- function(refresh = FALSE)
  build_or_cache("gene_twas_data", refresh = refresh, {
    gs  <- gene_sets()
    adj <- read_conditional_twas() %>%
      distinct(gene, coloc_trait, molecular_trait_id, .keep_all = TRUE) %>%
      mutate(inflation = abs(z_rare_cond) / pmax(abs(rare_zscore), 0.5),
             artifact_primary = coalesce(inflation > 3 & abs(Vcor) > 0.3, FALSE)) %>%
      filter(!artifact_primary, any_coloc == TRUE) %>%
      group_by(TWAS_path) %>%
      mutate(fdr_rare_cond   = p.adjust(p_rare_cond,   method = "BH"),
             fdr_common      = p.adjust(common_pvalue, method = "BH"),
             rare_zscore_fdr = p.adjust(rare_pvalue,   method = "BH"),
             fdr_full        = p.adjust(full_pvalue,   method = "BH")) %>%
      ungroup() %>%
      mutate(genuine_rare = p_rare_cond < 0.05 & rare_pvalue < 0.05,
             t_fdrcond    = fdr_rare_cond < 0.10 & rare_pvalue < 0.05,
             t_fdrboth    = fdr_rare_cond < 0.10 & rare_zscore_fdr < 0.10,
             gene = factor(gene), coloc_trait = factor(coloc_trait),
             druggable = toupper(as.character(Description)) %in% gs$druggable,
             approved  = toupper(as.character(Description)) %in% gs$approved)

    # Guard: the artifact rule must have removed all offending rows.
    stopifnot(nrow(filter(adj, abs(Vcor) > 0.3, inflation > 3)) == 0)

    adj %>%
      group_by(gene, Description) %>%
      summarise(max_rare    = max(abs(z_rare_cond),   na.rm = TRUE),
                max_common  = max(abs(common_zscore), na.rm = TRUE),
                TPM         = mean(mean_TPM, na.rm = TRUE),
                druggable   = first(druggable), approved = first(approved),
                ever_common = any(fdr_common < 0.1, na.rm = TRUE),
                ever_rare   = any(t_fdrcond),
                shet        = first(shet), .groups = "drop") %>%
      mutate(log_shet = log10(shet), log_tpm = log10(TPM + 1),
             MendelianGene  = as.character(gene) %in% gs$clingen_cleaned$targetId,
             HemeImmuneGene = as.character(gene) %in% gs$heme_immune) %>%
      filter(!is.na(log_shet))
  })

fit_enrich <- function(outcome, predictors, data)
  glm(reformulate(predictors, response = outcome), data = data, family = binomial()) %>%
    broom::tidy(conf.int = TRUE, exponentiate = TRUE)

ENRICH_OUTCOMES <- c(
  approved       = "Approved Drug Target",
  druggable      = "Druggable gene",
  MendelianGene  = "Mendelian Gene\n(Any)",
  HemeImmuneGene = "Mendelian Gene\n(Immune/Hematological)")

ENRICH_PREDICTORS <- list(
  continuous = c("max_rare", "max_common", "log_shet", "log2(TPM)"),
  binary     = c("ever_rare", "ever_common", "log_shet", "log2(TPM)"))
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/test_fig4_helpers.R`
Expected: PASS — the `fit_enrich` assertions pass; final line `test_fig4_helpers.R (Task 2): all assertions passed`

- [ ] **Step 5: Commit**

```bash
git add Fig4/00_setup.R tests/test_fig4_helpers.R
git commit -m "Add Fig4 gene-level TWAS enrichment builder

gene_twas_data keeps ONLY the clean (BH-adjusted) ColocAdjustedConditionalTWAS
and builds the per-gene GeneTWASData from it, joining the gene-set flags. This
deletes the original's first, duplicate definition, which referenced the
undefined fdr_rare_sig. fit_enrich plus the outcome/predictor constants complete
Panel D's upstream."
```

---

## Task 6: Panels A and B

**Files:**
- Create: `Fig4/panels/panelA.R`, `Fig4/panels/panelB.R`

**Interfaces:**
- Consumes: `coloc_by_trait`, `build_plot`, `input_path`, `read_any`, `OBJ_DIR`.
- Produces: `fig4/objects/A.rds`, `fig4/objects/B.rds`.

- [ ] **Step 1: Write panel A**

Create `Fig4/panels/panelA.R`:

```r
# Panel A — colocalizing-gene counts + credible-set coloc rate, by trait category
# Inputs: coloc_gene, coloc_rate
source("fig4/00_setup.R")

coloc_plot <- coloc_by_trait() %>%
  filter(trait_category != "neuropsychiatric", trait_category != "renal_urological") %>%
  filter(count > 5) %>%
  group_by(trait_category) %>%
  arrange(desc(count)) %>%
  filter(row_number() < 11) %>%
  build_plot() +
  theme(axis.text.x = element_text(size = 10))

saveRDS(coloc_plot, file.path(OBJ_DIR, "A.rds"))
message("Panel A built.")
```

- [ ] **Step 2: Write panel B**

Create `Fig4/panels/panelB.R`:

```r
# Panel B — R^2 increase (%) contributed by rare variants
# Inputs: r2_twas
source("fig4/00_setup.R")

R2_combined <- read_any(input_path("r2_twas")) %>%
  mutate(H2_upper = (SE * 1.96) + H2, adj_R2 = r.squared / H2) %>%
  filter(Pval < 1e-5) %>%
  select(phenotype_id, AF_threshold, adj_R2, H2, Pval) %>%
  pivot_wider(names_from = AF_threshold, values_from = adj_R2) %>%
  mutate(difference = (`0` - `0.01`))

R2_dist <- R2_combined %>%
  ggplot(aes(x = 100 * difference)) +
  geom_histogram(color = "white") +
  ylab("Number of genes") +
  xlim(-20, 40) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  labs(x = expression(R^2 ~ increase ~ "(%)"))

saveRDS(R2_dist, file.path(OBJ_DIR, "B.rds"))
message("Panel B built.")
```

- [ ] **Step 3: Verify no panel sets its own tag**

Run: `grep -n "labs(tag" Fig4/panels/*.R; test $? -eq 1 && echo "OK: no panel sets a tag"`
Expected: `OK: no panel sets a tag`

- [ ] **Step 4: Build Panel A against the fixture**

Run:
```bash
Rscript -e '
  source("tests/fixtures.R"); fx <- make_fig4_fixture(); use_fig4_fixture(fx)
  source("fig4/panels/panelA.R")
  o <- readRDS(file.path(Sys.getenv("FIG4_OBJ"), "A.rds"))
  stopifnot(inherits(o, "ggplot") || inherits(o, "gg"))
  cat("panel A ->", class(o)[1], "\n")'
```
Expected: `Panel A built.` then `panel A -> ...`
(The fixture has only hematological + lipid categories with count 8; Panel A's `count > 5` keeps them. Panel B needs an `r2_twas` fixture, added in Task 8's end-to-end step — building B standalone here is optional.)

- [ ] **Step 5: Commit**

```bash
git add Fig4/panels/panelA.R Fig4/panels/panelB.R
git commit -m "Add Fig4 panels A and B

Panel A (coloc counts + CS rate by trait category) via coloc_by_trait + build_plot;
Panel B (R^2 increase histogram) from the r2_twas table. Neither sets its own tag."
```

---

## Task 7: Panels C and D

**Files:**
- Create: `Fig4/panels/panelC.R`, `Fig4/panels/panelD.R`

**Interfaces:**
- Consumes: `qq_match_source`, `calculate_qq_data`, `plot_qq`, `BLOOD_TRAITS`, `gene_twas_data`, `fit_enrich`, `ENRICH_OUTCOMES`, `ENRICH_PREDICTORS`, `OBJ_DIR`.
- Produces: `fig4/objects/C.rds`, `fig4/objects/D.rds`.

- [ ] **Step 1: Write panel C**

Create `Fig4/panels/panelC.R`:

```r
# Panel C — conditional-TWAS QQ, colocalized vs matched non-colocalized (blood)
# Inputs: conditional_twas
source("fig4/00_setup.R")

qq_blood <- calculate_qq_data(qq_match_source(), traits = BLOOD_TRAITS)

p_blood_qq <- qq_blood %>%
  mutate(group = case_when(str_detect(group, "Matched") ~ "Null", TRUE ~ group)) %>%
  plot_qq() +
  theme(legend.position = "none")

saveRDS(p_blood_qq, file.path(OBJ_DIR, "C.rds"))
message("Panel C built.")
```

- [ ] **Step 2: Write panel D**

Create `Fig4/panels/panelD.R`:

```r
# Panel D — gene-set enrichment odds-ratio forest (rare vs common TWAS signal)
# Inputs: conditional_twas, clingen, drug_interactions
source("fig4/00_setup.R")

gtd <- gene_twas_data()

# Every outcome x predictor-set combination is fit; the panel plots only the
# continuous max_rare / max_common terms. The binary-predictor fits are computed
# but not plotted (see NOTES.md).
enrichment_results <- tidyr::expand_grid(outcome = names(ENRICH_OUTCOMES),
                                         model   = names(ENRICH_PREDICTORS)) %>%
  mutate(res = map2(outcome, model, ~ fit_enrich(.x, ENRICH_PREDICTORS[[.y]], gtd))) %>%
  unnest(res) %>%
  mutate(type = ENRICH_OUTCOMES[outcome]) %>%
  relocate(type, model)

plot_df <- enrichment_results %>%
  filter(term %in% c("max_rare", "max_common")) %>%
  select(type, model, term, OR = estimate, conf.low, conf.high, p.value) %>%
  mutate(term = str_remove_all(term, "max_|ever_|TRUE"),
         term = recode(term, rare = "Rare", common = "Common"),
         term = factor(term, levels = c("Rare", "Common")),
         type = factor(type, levels = c("Mendelian Gene\n(Any)",
                                        "Mendelian Gene\n(Immune/Hematological)",
                                        "Druggable gene", "Approved Drug Target")))

pal <- c(Rare = as.character(paletteer_d("futurevisions::mars"))[1], Common = "grey60")
p <- ggplot(plot_df, aes(OR, reorder(type, OR), colour = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  position = position_dodge(width = .55), linewidth = .8, fatten = 2.6) +
  scale_colour_manual(values = pal, name = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = c(0.9, 0.2), panel.spacing = unit(1.4, "lines")) +
  xlab("Odds Ratio (per unit increase in Z-score)") + ylab("")

saveRDS(p, file.path(OBJ_DIR, "D.rds"))
message("Panel D built.")
```

- [ ] **Step 3: Verify both parse and set no tag**

Run:
```bash
for f in Fig4/panels/panelC.R Fig4/panels/panelD.R; do Rscript -e "invisible(parse('$f'))"; done
grep -n "labs(tag" Fig4/panels/panelC.R Fig4/panels/panelD.R; test $? -eq 1 && echo "OK: no tags"
```
Expected: no parse error; `OK: no tags`

- [ ] **Step 4: Commit**

```bash
git add Fig4/panels/panelC.R Fig4/panels/panelD.R
git commit -m "Add Fig4 panels C and D

Panel C (blood-trait conditional-TWAS QQ, colocalized vs matched Null) via
qq_match_source + calculate_qq_data. Panel D (gene-set enrichment forest) via
gene_twas_data; both predictor sets are fit but only the continuous rare/common
terms are plotted, matching the original."
```

---

## Task 8: Assembly, runner, and the rendered-tag test

**Files:**
- Create: `Fig4/assemble.R`, `Fig4/run_all.R`, `tests/test_fig4_assembly.R`
- Modify: `tests/fixtures.R` (extend `make_fig4_fixture` with the remaining inputs)

**Interfaces:**
- Consumes: `fig4/objects/A.rds`–`D.rds`, `SZ_TAG`, `DPI`, `OUTDIR`.
- Produces: `Figure4.pdf`, `Figure4.png`.

- [ ] **Step 1: Write assemble.R**

Create `Fig4/assemble.R`:

```r
# =============================================================================
# Figure 4 — assembly. Loads the per-panel ggplot objects (fig4/objects/*.rds)
# and composes the final figure. Panels must be built first (fig4/panels/*.R),
# or just run fig4/run_all.R.
#
# Tagging happens HERE and only here. The layout is irregular and wraps the B/C
# row in wrap_elements, so tag_levels="A" traversal is unreliable (it mislettered
# Fig3's inset); tags are assigned explicitly per panel instead.
# =============================================================================
source("fig4/00_setup.R")

P <- setNames(lapply(LETTERS[1:4],
                     function(x) readRDS(file.path(OBJ_DIR, paste0(x, ".rds")))),
              LETTERS[1:4])

top_row <- ((P$B + labs(tag = "B")) | (P$C + labs(tag = "C"))) +
  plot_layout(widths = c(0.8, 1.2))
top_row <- wrap_elements(full = top_row)   # isolate B/C from D's axis-label spacing

right_side <- top_row / (P$D + labs(tag = "D")) +
  plot_layout(heights = c(1, 0.75))

Fig4 <- ((P$A + labs(tag = "A") +
            theme(plot.margin = margin(r = 55, unit = "pt"),
                  strip.text.y.left = element_text(angle = 0, hjust = 1, size = 13))) |
           right_side) +
  plot_layout(widths = c(7, 8)) &
  theme(plot.tag = element_text(size = SZ_TAG), plot.tag.position = c(0, 1))

ggsave(file.path(OUTDIR, "Figure4.pdf"), Fig4, width = 15, height = 10, dpi = DPI)
showtext_opts(dpi = 100)
ggsave(file.path(OUTDIR, "Figure4.png"), Fig4, width = 15, height = 10, dpi = 100)
showtext_opts(dpi = DPI)

message("Assembled: ", file.path(OUTDIR, "Figure4.pdf"), " and Figure4.png")
```

- [ ] **Step 2: Write run_all.R**

Create `Fig4/run_all.R`:

```r
# =============================================================================
# Figure 4 — end-to-end runner. From the project root:  Rscript fig4/run_all.R
#   1. panels  : each stages the inputs it needs (fig4/inputs.tsv) and builds the
#                derived tables it needs (memoized in fig4/data/) -> fig4/objects/
#   2. assemble: compose + save Figure4.pdf / Figure4.png
# The derived-data cache is keyed by NAME, not content: clear fig4/data/ after
# changing an input file.
# =============================================================================
for (p in LETTERS[1:4]) source(sprintf("fig4/panels/panel%s.R", p))
source("fig4/assemble.R")
message("run_all complete.")
```

- [ ] **Step 3: Extend the fixture with the remaining inputs**

In `tests/fixtures.R`, inside `make_fig4_fixture`, before the final `list(root = root)`, add the `r2_twas`, `conditional_twas`, and ClinGen/drug inputs so the whole pipeline can build. Add:

```r
  # R^2 TWAS: two AF thresholds per phenotype so the pivot has `0` and `0.01`
  r2 <- do.call(rbind, lapply(1:12, function(i) data.frame(
    phenotype_id = paste0("g", i), AF_threshold = c(0, 0.01),
    r.squared = c(0.30 + i/100, 0.20 + i/100), H2 = 0.5, SE = 0.05, Pval = 1e-8)))
  wtsv(r2, file.path(root, "R2_TWAS_combined.tsv"))

  # Conditional TWAS: enough colocalized + non-colocalized rows for MatchIt and glm
  ct <- do.call(rbind, lapply(1:60, function(i) data.frame(
    status = "ok", gene = paste0("ENSG", sprintf("%08d", i)),
    coloc_trait = "eosinophil", molecular_trait_id = paste0("m", i),
    any_coloc = TRUE, TWAS_path = "p1",
    z_rare_cond = stats::rnorm(1, 2), rare_zscore = stats::rnorm(1, 2),
    common_zscore = stats::rnorm(1, 1), Vcor = 0.1, GLCP = 0.9,
    p_rare_cond = stats::runif(1), rare_pvalue = stats::runif(1),
    common_pvalue = stats::runif(1), full_pvalue = stats::runif(1),
    mean_TPM = stats::runif(1, 1, 50), shet = stats::runif(1, 0.001, 0.1),
    n_rare = sample(1:10, 1), Description = paste0("GENE", i))))
  # data.table::fwrite handles the .gz extension by gzipping
  data.table::fwrite(ct, file.path(root, "ConditionalTWAS.cleaned.tsv.gz"), sep = "\t")

  # Drug interactions
  wtsv(data.frame(gene_name = paste0("GENE", 1:60),
                  approved = rep(c(TRUE, FALSE), 30)),
       file.path(root, "reference/interactions.tsv"))

  # ClinGen parquet (only if arrow can write it; gene_sets/Panel D need it)
  if (requireNamespace("arrow", quietly = TRUE)) {
    dir.create(file.path(root, "reference/evidence_clingen"), showWarnings = FALSE)
    clingen <- data.frame(
      targetFromSourceId = paste0("S", 1:60),
      allelicRequirements = "autosomal recessive",
      confidence = rep(c("Definitive", "Strong"), 30),
      score = 1, studyId = paste0("st", 1:60),
      targetId = paste0("ENSG", sprintf("%08d", 1:60)),
      diseaseFromSource = rep(c("hereditary anemia", "cardiomyopathy"), 30))
    arrow::write_parquet(clingen,
      file.path(root, "reference/evidence_clingen/part-0.parquet"))
  }
```

- [ ] **Step 4: Write the assembly tag test**

Create `tests/test_fig4_assembly.R`. It reuses the **proven** tag-extraction from
`tests/test_fig3_assembly.R` verbatim (render to a null device, grab the grid
display list, collect the bare single-letter `plot.tag` grobs) rather than a new
grob-walk:

```r
# Asserts the assembled Figure 4 carries exactly the four tags A, B, C, D, by
# rendering and inspecting the tag grobs (not by reading the code) -- the same
# check that caught Fig3's mislettered inset. Requires MatchIt (Panel C) and
# arrow (Panel D). Run from the project root.
source("tests/fixtures.R")
source("tests/helpers.R")
if (!requireNamespace("MatchIt", quietly = TRUE) ||
    !requireNamespace("arrow", quietly = TRUE)) {
  cat("  (skipped test_fig4_assembly: MatchIt/arrow not installed)\n")
  quit(status = 0)
}
fx <- make_fig4_fixture(); use_fig4_fixture(fx)

# Proven extractor, copied from tests/test_fig3_assembly.R: render the built
# figure to a null device, grab the grid display list, and collect the bare
# single-uppercase-letter text grobs (the plot.tag grobs). Nothing else in the
# figure draws a bare A-Z character.
extract_tag_labels <- function(patchwork_obj) {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  grid::grid.newpage(); print(patchwork_obj)
  grabbed <- grid::grid.grab()
  labels <- character(0)
  walk <- function(x) {
    if (is.null(x)) return(invisible())
    if (inherits(x, "text")) {
      lab <- paste(x$label, collapse = "")
      if (nchar(lab) == 1L && lab %in% LETTERS) labels <<- c(labels, lab)
    } else if (inherits(x, "gtable")) { for (gr in x$grobs) walk(gr) }
    else if (!is.null(x$children)) { for (ch in x$children) walk(ch) }
  }
  walk(grabbed); labels
}

for (p in LETTERS[1:4]) source(sprintf("fig4/panels/panel%s.R", p))
source("fig4/assemble.R")

tags <- extract_tag_labels(Fig4)
# The layout wraps the B/C row in wrap_elements, so draw order is not
# necessarily A,B,C,D; assert the SET is exactly the four tags, four in total,
# with none missing, extra, or duplicated.
ok("Figure 4 carries exactly 4 tags", length(tags) == 4L)
ok("Figure 4's tags are exactly A, B, C, D",
   identical(sort(unique(tags)), c("A", "B", "C", "D")))

cat("test_fig4_assembly.R: all assertions passed\n")
```

- [ ] **Step 5: Run the whole pipeline end-to-end against the fixture**

Run:
```bash
Rscript -e '
  source("tests/fixtures.R"); fx <- make_fig4_fixture(); use_fig4_fixture(fx)
  source("fig4/run_all.R")
  out <- Sys.getenv("FIG4_OUTDIR")
  for (f in c("Figure4.pdf","Figure4.png")) {
    p <- file.path(out, f); stopifnot(file.exists(p), file.size(p) > 1000)
    cat(sprintf("%-14s %8d bytes\n", f, file.size(p)))
  }
  cat("end-to-end OK\n")'
```
Expected: panel build messages, `run_all complete.`, two size lines > 1000, `end-to-end OK`.
(If MatchIt or arrow is absent, Panels C/D error — install them first: `Rscript -e 'install.packages(c("MatchIt","arrow"))'`.)

- [ ] **Step 6: Run the test suite**

Run: `Rscript tests/run_all.R`
Expected: all test files pass (or skip C/D-dependent assertions with the printed skip line).

- [ ] **Step 7: Commit**

```bash
git add Fig4/assemble.R Fig4/run_all.R tests/test_fig4_assembly.R tests/fixtures.R
git commit -m "Add Fig4 assembly, runner and rendered-tag test

assemble.R reproduces the original nested layout and assigns A/B/C/D explicitly
(tag_levels traversal is unreliable across the wrap_elements boundary, as it was
for Fig3). test_fig4_assembly.R asserts the rendered tag set. The fixture gains
r2_twas / conditional_twas / drug / ClinGen-parquet inputs for the end-to-end run."
```

---

## Task 9: Documentation, and retire Fig4.R

**Files:**
- Create: `Fig4/README.md`, `Fig4/NOTES.md`
- Remove: `Fig4/Fig4.R`

- [ ] **Step 1: Write NOTES.md**

Create `Fig4/NOTES.md`:

```markdown
# Figure 4 — notes and deferred work

## Result-affecting change: the duplicated conditional-TWAS table

`Fig4.R` computed `ColocAdjustedConditionalTWAS` **twice**. The first copy
(`Fig4.R:590`) referenced `fdr_rare_sig`, which is never defined
(`t_fdrboth = fdr_rare_cond < .1 & fdr_rare_sig`) — a typo for `rare_zscore_fdr`
— and used `p.adjust(method = "fdr")`. `GeneTWASData`, which Panel D depends on,
was built from that first copy. The second copy (`Fig4.R:636`) recomputed the
table cleanly with BH adjustment and a correct
`t_fdrboth = fdr_rare_cond < 0.10 & rare_zscore_fdr < 0.10`, but nothing used it.

`gene_twas_data()` keeps **only the clean (second) version** and builds
`GeneTWASData` from it. Panel D may therefore differ from a figure produced by
running the original top to bottom.

## GeneBayes was read but never used

`Fig4.R:91` read a GeneBayes s_het table from
`41588_2024_1820_MOESM4_ESM.xlsx` (sheet "Supplementary Table 1"), but nothing
downstream used it — every `shet` / `log_shet` comes from the `shet` column
already in `ConditionalTWAS`. The read is dropped, which also removes the
`readxl` dependency from Fig4. (This is the same supplement Fig1 (deferred) and
Fig3 reference.)

## Binary-predictor enrichment is computed but not plotted

Panel D fits enrichment for two predictor sets — `continuous`
(`max_rare`, `max_common`, `log_shet`, `log2(TPM)`) and `binary`
(`ever_rare`, `ever_common`, …) — but plots only the continuous `max_rare` /
`max_common` terms. Both are computed (`ENRICH_PREDICTORS` in `00_setup.R`); the
binary fits are available for inspection.

## Hand-curated heme/immune keyword list

The heme/immune gene set is defined by a regex of disease-name keywords
(`HEME_IMMUNE_KW` in `00_setup.R`) matched against ClinGen `diseaseFromSource`.
It is a curated list, not an ontology mapping.

## Inputs still missing a `gs://` origin

`inputs.tsv` records these local-only. The two `high_level_coloc_outputs/`
tables are extracted from `gs://fc-secure-…/high_level_coloc_outputs.tar.gz`
(README documents the staging command); `R2_TWAS_combined.tsv`,
`ConditionalTWAS.cleaned.tsv.gz`, `reference/interactions.tsv` and the ClinGen
parquet need URIs or pre-staging.
```

- [ ] **Step 2: Write README.md**

Create `Fig4/README.md` modelled on `Fig3/README.md`: the four-panel description
table (A coloc-by-trait, B R² increase, C conditional-TWAS QQ, D enrichment
forest), the manifest/staging section (including the tar.gz `gsutil cp … && tar
-xf` step and the ClinGen parquet), the derived-data builder table and the
name-keyed cache caveat, run/rebuild commands (`Rscript fig4/run_all.R`), the
requirements (add `MatchIt`, `arrow`, `broom`, `ggh4x`; note `readxl` is *not*
needed), the tests section, and a Data & inputs section pointing at `NOTES.md`
for the result-affecting change. Provide the exact content:

```markdown
# Figure 4 — QTL–GWAS colocalization and rare-variant TWAS

Code to generate **Figure 4** for the manuscript. The figure connects molecular
QTL colocalization with GWAS traits to downstream rare-variant signal: how many
genes colocalize per trait category and at what credible-set rate, the R² that
rare variants add to TWAS, whether colocalized genes carry inflated conditional
TWAS signal versus matched controls, and whether that signal enriches in
clinically important gene sets.

> **Data note.** The inputs derive from *All of Us* **controlled-access** data
> and are **not** in this repository. Run this on the Researcher Workbench.

## The figure

| Panel | Content | Source |
|-------|---------|--------|
| **A** | Colocalizing-gene counts + credible-set coloc rate, by trait category | `harmonized_coloc.gene.tsv.gz` + `coloc_rate_by_trait.tsv` |
| **B** | Histogram of the R² increase (%) contributed by rare variants | `R2_TWAS_combined.tsv` |
| **C** | Common/rare conditional-TWAS QQ, colocalized vs. matched non-colocalized (blood traits) | `ConditionalTWAS.cleaned.tsv.gz` |
| **D** | Odds-ratio forest: rare/common TWAS-signal enrichment across gene sets | conditional TWAS + ClinGen + drug interactions |

## Repository layout, running, tests

Run from the **project root**:

    Rscript fig4/run_all.R          # panels -> assemble -> Figure4.pdf / .png
    Rscript fig4/panels/panelD.R    # rebuild one panel
    Rscript fig4/assemble.R         # recompose
    rm -rf fig4/data && Rscript fig4/run_all.R   # force-rebuild the derived cache
    Rscript tests/run_all.R         # unit tests (no controlled-access data)

The panels share upstream work through memoized builders in `00_setup.R`
(`read_conditional_twas`, `gene_sets`, `coloc_by_trait`, `qq_match_source`,
`gene_twas_data`), cached to `fig4/data/` **by name** — clear that directory
after changing an input.

## Requirements

R ≥ 4.3, and CRAN packages `tidyverse`, `data.table`, `patchwork`, `ggh4x`,
`MatchIt`, `arrow`, `broom`, `scales`, `showtext`, `sysfonts`, `paletteer`.
(`readxl` is **not** needed — the GeneBayes read was unused and removed.)

## Data & inputs

Every input is declared in `fig4/inputs.tsv`. The two `high_level_coloc_outputs/`
tables come from a single archive:

    gsutil cp gs://fc-secure-…/high_level_coloc_outputs.tar.gz .
    tar -xf high_level_coloc_outputs.tar.gz

The ClinGen `reference/evidence_clingen/` parquet, `reference/interactions.tsv`,
`R2_TWAS_combined.tsv` and `ConditionalTWAS.cleaned.tsv.gz` are local-only and
must be pre-staged. See [`NOTES.md`](NOTES.md) — in particular the
result-affecting fix to the duplicated conditional-TWAS table (Panel D).

## Provenance

This directory replaces a single 774-line notebook export (`Fig4.R`, recoverable
from git history). The restructuring and its fixes are recorded in
[`docs/superpowers/specs/2026-07-24-fig4-cleanup-design.md`](../docs/superpowers/specs/2026-07-24-fig4-cleanup-design.md).
```

- [ ] **Step 3: Delete the superseded script**

Run:
```bash
git rm Fig4/Fig4.R
test ! -e Fig4/Fig4.R && echo "Fig4.R removed"
```
Expected: `Fig4.R removed`

- [ ] **Step 4: Full suite + parse check + end-to-end**

Run:
```bash
Rscript tests/run_all.R
for f in Fig4/00_setup.R Fig4/panels/*.R Fig4/assemble.R Fig4/run_all.R; do
  Rscript -e "invisible(parse('$f'))" || echo "PARSE FAIL: $f"
done
Rscript -e '
  source("tests/fixtures.R"); fx <- make_fig4_fixture(); use_fig4_fixture(fx)
  source("fig4/run_all.R"); cat("final end-to-end OK\n")'
```
Expected: `All test files passed.`, no `PARSE FAIL`, `final end-to-end OK`.

- [ ] **Step 5: Confirm nothing sensitive is staged**

Run:
```bash
git status --short
git check-ignore Fig4/objects/A.rds Fig4/data/gene_sets.rds && echo "caches ignored"
```
Expected: only the intended docs; `caches ignored`.

- [ ] **Step 6: Commit**

```bash
git add Fig4/README.md Fig4/NOTES.md
git add -u Fig4
git commit -m "Document Fig4 and retire the notebook export

README covers the four panels, the manifest + tar.gz/parquet staging, the two
caches, and the requirements (readxl no longer needed). NOTES records the
result-affecting conditional-TWAS dedup, the dropped GeneBayes read, the unplotted
binary-predictor enrichment and the heme/immune keyword list. Deletes Fig4/Fig4.R."
```

---

## Final verification

- [ ] **Push and open the PR**

```bash
git push -u origin fig4-cleanup
gh pr create --base main --head fig4-cleanup \
  --title "Clean up Figure 4: modular QTL-GWAS colocalization pipeline" \
  --body "$(cat <<'BODY'
Replaces the 774-line notebook export `Fig4/Fig4.R` with a Fig3-style modular
pipeline reusing `common/manifest.R`.

Spec: `docs/superpowers/specs/2026-07-24-fig4-cleanup-design.md`

## Results can change — see `Fig4/NOTES.md`

- `ColocAdjustedConditionalTWAS` was computed twice; the first copy referenced the
  undefined `fdr_rare_sig` and fed Panel D. The cleanup keeps only the clean
  (BH-adjusted) version and rebuilds `GeneTWASData` from it.
- `GeneBayes` was read but never used (shet comes from `ConditionalTWAS`); the read
  and the `readxl` dependency are removed.

## Structure

- `Fig4/00_setup.R` — paths, theme, palettes, folded helpers, and six memoized
  builders cached to `fig4/data/`.
- `Fig4/panels/panelA–D.R`, `assemble.R` (deterministic A/B/C/D tags, verified
  against rendered grobs), `run_all.R`.
- Tests: pure-transform unit tests + a rendered-tag assembly test; no
  controlled-access data required.

## Not verified here

Inputs are controlled-access; validation is unit tests + a full end-to-end build
against mock fixtures. The real run is `Rscript fig4/run_all.R` on the Workbench.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
BODY
)"
```

---

## Traceability — spec section → task

| Spec § | Requirement | Task |
|---|---|---|
| §3 | Target structure | 1, 6–9 |
| §4 | Folded helpers, library changes, DPI/theme | 1–5 |
| §5 | `build_or_cache` + six builders | 1–5 |
| §6 | Panel decomposition + layout + tagging | 6–8 |
| §7.1.1 | Dedup `ColocAdjustedConditionalTWAS`, fix `fdr_rare_sig` | 5 |
| §7.1.2 | Duplicate `genuine_rare`/`t_nominal` | 5 |
| §7.2.3 | Unused reads (`gene_dat`, `gene_data`, `blank_plot`) | not carried over (never written) |
| §7.2.3a | Drop dead `GeneBayes` read + `readxl` | 1, 3, 9 |
| §7.2.4 | Library block: +arrow +broom, −optparse −readxl | 1 |
| §7.2.5 | Reworded `build_plot_data` error messages | 2 |
| §7.2.6 | `heme_immune` column resolved | 3 |
| §7.2.7 | Remove `options(repr.plot.*)` | 1 (never carried over) |
| §7.2.8 | One `DPI` constant | 1, 8 |
| §7.3.9 | Binary predictor set documented | 7, 9 |
| §7.3.10 | GeneBayes / keyword list documented | 9 |
| §8 | Input manifest | 1 |
| §9 | Tests | 2–5, 8 |
