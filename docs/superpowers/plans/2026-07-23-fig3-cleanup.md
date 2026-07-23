# Figure 3 Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the 858-line flat notebook export `Fig3/Fig3.R` with a modular, manifest-driven, Fig2-style pipeline, fixing 21 recorded defects — including two that currently build panel B's colocalization graph from a single assay pair.

**Architecture:** The generic manifest/gsutil staging layer is extracted from `Fig1/00_setup.R` into a shared `common/manifest.R`. `Fig3/00_setup.R` sources it and adds Fig3's paths, theme, FastENLOC readers, and six *memoized derived-data builders* that cache to `fig3/data/*.rds`. Each panel script sources setup, calls only the builders it needs, and saves one ggplot object to `fig3/objects/<X>.rds`. `assemble.R` composes those objects; `run_all.R` runs everything.

**Tech Stack:** R 4.3, tidyverse, data.table, patchwork, igraph, broom, showtext/sysfonts, paletteer. Tests are plain `Rscript` files using `stopifnot`-style assertions — no test framework is added.

**Spec:** `docs/superpowers/specs/2026-07-23-fig3-cleanup-design.md`

## Global Constraints

- **Branch:** `fig3-cleanup` (already created, off `main`; spec committed as `7406a97`).
- **Working directory:** every script is run from the **project root** and refers to itself with a lowercase `fig3/...` path (macOS is case-insensitive; this matches Fig1).
- **Never commit inputs or figure outputs.** Inputs derive from *All of Us* controlled-access data.
- **`showtext` DPI must equal every `ggsave()` DPI.** Fig3 uses `DPI <- 400`.
- **`element_text()` sizes are pt; `geom_text()`/`geom_label()`/`annotate()` sizes are mm.** Never reuse a constant across the two.
- **Tagging is `assemble.R`'s job.** No panel script may call `labs(tag = ...)`.
- **Reserved manifest columns are exactly `key`, `local_subpath`, `gs_uri`.** Every other column is metadata.
- **Derived-data cache is keyed by name only**, not by input content.
- Commit after every task. Never `git add -A` (untracked staged inputs may be present).

---

## File Structure

| File | Responsibility |
|---|---|
| `common/manifest.R` | **Create.** Generic manifest reading + lazy gsutil staging. Figure-agnostic. |
| `Fig1/00_setup.R` | **Modify** (lines 74–142). Replace the inlined loader block with `source("common/manifest.R")`. |
| `Fig1/inputs.tsv` | **Modify** (header only). `assay` → `Assay`. |
| `Fig1/panels/panelH.R` | **Modify** (line 8). `annot$assay` → `annot$Assay`. |
| `Fig3/00_setup.R` | **Create.** Libraries, paths, theme, size constants, palettes, FastENLOC readers, `add_fdr_threshold_label`, `build_or_cache`, six derived-data builders. |
| `Fig3/inputs.tsv` | **Create.** 10 keys / 16 rows. |
| `Fig3/panels/panelA.R` … `panelE.R` | **Create.** One ggplot object each. |
| `Fig3/assemble.R` | **Create.** Compose → `Figure3.pdf` / `Figure3.png`. |
| `Fig3/run_all.R` | **Create.** Panels → assemble. |
| `Fig3/supp_coloc_qc.R` | **Create.** Supplemental S1. |
| `Fig3/supp_prop_gene_coloc.R` | **Create.** Supplemental S2. |
| `Fig3/README.md`, `Fig3/NOTES.md`, `Fig3/.gitignore` | **Create.** |
| `Fig3/Fig3.R` | **Delete** (last task). |
| `.gitignore` | **Modify.** Add `/enloc/`, `/Figure3.pdf`, `/Figure3.png`. |
| `tests/helpers.R` | **Create.** `ok()` assertion + `expect_error()`. |
| `tests/fixtures.R` | **Create.** Builds a mock Fig3 project (manifest + tiny input files) in a temp dir. |
| `tests/test_manifest.R` | **Create.** Task 1. |
| `tests/test_fig3_readers.R` | **Create.** Task 3. |
| `tests/test_fig3_edges.R` | **Create.** Task 6 — regression test for the two data-corrupting defects. |
| `tests/test_fig3_builders.R` | **Create.** Task 7. |
| `tests/run_all.R` | **Create.** Runs every `tests/test_*.R`. |

---

## Task 1: Extract the shared manifest loader

**Files:**
- Create: `common/manifest.R`
- Create: `tests/helpers.R`
- Create: `tests/test_manifest.R`
- Create: `tests/run_all.R`
- Modify: `Fig1/00_setup.R:74-142`
- Modify: `Fig1/inputs.tsv:1`
- Modify: `Fig1/panels/panelH.R:8`

**Interfaces:**
- Produces: `input_manifest()`, `manifest_rows(key)`, `input_files(key) -> tibble` (manifest rows + resolved `path` column), `input_path(key) -> character(1)`, `load_input(key, reader = data.table::fread, ..., attach_meta = TRUE) -> tibble`, `manifest_meta_cols(f) -> character`, and the constant `.MANIFEST_RESERVED`. All require the sourcing script to have defined the globals `MANIFEST` (manifest TSV path) and `INPUT` (staging root) first.

- [ ] **Step 1: Write the failing test helper**

Create `tests/helpers.R`:

```r
# Minimal assertion helpers shared by the test scripts.
# Each test script sources this, then calls ok() / expect_error().
ok <- function(label, condition) {
  if (!isTRUE(condition)) stop("FAIL: ", label, call. = FALSE)
  cat("  ok -", label, "\n")
}

expect_error <- function(label, expr) {
  ok(label, inherits(try(expr, silent = TRUE), "try-error"))
}
```

- [ ] **Step 2: Write the failing test**

Create `tests/test_manifest.R`:

```r
# Tests for common/manifest.R. Run from the project root:
#   Rscript tests/test_manifest.R
suppressPackageStartupMessages({
  library(data.table); library(tibble); library(purrr)
})
source("tests/helpers.R")

tmp <- file.path(tempdir(), "manifest_test")
unlink(tmp, recursive = TRUE)
dir.create(file.path(tmp, "sub"), recursive = TRUE)

wtsv <- function(df, p) write.table(df, p, sep = "\t", row.names = FALSE, quote = FALSE)
wtsv(data.frame(v = 1:2), file.path(tmp, "sub/a.tsv"))
wtsv(data.frame(v = 3:5), file.path(tmp, "sub/b.tsv"))
wtsv(data.frame(v = 9),   file.path(tmp, "sub/solo.tsv"))

writeLines(c(
  "key\ttype\tlocal_subpath\tgs_uri",
  "multi\teQTL-pQTL\tsub/a.tsv\t-",
  "multi\tsQTL-pQTL\tsub/b.tsv\t-",
  "solo\t-\tsub/solo.tsv\t-",
  "nofile\t-\tsub/missing.tsv\t-"
), file.path(tmp, "inputs.tsv"))

MANIFEST <- file.path(tmp, "inputs.tsv")
INPUT    <- tmp
source("common/manifest.R")

# A manifest column named `key` must not shadow the `key` argument. This is the
# data.table non-standard-evaluation bug found during the Fig1 work: on a
# data.table, `m[m$key == key, ]` compares the column to itself and matches
# every row.
ok("multi resolves to exactly 2 rows", nrow(input_files("multi")) == 2L)
ok("solo resolves to exactly 1 row",   nrow(input_files("solo")) == 1L)
ok("input_path returns the solo path",
   input_path("solo") == file.path(tmp, "sub/solo.tsv"))

d <- load_input("multi")
ok("load_input row-binds every file for the key", nrow(d) == 5L)
ok("metadata attaches under its own manifest column name", "type" %in% names(d))
ok("metadata values follow their own file",
   identical(d$type, c("eQTL-pQTL", "eQTL-pQTL", "sQTL-pQTL", "sQTL-pQTL", "sQTL-pQTL")))
ok("reserved columns are never attached",
   !any(c("key", "local_subpath", "gs_uri", "path") %in% names(d)))

ok("'-' metadata is read as NA and not attached",
   !"type" %in% names(load_input("solo")))

expect_error("unknown key errors",                  input_files("nope"))
expect_error("multi-file key errors in input_path", input_path("multi"))
expect_error("missing file with no gs_uri errors",  input_files("nofile"))

cat("test_manifest.R: all assertions passed\n")
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `Rscript tests/test_manifest.R`
Expected: FAIL — `cannot open file 'common/manifest.R': No such file or directory`

- [ ] **Step 4: Write the implementation**

Create `common/manifest.R`:

```r
# =============================================================================
# Shared input manifest and lazy gsutil staging. Sourced by every figure's
# 00_setup.R, which must define these two globals FIRST:
#   MANIFEST  path to the figure's TSV manifest
#   INPUT     staging root (the cwd by default)
#
# Manifest columns `key`, `local_subpath` and `gs_uri` are reserved. EVERY OTHER
# COLUMN is metadata, attached by load_input() to the rows it returns, under its
# own name. Empty or "-" is read as NA and never attached.
#
# Files are cached at INPUT/<local_subpath>. A gs_uri is copied down with gsutil
# ONLY when the local copy is missing, so a pre-staged file is never re-fetched.
# =============================================================================
.MANIFEST_RESERVED <- c("key", "local_subpath", "gs_uri")

.manifest_env <- new.env(parent = emptyenv())

input_manifest <- function() {
  # Re-read when MANIFEST changes so several figures (or tests) can share a session.
  if (is.null(.manifest_env$m) || !identical(.manifest_env$path, MANIFEST)) {
    # as_tibble() avoids data.table non-standard evaluation: the manifest has a
    # column literally named "key", which would shadow a `key` filter argument.
    .manifest_env$m <- tibble::as_tibble(
      data.table::fread(MANIFEST, sep = "\t", na.strings = c("", "-"),
                        colClasses = "character"))
    .manifest_env$path <- MANIFEST
  }
  .manifest_env$m
}

manifest_rows <- function(key) {
  m <- input_manifest()
  rows <- m[m$key == key, ]
  if (nrow(rows) == 0L) stop("Unknown input key '", key, "' in ", MANIFEST)
  rows
}

# Ensure every file for `key` exists locally (staging from gs_uri if needed) and
# return the manifest rows plus a resolved local `path` column.
input_files <- function(key) {
  rows  <- manifest_rows(key)
  paths <- file.path(INPUT, rows$local_subpath)
  for (i in seq_len(nrow(rows))) {
    lp <- paths[i]
    if (file.exists(lp)) next
    gs <- rows$gs_uri[i]
    if (is.na(gs))
      stop("Missing local input '", lp, "' and no gs_uri to stage it (key '", key, "').")
    dir.create(dirname(lp), showWarnings = FALSE, recursive = TRUE)
    message("  staging ", gs, " -> ", lp)
    status <- system2("gsutil", c("cp", gs, lp))
    if (status != 0L || !file.exists(lp)) stop("gsutil cp failed for ", gs)
  }
  rows$path <- paths
  tibble::as_tibble(rows)
}

# Single-file key -> its local path (staged if needed).
input_path <- function(key) {
  f <- input_files(key)
  if (nrow(f) != 1L)
    stop("Key '", key, "' maps to ", nrow(f), " files; use input_files()/load_input().")
  f$path
}

# Metadata columns of a resolved input_files() tibble.
manifest_meta_cols <- function(f) setdiff(names(f), c(.MANIFEST_RESERVED, "path"))

# Read every file for `key` and row-bind, attaching each metadata column under
# its own name. `reader` defaults to fread; override it (and pass extra args via
# ...) for other formats.
load_input <- function(key, reader = data.table::fread, ..., attach_meta = TRUE) {
  f    <- input_files(key)
  meta <- manifest_meta_cols(f)
  purrr::map_dfr(seq_len(nrow(f)), function(i) {
    d <- tibble::as_tibble(reader(f$path[i], ...))
    if (attach_meta) {
      for (mc in meta) {
        v <- f[[mc]][i]
        if (!is.na(v)) d[[mc]] <- v
      }
    }
    d
  })
}
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `Rscript tests/test_manifest.R`
Expected: PASS — 11 `ok - ...` lines then `test_manifest.R: all assertions passed`

- [ ] **Step 6: Migrate Fig1 to the shared loader**

Replace `Fig1/00_setup.R` lines 74–142 (the banner comment through the closing brace of `load_input`) with exactly:

```r
# =============================================================================
# Input manifest & lazy staging (shared implementation)
# -----------------------------------------------------------------------------
# fig1/inputs.tsv columns: key, group, Assay, local_subpath, gs_uri
# `key`, `local_subpath` and `gs_uri` are reserved; `group` and `Assay` are
# metadata attached to the rows returned by load_input(). Each script loads only
# the keys it needs (see the panel scripts).
# =============================================================================
source("common/manifest.R")
```

Change the header line of `Fig1/inputs.tsv` (line 1) from:

```
key	group	assay	local_subpath	gs_uri
```

to:

```
key	group	Assay	local_subpath	gs_uri
```

Change `Fig1/panels/panelH.R` line 8 from:

```r
  list(annot$path, annot$assay),
```

to:

```r
  list(annot$path, annot$Assay),
```

- [ ] **Step 7: Verify Fig1 still parses and its loaders still resolve**

Run:
```bash
for f in common/manifest.R Fig1/00_setup.R Fig1/panels/*.R Fig1/assemble.R Fig1/run_all.R Fig1/supp_metasoft.R; do
  Rscript -e "invisible(parse('$f'))" || echo "PARSE FAIL: $f"
done
Rscript -e '
  MANIFEST <- "Fig1/inputs.tsv"; INPUT <- "."
  suppressPackageStartupMessages({library(data.table);library(tibble);library(purrr)})
  source("common/manifest.R")
  m <- input_manifest()
  stopifnot(identical(manifest_meta_cols(cbind(m, path = "")), c("group", "Assay")))
  stopifnot(nrow(manifest_rows("perm")) == 21L)
  stopifnot(nrow(manifest_rows("annot")) == 3L)
  cat("Fig1 manifest OK: metadata cols are group/Assay; perm=21 rows; annot=3 rows\n")'
```
Expected: no `PARSE FAIL` lines, then `Fig1 manifest OK: ...`

- [ ] **Step 8: Add the test runner**

Create `tests/run_all.R`:

```r
# Runs every tests/test_*.R in a fresh child process. From the project root:
#   Rscript tests/run_all.R
files <- sort(list.files("tests", pattern = "^test_.*\\.R$", full.names = TRUE))
failed <- character(0)
for (f in files) {
  cat("\n==", f, "==\n")
  if (system2("Rscript", f) != 0L) failed <- c(failed, f)
}
cat("\n")
if (length(failed)) {
  cat("FAILED:\n"); cat(paste0("  ", failed, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("All test files passed.\n")
```

Run: `Rscript tests/run_all.R`
Expected: `All test files passed.`

- [ ] **Step 9: Commit**

```bash
git add common/manifest.R tests/helpers.R tests/test_manifest.R tests/run_all.R \
        Fig1/00_setup.R Fig1/inputs.tsv Fig1/panels/panelH.R
git commit -m "Extract shared manifest loader into common/manifest.R

The manifest/gsutil staging layer was inlined in Fig1/00_setup.R and is needed
byte-identically by Fig3 (and later Fig4/Fig5). Move it to common/manifest.R.

The metadata attach is generalised: every manifest column other than key,
local_subpath and gs_uri is attached to loaded rows under its own name, instead
of the two hardcoded group/assay columns. That moves Fig1's assay -> Assay
rename into its manifest header, so Fig1's loaded frames are unchanged.

Adds tests/ with a regression test for the data.table NSE bug where a manifest
column named 'key' shadowed the 'key' filter argument."
```

---

## Task 2: Fig3 setup skeleton and manifest

**Files:**
- Create: `Fig3/00_setup.R`
- Create: `Fig3/inputs.tsv`
- Create: `Fig3/.gitignore`
- Modify: `.gitignore`

**Interfaces:**
- Consumes: `common/manifest.R` from Task 1.
- Produces: globals `INPUT`, `DATA_DIR`, `OBJ_DIR`, `OUTDIR`, `MANIFEST`, `FONT`, `DPI`, `SZ_AXIS_TITLE`, `SZ_AXIS_TEXT`, `SZ_STRIP`, `SZ_LEG_TEXT`, `SZ_LEG_TITLE`, `SZ_TAG`, `SZ_TITLE`, `SZ_ANNOT`, `blank_plot`, `COLOC_PAIRS`, `LAYER_FILL`, `SIGNAL_FILL`, `AF_POP_COLS`.

- [ ] **Step 1: Create the manifest**

Create `Fig3/inputs.tsv` — **tab-separated**, 1 header + 16 rows:

```
key	type	local_subpath	gs_uri
susie_eqtl	-	finemapping/coloc/eQTL_susie.COMB.annotated.cleaned.tsv.gz	-
susie_sqtl	-	finemapping/coloc/sQTL_susie.COMB.annotated.cleaned.tsv.gz	-
susie_pqtl	-	finemapping/coloc/pQTL_susie.COMB.annotated.cleaned.tsv.gz	-
enloc_gene	eQTL-pQTL	enloc/eQTL_pQTL.enloc.gene.out	-
enloc_gene	sQTL-pQTL	enloc/sQTL_pQTL.enloc.gene.out	-
enloc_gene	eQTL-sQTL	enloc/eQTL_sQTL.enloc.gene.out	-
enloc_cs	eQTL-pQTL	enloc/eQTL_pQTL.enloc.sig.out	-
enloc_cs	sQTL-pQTL	enloc/sQTL_pQTL.enloc.sig.out	-
enloc_cs	eQTL-sQTL	enloc/eQTL_sQTL.enloc.sig.out	-
enloc_variant	eQTL-pQTL	enloc/eQTL_pQTL.enloc.snp.out	-
enloc_variant	sQTL-pQTL	enloc/sQTL_pQTL.enloc.snp.out	-
enloc_variant	eQTL-sQTL	enloc/eQTL_sQTL.enloc.snp.out	-
annotations	-	annotated/AoU.COMB.AC10.AN95.biallelic.v9.annotations.tsv.bgz	gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/checkpoints/AoU.COMB.AC10.AN95.biallelic.v9.annotations.tsv.bgz
protein_lod	-	reference/LOD_proportions_per_proteins.tsv	gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/LOD_proportions_per_proteins.tsv
median_tpm	-	reference/median_tpms.tsv.gz	-
genebayes_shet	-	reference/41588_2024_1820_MOESM4_ESM.xlsx	-
```

- [ ] **Step 2: Create the setup skeleton**

Create `Fig3/00_setup.R`:

```r
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
```

- [ ] **Step 3: Create Fig3/.gitignore**

```
# Per-panel ggplot object cache (regenerated by the pipeline)
objects/

# Memoized derived tables (regenerated by the builders in 00_setup.R)
data/

# macOS
.DS_Store
```

- [ ] **Step 4: Extend the root .gitignore**

In `.gitignore`, add `/enloc/` to the controlled-access inputs block (after `/annotated/`), and add these two lines to the "Figure outputs" block:

```
/Figure3.pdf
/Figure3.png
```

`/finemapping/`, `/annotated/`, `/reference/` and `/FigS_*.pdf` are already present from the Fig1 work.

- [ ] **Step 5: Verify setup sources and the manifest resolves**

Run:
```bash
Rscript -e '
  source("fig3/00_setup.R")
  stopifnot(nrow(input_manifest()) == 16L)
  stopifnot(length(unique(input_manifest()$key)) == 10L)
  stopifnot(nrow(manifest_rows("enloc_cs")) == 3L)
  stopifnot(identical(manifest_meta_cols(cbind(input_manifest(), path = "")), "type"))
  cat("Fig3 setup OK: 16 rows, 10 keys, metadata col = type\n")'
```
Expected: `Fig3 setup OK: 16 rows, 10 keys, metadata col = type`

- [ ] **Step 6: Verify the staged paths and outputs are git-ignored**

Run:
```bash
git check-ignore -v enloc/eQTL_pQTL.enloc.sig.out \
  finemapping/coloc/eQTL_susie.COMB.annotated.cleaned.tsv.gz \
  annotated/AoU.COMB.AC10.AN95.biallelic.v9.annotations.tsv.bgz \
  reference/median_tpms.tsv.gz Figure3.pdf Figure3.png \
  Fig3/objects/A.rds Fig3/data/susie_finemapping.rds
git check-ignore Fig2/objects/A.rds && echo "REGRESSION: Fig2 objects are ignored" || echo "Fig2 objects still tracked - OK"
```
Expected: one match line per path in the first command; second command prints `Fig2 objects still tracked - OK`

- [ ] **Step 7: Commit**

```bash
git add Fig3/00_setup.R Fig3/inputs.tsv Fig3/.gitignore .gitignore
git commit -m "Add Fig3 setup skeleton and input manifest

Fig3/00_setup.R carries the paths, font/DPI, theme and design constants; the
manifest declares 10 keys over 16 rows. Only the AoU annotation table and the
protein LOD table have gs:// origins today; the enloc and SuSiE inputs are
recorded as local-only and are listed as needing a URI in the README."
```

---

## Task 3: FastENLOC readers and FDR labelling

**Files:**
- Modify: `Fig3/00_setup.R` (append)
- Create: `tests/test_fig3_readers.R`

**Interfaces:**
- Consumes: nothing from earlier tasks beyond the setup globals.
- Produces: `ENLOC_CS_COLS`, `ENLOC_VARIANT_COLS`, `LoadEnlocTable(path, cols) -> tibble`, `add_fdr_threshold_label(df, prob_col, control_level = 0.05, label_col = "Coloc", quiet = FALSE) -> df` with an added logical column named by `label_col`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_fig3_readers.R`:

```r
# Tests for Fig3's FastENLOC readers and FDR labelling.
# Run from the project root:  Rscript tests/test_fig3_readers.R
source("fig3/00_setup.R")
source("tests/helpers.R")

tmp <- file.path(tempdir(), "fig3_readers")
unlink(tmp, recursive = TRUE); dir.create(tmp, recursive = TRUE)

# enloc .sig.out: a header line the parser discards, then whitespace-aligned rows.
sig <- file.path(tmp, "pair.enloc.sig.out")
writeLines(c(
  "Signal  Num_SNP  CPIP_qtl  CPIP_gwas_marginal  CPIP_gwas_qtl_prior  RCP  LCP",
  "ENSG001:chr1_100(@)ENSG001:chr1_200   5  0.90  0.10  0.80  0.95  0.90",
  "ENSG002:chr2_100(@)ENSG002:chr2_200   3  0.50  0.20  0.40  0.10  0.05"
), sig)

d <- LoadEnlocTable(sig, ENLOC_CS_COLS)
ok("reader drops the header line",        nrow(d) == 2L)
ok("reader applies the column names",     identical(names(d), ENLOC_CS_COLS))
ok("numeric columns are converted",       is.numeric(d$RCP))
ok("integer columns are converted",       is.numeric(d$Num_SNP))
ok("Signal stays character",              is.character(d$Signal))
ok("values are parsed positionally",      d$RCP[1] == 0.95)

# A .snp.out has six columns, not seven.
snp <- file.path(tmp, "pair.enloc.snp.out")
writeLines(c(
  "Signal  SNP  PIP_qtl  PIP_gwas_marginal  PIP_gwas_qtl_prior  SCP",
  "ENSG001:chr1_100  chr1_100_A_G  0.9  0.1  0.8  0.99",
  "ENSG002:chr2_100  chr2_100_C_T  0.2  0.1  0.1  0.01"
), snp)

v <- LoadEnlocTable(snp, ENLOC_VARIANT_COLS)
ok("variant-level reader parses 6 columns", identical(names(v), ENLOC_VARIANT_COLS))
ok("variant-level SCP is numeric",          is.numeric(v$SCP))

expect_error("wrong column count errors", LoadEnlocTable(snp, ENLOC_CS_COLS))

# FDR labelling: with one near-certain and one near-zero posterior, only the
# former should be labelled at a 5% Bayesian FDR.
lab <- add_fdr_threshold_label(d, "RCP", quiet = TRUE)
ok("adds the Coloc column",              "Coloc" %in% names(lab))
ok("Coloc is logical",                   is.logical(lab$Coloc))
ok("high posterior is labelled TRUE",    isTRUE(lab$Coloc[1]))
ok("low posterior is labelled FALSE",    isFALSE(lab$Coloc[2]))

lab2 <- add_fdr_threshold_label(d, "RCP", label_col = "Sig", quiet = TRUE)
ok("label_col is honoured",              "Sig" %in% names(lab2))

expect_error("an unknown column errors", add_fdr_threshold_label(d, "NOPE", quiet = TRUE))

cat("test_fig3_readers.R: all assertions passed\n")
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript tests/test_fig3_readers.R`
Expected: FAIL — `object 'ENLOC_CS_COLS' not found`

- [ ] **Step 3: Write the implementation**

Append to `Fig3/00_setup.R`:

```r
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
  x <- x[-1]                                   # drop the header line
  x <- x[nzchar(trimws(x))]                    # drop blank lines
  mat <- do.call(rbind, strsplit(trimws(x), "\\s+"))
  if (ncol(mat) != length(cols))
    stop("Expected ", length(cols), " columns in ", path, ", found ", ncol(mat))
  out <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(out) <- cols
  # strsplit() yields character; convert here so downstream code needs no
  # scattered as.numeric() / across(contains('PIP'), as.numeric) calls.
  tibble::as_tibble(type.convert(out, as.is = TRUE))
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

  lfdr  <- sort(1 - probs, na.last = NA)
  FDR   <- cumsum(lfdr) / seq_along(lfdr)
  valid <- which(FDR <= control_level)

  thresh <- if (length(valid) == 0) Inf else 1 - lfdr[max(valid)]

  if (!quiet)
    message(sprintf("    %s >= %.4g at FDR %.0f%%  (%d / %d rows labelled)",
                    prob_col, thresh, 100 * control_level,
                    sum(probs >= thresh, na.rm = TRUE), length(probs)))

  df[[label_col]] <- probs >= thresh
  df
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/test_fig3_readers.R`
Expected: PASS — 14 `ok - ...` lines then `test_fig3_readers.R: all assertions passed`

- [ ] **Step 5: Commit**

```bash
git add Fig3/00_setup.R tests/test_fig3_readers.R
git commit -m "Add Fig3 FastENLOC readers and FDR labelling

Collapses the original's two near-identical loaders (LoadEnlocCSLevel and
LoadEnlocVariantLevel, differing only in column names) into one
LoadEnlocTable(path, cols), and deletes LoadEnlocGeneLevel, which was never
called and declared credible-set column names for gene-level files that are
read with fread.

type.convert() in the reader removes ~10 scattered as.numeric() calls
downstream. The reader now returns explicitly rather than relying on
setnames()'s invisible return, and errors on an unexpected column count.

add_fdr_threshold_label reports its threshold via message(); the original
stored it in an attribute that the following dplyr verbs silently dropped."
```

---

## Task 4: Derived-data cache and the SuSiE builder

**Files:**
- Modify: `Fig3/00_setup.R` (append)
- Create: `tests/fixtures.R`

**Interfaces:**
- Consumes: `input_path()` (Task 1), setup globals (Task 2).
- Produces: `build_or_cache(name, expr, refresh = FALSE)`; `susie_finemapping(refresh = FALSE) -> list(eQTL, sQTL, pQTL, cs_list, gene_map, genes_all_three_assays)`. Each of `$eQTL`/`$sQTL`/`$pQTL` has columns `molecular_trait_id`, `variant`, `pip`, `cs_index`, `posterior_mean`, `posterior_sd`, `cs_id`, `cpip`, `number_variants` (and `gene_name` on `$eQTL`). `$cs_list` has `molecular_trait_id`, `cs_id`, `type`.

- [ ] **Step 1: Write the mock-project fixture**

Create `tests/fixtures.R`:

```r
# Builds a mock Fig3 project in a temp directory: a manifest plus tiny input
# files with the same schemas as the real ones. Sourcing scripts set the FIG3_*
# environment variables from the returned paths BEFORE sourcing fig3/00_setup.R.
#
# Three genes:
#   ENSG001 - fine-mapped in all three assays, colocalized in all three pairs
#   ENSG002 - fine-mapped in all three assays, colocalized in one pair
#   ENSG003 - eQTL only (so it is excluded from genes_all_three_assays)
#
# CREDIBLE-SET ID CONVENTION. A FastENLOC Signal is "<gene>:<cs_id>", and the
# pipeline recovers cs_id by stripping the "^ENSG[0-9]+:" prefix, then joins it
# back to the SuSiE cs_id column. So the SuSiE cs_id must be the BARE id, with
# no gene prefix, or the component join silently matches nothing. The pQTL cs_id
# additionally carries a "<UniProt>_" prefix, which is where the gene/protein map
# comes from.
make_fig3_fixture <- function(root = file.path(tempdir(), "fig3_fixture")) {
  unlink(root, recursive = TRUE)
  for (d in c("finemapping/coloc", "enloc", "annotated", "reference"))
    dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)

  wtsv <- function(df, p) write.table(df, p, sep = "\t", row.names = FALSE, quote = FALSE)

  # Nine variants, alternating between the two shared genes. Enough rows for the
  # per-MAF-bin regressions in panels C and D to actually fit.
  nv     <- 9
  pos    <- seq_len(nv) * 100
  snps   <- sprintf("chr1_%d_A_G", pos)
  genes  <- rep(c("ENSG001", "ENSG002"), length.out = nv)
  uni    <- c(ENSG001 = "P00001", ENSG002 = "P00002")

  # Bare credible-set ids, one per gene per assay.
  e_cs <- paste0("E_", genes, "_L1")
  s_cs <- paste0("S_", genes, "_L1")
  p_cs <- paste0(uni[genes], "_L1")

  # --- SuSiE fine-mapping ---------------------------------------------------
  eqtl <- data.frame(
    molecular_trait_id = paste0(genes, ".1"),        # version suffix is stripped
    variant        = snps,
    pip            = seq(0.9, 0.1, length.out = nv),
    cs_index       = "L1",
    posterior_mean = seq(0.8, 0.1, length.out = nv),
    posterior_sd   = 0.1,
    cs_id          = e_cs,
    gene_name      = paste0("GENE_", genes)
  )
  # An eQTL-only gene, so genes_all_three_assays is a strict subset.
  eqtl <- rbind(eqtl, data.frame(
    molecular_trait_id = "ENSG003.1", variant = "chr3_100_G_A", pip = 0.6,
    cs_index = "L1", posterior_mean = 0.3, posterior_sd = 0.1,
    cs_id = "E_ENSG003_L1", gene_name = "GENE_ENSG003"))
  wtsv(eqtl, file.path(root, "finemapping/coloc/eQTL_susie.COMB.annotated.cleaned.tsv.gz"))

  sqtl <- data.frame(
    molecular_trait_id = paste0(genes, ":clu_1"),    # gene is extracted from this
    variant        = snps,
    pip            = seq(0.8, 0.2, length.out = nv),
    posterior_mean = seq(-0.4, 0.4, length.out = nv),   # signed, to exercise abs()
    cs_index       = "L1",
    posterior_sd   = 0.1,
    cs_id          = s_cs
  )
  wtsv(sqtl, file.path(root, "finemapping/coloc/sQTL_susie.COMB.annotated.cleaned.tsv.gz"))

  pqtl <- data.frame(
    molecular_trait_id = genes,
    variant        = snps,
    pip            = seq(0.9, 0.3, length.out = nv),
    cs_index       = "L1",
    posterior_mean = seq(0.7, 0.2, length.out = nv),
    posterior_sd   = 0.1,
    cs_id          = p_cs,
    group          = "COMB"
  )
  # A non-COMB row that filter(group == "COMB") must drop.
  pqtl <- rbind(pqtl, transform(pqtl[1, ], group = "AFR"))
  wtsv(pqtl, file.path(root, "finemapping/coloc/pQTL_susie.COMB.annotated.cleaned.tsv.gz"))

  # --- FastENLOC gene level (a normal TSV) ----------------------------------
  gene_out <- function(p, glcp) wtsv(
    data.frame(Gene = c("ENSG001", "ENSG002"), GLCP = glcp), p)
  gene_out(file.path(root, "enloc/eQTL_pQTL.enloc.gene.out"), c(0.99, 0.98))
  gene_out(file.path(root, "enloc/sQTL_pQTL.enloc.gene.out"), c(0.99, 0.01))
  gene_out(file.path(root, "enloc/eQTL_sQTL.enloc.gene.out"), c(0.99, 0.01))

  # --- FastENLOC signal level (whitespace-aligned; the header is discarded) --
  # Signal = "<gene>:<left cs_id>(@)<gene>:<right cs_id>".
  sig_out <- function(p, left, right, rcp) writeLines(c(
    "Signal Num_SNP CPIP_qtl CPIP_gwas_marginal CPIP_gwas_qtl_prior RCP LCP",
    sprintf("%s(@)%s 5 0.9 0.1 0.8 %s 0.9", left, right, rcp)), p)
  g2 <- c("ENSG001", "ENSG002")
  sig_out(file.path(root, "enloc/eQTL_pQTL.enloc.sig.out"),
          paste0(g2, ":E_", g2, "_L1"), paste0(g2, ":", uni[g2], "_L1"),
          c("0.99", "0.98"))
  sig_out(file.path(root, "enloc/sQTL_pQTL.enloc.sig.out"),
          paste0(g2, ":S_", g2, "_L1"), paste0(g2, ":", uni[g2], "_L1"),
          c("0.99", "0.01"))
  sig_out(file.path(root, "enloc/eQTL_sQTL.enloc.sig.out"),
          paste0(g2, ":E_", g2, "_L1"), paste0(g2, ":S_", g2, "_L1"),
          c("0.99", "0.01"))

  # --- FastENLOC variant level ---------------------------------------------
  snp_out <- function(p, cs) writeLines(c(
    "Signal SNP PIP_qtl PIP_gwas_marginal PIP_gwas_qtl_prior SCP",
    sprintf("%s:%s %s 0.9 0.1 0.8 0.99", genes, cs, snps)), p)
  snp_out(file.path(root, "enloc/eQTL_pQTL.enloc.snp.out"), e_cs)
  snp_out(file.path(root, "enloc/sQTL_pQTL.enloc.snp.out"), s_cs)
  snp_out(file.path(root, "enloc/eQTL_sQTL.enloc.snp.out"), e_cs)

  # --- AoU annotation table -------------------------------------------------
  # Three variants per MAF bin so each panel-D facet can fit a slope.
  # Variant 1 has gvs_eur_af = 0.90, which must FOLD to 0.10, leaving max_af at
  # 0.30 from gvs_afr_af — the assertion in tests/test_fig3_builders.R.
  afr <- c(0.30, 0.20, 0.10, 0.030, 0.020, 0.015, 0.005, 0.003, 0.002)
  eur <- afr; eur[1] <- 0.90
  ann <- data.frame(
    chrom = "chr1", pos = pos, ref = "A", alt = "G",
    gvs_afr_af = afr,
    gvs_amr_af = afr * 0.9,
    gvs_eas_af = afr * 0.8,
    gvs_eur_af = eur,
    gvs_mid_af = afr * 0.7,
    gvs_sas_af = afr * 0.6,
    gvs_oth_af = afr * 0.5,
    gvs_all_af = afr,
    gvs_max_af = pmax(afr, eur)
  )
  wtsv(ann, file.path(root, "annotated/AoU.COMB.AC10.AN95.biallelic.v9.annotations.tsv.bgz"))

  # --- Manifest -------------------------------------------------------------
  manifest <- file.path(root, "inputs.tsv")
  rows <- read.delim("Fig3/inputs.tsv", colClasses = "character")
  rows <- rows[!rows$key %in% c("protein_lod", "median_tpm", "genebayes_shet"), ]
  write.table(rows, manifest, sep = "\t", row.names = FALSE, quote = FALSE)

  list(root = root, manifest = manifest)
}

# Point fig3/00_setup.R at the fixture. Call BEFORE sourcing it.
use_fig3_fixture <- function(fx) {
  Sys.setenv(FIG3_INPUT    = fx$root,
             FIG3_MANIFEST = fx$manifest,
             FIG3_DATA     = file.path(fx$root, "data"),
             FIG3_OBJ      = file.path(fx$root, "objects"),
             FIG3_OUTDIR   = fx$root)
}
```

- [ ] **Step 2: Write the failing test**

Create `tests/test_fig3_builders.R` (extended again in Task 7):

```r
# Tests for Fig3's derived-data builders against the mock fixture.
# Run from the project root:  Rscript tests/test_fig3_builders.R
source("tests/fixtures.R")
fx <- make_fig3_fixture()
use_fig3_fixture(fx)

source("fig3/00_setup.R")
source("tests/helpers.R")

# --- build_or_cache -------------------------------------------------------
calls <- 0
counted <- function(refresh = FALSE)
  build_or_cache("counted", refresh = refresh, { calls <<- calls + 1; calls })

ok("first call evaluates the expression", counted() == 1L)
ok("second call reads the cache",         counted() == 1L)
ok("the expression ran exactly once",     calls == 1L)
ok("refresh = TRUE rebuilds",             counted(refresh = TRUE) == 2L)
ok("cache file is written to DATA_DIR",   file.exists(file.path(DATA_DIR, "counted.rds")))

# --- susie_finemapping ----------------------------------------------------
fm <- susie_finemapping()

ok("returns all three assays plus derived tables",
   all(c("eQTL", "sQTL", "pQTL", "cs_list", "gene_map", "genes_all_three_assays") %in% names(fm)))
ok("eQTL trait ids are version-stripped",
   setequal(fm$eQTL$molecular_trait_id, c("ENSG001", "ENSG002", "ENSG003")))
ok("sQTL trait ids reduce to the gene",
   setequal(fm$sQTL$molecular_trait_id, c("ENSG001", "ENSG002")))
ok("pQTL drops the non-COMB row",          nrow(fm$pQTL) == 9L)
# Row 1 is ENSG001, whose credible set holds the 5 odd-indexed variants
# (pip 0.9, 0.7, 0.5, 0.3, 0.1).
ok("cpip sums pip within a credible set",  abs(fm$eQTL$cpip[1] - 2.5) < 1e-9)
ok("number_variants counts the credible set", fm$eQTL$number_variants[1] == 5L)
ok("cs_list labels every assay",           setequal(fm$cs_list$type, c("eQTL", "pQTL", "sQTL")))

# genes_all_three_assays is an INTERSECTION; the original misnamed it QTL_union_gene_list.
ok("genes_all_three_assays is the three-way intersection",
   setequal(fm$genes_all_three_assays, c("ENSG001", "ENSG002")))
ok("eQTL-only genes are excluded",         !"ENSG003" %in% fm$genes_all_three_assays)
ok("gene_map carries the UniProt id",      "P00001" %in% fm$gene_map$UniProt)

cat("test_fig3_builders.R: all assertions passed\n")
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `Rscript tests/test_fig3_builders.R`
Expected: FAIL — `could not find function "build_or_cache"`

- [ ] **Step 4: Write the implementation**

Append to `Fig3/00_setup.R`:

```r
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
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `Rscript tests/test_fig3_builders.R`
Expected: PASS — 15 `ok - ...` lines then `test_fig3_builders.R: all assertions passed`

- [ ] **Step 6: Commit**

```bash
git add Fig3/00_setup.R tests/fixtures.R tests/test_fig3_builders.R
git commit -m "Add Fig3 derived-data cache and SuSiE builder

build_or_cache() memoizes each derived table to fig3/data/<name>.rds, so panels
can declare only the inputs they need without re-running the shared upstream.
Its expr argument is a lazy promise, so a cache hit never touches the inputs.

susie_finemapping() folds the three cleaned SuSiE tables plus the credible-set
list, gene/UniProt map and three-way gene intersection into one builder. That
intersection was named QTL_union_gene_list in the original despite being an
intersection; it is now genes_all_three_assays.

tests/fixtures.R builds a mock Fig3 project so the builders are testable without
controlled-access data."
```

---

## Task 5: The three FastENLOC colocalization builders

**Files:**
- Modify: `Fig3/00_setup.R` (append)

**Interfaces:**
- Consumes: `load_input()`, `input_path()`, `LoadEnlocTable`, `ENLOC_CS_COLS`, `ENLOC_VARIANT_COLS`, `add_fdr_threshold_label`, `build_or_cache`.
- Produces:
  - `coloc_gene_level(refresh = FALSE)` → columns `Gene`, `GLCP`, `type`, `Coloc`
  - `coloc_cs_level(refresh = FALSE)` → `ENLOC_CS_COLS` + `type`, `Coloc`, `cs_1`, `cs_2` (`Signal` retained)
  - `coloc_variant_level(refresh = FALSE)` → `ENLOC_VARIANT_COLS` + `type`, `Coloc`, `Gene`, and every annotation column

- [ ] **Step 1: Write the implementation**

Append to `Fig3/00_setup.R`:

```r
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
    load_input("enloc_cs", reader = LoadEnlocTable, cols = ENLOC_CS_COLS) %>%
      .label_by_pair("RCP", "signal level") %>%
      separate(Signal, into = c("cs_1", "cs_2"), sep = "\\(@\\)",
               remove = FALSE, fill = "right", extra = "merge") %>%
      mutate(cs_1 = str_remove(cs_1, "^ENSG[0-9]+:"),
             cs_2 = str_remove(cs_2, "^ENSG[0-9]+:"))
  })

coloc_variant_level <- function(refresh = FALSE)
  build_or_cache("coloc_variant_level", refresh = refresh, {
    annotations <- fread(input_path("annotations")) %>%
      mutate(SNP = paste(chrom, pos, ref, alt, sep = "_"))

    load_input("enloc_variant", reader = LoadEnlocTable, cols = ENLOC_VARIANT_COLS) %>%
      .label_by_pair("SCP", "variant level") %>%
      mutate(Gene = str_remove(Signal, ":.*")) %>%
      left_join(annotations, by = "SNP")
  })
```

- [ ] **Step 2: Verify the builders run against the fixture**

Run:
```bash
Rscript -e '
  source("tests/fixtures.R"); fx <- make_fig3_fixture(); use_fig3_fixture(fx)
  source("fig3/00_setup.R")

  g <- coloc_gene_level()
  stopifnot(nrow(g) == 6L, all(c("Gene","GLCP","type","Coloc") %in% names(g)))
  stopifnot(setequal(g$type, c("eQTL-pQTL","sQTL-pQTL","eQTL-sQTL")))

  cs <- coloc_cs_level()
  stopifnot(nrow(cs) == 6L, all(c("cs_1","cs_2","Coloc","type") %in% names(cs)))
  stopifnot(all(!grepl("^ENSG[0-9]+:", cs$cs_1)))   # gene prefix stripped
  stopifnot(all(!grepl("^ENSG[0-9]+:", cs$cs_2)))

  v <- coloc_variant_level()
  stopifnot(nrow(v) == 27L, "Gene" %in% names(v), "gvs_eur_af" %in% names(v))
  stopifnot(all(v$Gene %in% c("ENSG001","ENSG002")))
  stopifnot(sum(is.na(v$gvs_eur_af)) == 0)          # every SNP joined
  cat("coloc builders OK: gene=6 cs=6 variant=27 rows, annotation joined\n")'
```
Expected: `coloc builders OK: gene=6 cs=6 variant=27 rows, annotation joined`

- [ ] **Step 3: Commit**

```bash
git add Fig3/00_setup.R
git commit -m "Add Fig3 gene/signal/variant colocalization builders

Each reads its three per-pair FastENLOC outputs through the manifest, labels
colocalization at a 5% Bayesian FDR per assay pair (matching the original,
which labelled each file before binding), and caches the result.

The assay-pair label now arrives attached from the manifest, removing nine
mutate(type = '...') calls."
```

---

## Task 6: Credible-set edges and components — the two data-corrupting fixes

**Files:**
- Modify: `Fig3/00_setup.R` (append)
- Create: `tests/test_fig3_edges.R`

**Interfaces:**
- Consumes: `coloc_cs_level()`, `susie_finemapping()`, `build_or_cache`.
- Produces: `coloc_cs_edges(cs_colocs) -> tibble(gene, RCP, LCP, type, cs_1, cs_2)` (a **pure function** of its argument, so it is testable without any manifest); `coloc_cs_components(refresh = FALSE) -> tibble(Gene, cs_id, type, coloc_component, has_coloc)`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_fig3_edges.R`:

```r
# Regression tests for the two data-corrupting defects in Fig3.R.
#
#   Fig3.R:462  str_extract(sig_1, "^EN[<43;65;42MSG[0-9]+") - a pasted terminal
#               escape sequence inside the regex. gene_1 was NA for every
#               eQTL-pQTL row, so filter(gene_1 == gene_2) dropped all of them.
#   Fig3.R:477  bind_rows() listed cleaned_eqtl_pqtl twice and never bound
#               cleaned_eqtl_sqtl.
#
# Together those meant the credible-set graph behind panel B was built from
# sQTL-pQTL edges only.
#
# coloc_cs_edges() is a pure function, so no manifest or fixture is needed.
# Run from the project root:  Rscript tests/test_fig3_edges.R
source("fig3/00_setup.R")
source("tests/helpers.R")

cs_colocs <- tibble::tribble(
  ~type,       ~Signal,                                                ~cs_1,      ~cs_2,      ~RCP, ~LCP, ~Coloc,
  "eQTL-pQTL", "ENSG001:chr1_100(@)ENSG001:chr1_200",                  "chr1_100", "chr1_200", 0.90, 0.90, TRUE,
  "eQTL-sQTL", "ENSG002:chr2_100(@)ENSG002:chr2_200",                  "chr2_100", "chr2_200", 0.80, 0.80, TRUE,
  "sQTL-pQTL", "ENSG003:chr3_100(@)ENSG003:chr3_200",                  "chr3_100", "chr3_200", 0.70, 0.70, TRUE,
  "eQTL-pQTL", "ENSG004:chr4_100(@)ENSG005:chr4_200",                  "chr4_100", "chr4_200", 0.60, 0.60, TRUE,
  "eQTL-pQTL", "ENSG006:chr6_100(@)ENSG006:chr6_200",                  "chr6_100", "chr6_200", 0.50, 0.50, FALSE,
  "eQTL-sQTL", "ENSG007:chr7_100(@)ENSG007:chr7_200_ENSG007:chr7_300", "chr7_100", NA,         0.95, 0.95, TRUE
)

e <- coloc_cs_edges(cs_colocs)

# Defect 2: every assay pair must contribute, exactly once.
ok("all three assay pairs are represented",
   setequal(e$type, c("eQTL-pQTL", "eQTL-sQTL", "sQTL-pQTL")))
ok("eQTL-sQTL is present (was dropped by the duplicated bind_rows)",
   "ENSG002" %in% e$gene)
ok("no pair is bound twice", sum(e$gene == "ENSG003") == 1L)

# Defect 1: the eQTL-pQTL gene regex must match.
ok("eQTL-pQTL edges survive the gene regex", "ENSG001" %in% e$gene)

# Behaviour that must be preserved.
ok("signals spanning two genes are dropped",  !"ENSG004" %in% e$gene)
ok("non-colocalized signals are dropped",     !"ENSG006" %in% e$gene)
ok("a multi-signal second side is unnested",  sum(e$gene == "ENSG007") == 2L)
ok("unnested credible sets lose the gene prefix",
   setequal(e$cs_2[e$gene == "ENSG007"], c("chr7_200", "chr7_300")))
ok("edges carry gene/type/cs columns",
   identical(sort(names(e)), sort(c("gene", "RCP", "LCP", "type", "cs_1", "cs_2"))))
ok("total edge count is as expected", nrow(e) == 5L)

cat("test_fig3_edges.R: all assertions passed\n")
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript tests/test_fig3_edges.R`
Expected: FAIL — `could not find function "coloc_cs_edges"`

- [ ] **Step 3: Write the implementation**

Append to `Fig3/00_setup.R`:

```r
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
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/test_fig3_edges.R`
Expected: PASS — 10 `ok - ...` lines then `test_fig3_edges.R: all assertions passed`

- [ ] **Step 5: Verify components build against the fixture**

Run:
```bash
Rscript -e '
  source("tests/fixtures.R"); fx <- make_fig3_fixture(); use_fig3_fixture(fx)
  source("fig3/00_setup.R")
  cs <- coloc_cs_components()
  stopifnot(all(c("Gene","cs_id","type","coloc_component","has_coloc") %in% names(cs)))
  # 7 credible sets in the fixture; 5 join into 2 components. If this is 0, the
  # component ids did not join to the SuSiE cs_id column.
  stopifnot(nrow(cs) == 7L, sum(cs$has_coloc) == 5L)
  stopifnot(length(unique(na.omit(cs$coloc_component))) == 2L)
  cat("components OK:", sum(cs$has_coloc), "of", nrow(cs), "credible sets colocalized\n")'
```
Expected: `components OK: 5 of 7 credible sets colocalized`

- [ ] **Step 6: Commit**

```bash
git add Fig3/00_setup.R tests/test_fig3_edges.R
git commit -m "Fix the two defects that reduced Fig3's coloc graph to one pair

Fig3.R:462 held a pasted terminal escape sequence inside a regex,
'^EN[<43;65;42MSG[0-9]+', so gene_1 was NA for every eQTL-pQTL row and the
following filter(gene_1 == gene_2) dropped all of them. Fig3.R:477 then listed
cleaned_eqtl_pqtl twice in bind_rows() and never bound cleaned_eqtl_sqtl.

Together the credible-set graph behind panel B and its inset was built from
sQTL-pQTL edges only. Panel B is expected to change materially.

coloc_cs_edges() is now a pure function of its argument, so tests/test_fig3_edges.R
regression-tests both defects without needing controlled-access data."
```

---

## Task 7: The variant-level model input builder

**Files:**
- Modify: `Fig3/00_setup.R` (append)
- Modify: `tests/test_fig3_builders.R` (append)

**Interfaces:**
- Consumes: `coloc_variant_level()`, `susie_finemapping()`, `AF_POP_COLS`, `build_or_cache`.
- Produces: `coloc_model_input(refresh = FALSE) -> named list` with elements `"eQTL-sQTL"`, `"eQTL-pQTL"`, `"sQTL-pQTL"`. Each element adds `max_af`, `rare`, `log_maf`, and the two joined SuSiE posterior pairs `posterior_mean.x`/`posterior_sd.x` (left assay) and `posterior_mean.y`/`posterior_sd.y` (right assay). Rows are filtered to `Coloc == TRUE`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_fig3_builders.R`, immediately before the final `cat(...)` line:

```r
# --- coloc_model_input ----------------------------------------------------
mi <- coloc_model_input()

ok("returns one element per assay pair",
   setequal(names(mi), c("eQTL-sQTL", "eQTL-pQTL", "sQTL-pQTL")))

ep <- mi[["eQTL-pQTL"]]
ok("adds max_af / rare / log_maf",
   all(c("max_af", "rare", "log_maf") %in% names(ep)))
ok("adds both joined posterior columns",
   all(c("posterior_mean.x", "posterior_mean.y") %in% names(ep)))
ok("keeps only colocalized rows",     all(ep$Coloc))
ok("log_maf is a plain numeric vector, not a matrix", is.null(dim(ep$log_maf)))

# chr1_100 has gvs_eur_af = 0.90, which must fold to a MAF of 0.10 before the
# per-variant maximum is taken.
chr1 <- ep[ep$SNP == "chr1_100_A_G", ]
ok("allele frequencies fold to MAF before pmax", abs(chr1$max_af[1] - 0.30) < 1e-9)
ok("common variants are not flagged rare",       isFALSE(chr1$rare[1]))

sp <- mi[["sQTL-pQTL"]]
ok("rare variants are flagged",
   nrow(sp) == 0L || all(sp$max_af[sp$rare] < 0.01))
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript tests/test_fig3_builders.R`
Expected: FAIL — `could not find function "coloc_model_input"`

- [ ] **Step 3: Write the implementation**

Append to `Fig3/00_setup.R`:

```r
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
      stop("Annotation table is missing allele-frequency columns: ",
           paste(missing_af, collapse = ", "))

    folded <- v %>%
      select(all_of(AF_POP_COLS)) %>%
      mutate(across(everything(), ~ ifelse(. > .5, 1 - ., .)))

    cleaned <- v %>%
      mutate(max_af  = do.call(pmax, c(as.list(folded), na.rm = TRUE)),
             rare    = max_af < .01,
             log_maf = as.numeric(scale(-log10(max_af)))) %>%
      filter(Coloc == TRUE)

    posteriors <- function(d)
      d %>% select(variant, molecular_trait_id, posterior_mean, posterior_sd)

    pair_input <- function(pair, left, right) {
      cleaned %>%
        filter(type == pair) %>%
        left_join(posteriors(left),
                  by = c("SNP" = "variant", "Gene" = "molecular_trait_id")) %>%
        left_join(posteriors(right),
                  by = c("SNP" = "variant", "Gene" = "molecular_trait_id"))
    }

    list(
      `eQTL-sQTL` = pair_input("eQTL-sQTL", fm$eQTL, fm$sQTL),
      `eQTL-pQTL` = pair_input("eQTL-pQTL", fm$eQTL, fm$pQTL),
      `sQTL-pQTL` = pair_input("sQTL-pQTL", fm$sQTL, fm$pQTL)
    )
  })
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/test_fig3_builders.R`
Expected: PASS — 23 `ok - ...` lines then `test_fig3_builders.R: all assertions passed`

- [ ] **Step 5: Commit**

```bash
git add Fig3/00_setup.R tests/test_fig3_builders.R
git commit -m "Add Fig3 variant-level model input builder

Folds annotation allele frequencies to MAF, takes the per-variant maximum across
ancestry groups, and joins each pair's two sets of SuSiE posteriors.

The per-ancestry AF columns are now named explicitly (AF_POP_COLS) rather than
selected by contains('gvs') & contains('af') and then subtracted from, so a
schema change fails loudly. max_af is computed with mutate() instead of the
original's positional bind_cols(), the dead first 'rare' assignment on the
uncorrected frequency is gone, and log_maf is a plain numeric rather than the
1-column matrix scale() returns."
```

---

## Task 8: Panels A and B

**Files:**
- Create: `Fig3/panels/panelA.R`
- Create: `Fig3/panels/panelB.R`

**Interfaces:**
- Consumes: `coloc_gene_level()`, `susie_finemapping()`, `coloc_cs_components()`, `LAYER_FILL`, `SIGNAL_FILL`, `SZ_TITLE`, `SZ_ANNOT`, `OBJ_DIR`.
- Produces: `fig3/objects/A.rds`, `fig3/objects/B.rds`.

- [ ] **Step 1: Write panel A**

Create `Fig3/panels/panelA.R`:

```r
# Panel A — gene-level colocalization: percent of genes by the number of
# colocalized molecular layers, among genes fine-mapped in all three assays.
# Inputs: enloc_gene, susie_eqtl, susie_sqtl, susie_pqtl
source("fig3/00_setup.R")

genes_all_three_assays <- susie_finemapping()$genes_all_three_assays

layer_counts <- coloc_gene_level() %>%
  filter(Gene %in% genes_all_three_assays) %>%
  group_by(Gene) %>%
  summarize(NumColocs = sum(Coloc), .groups = "drop") %>%
  dplyr::count(NumColocs) %>%
  mutate(prop = n / sum(n))

# Computed, never hardcoded: the original pasted "86.2%" into this annotation
# while computing the same value 130 lines earlier.
shared_pct <- layer_counts %>%
  filter(NumColocs > 0) %>%
  summarize(p = sum(prop)) %>%
  pull(p)

message(sprintf("Panel A: %.1f%% of genes share a signal across >=1 layer",
                100 * shared_pct))

pA <- layer_counts %>%
  mutate(NumColocs = factor(NumColocs, levels = 0:3)) %>%
  ggplot(aes(x = NumColocs, y = prop, fill = NumColocs)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = LAYER_FILL, guide = "none") +
    annotate("segment", x = 1.5, xend = 4.5, y = 0.455, yend = 0.455,
             color = "grey50", linewidth = 0.25) +
    annotate("text", x = 3, y = 0.465,
             label = sprintf("%.1f%% shared across ≥1 molecular layer",
                             100 * shared_pct),
             size = SZ_ANNOT, color = "grey20") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.55)) +
    labs(x = "Number of Colocalized Molecular Layers", y = "Percent of Genes",
         title = "Gene-level Colocalization") +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = SZ_TITLE, face = "bold"),
          plot.title.position = "plot")

saveRDS(pA, file.path(OBJ_DIR, "A.rds"))
message("Panel A built.")
```

- [ ] **Step 2: Write panel B**

Create `Fig3/panels/panelB.R`:

```r
# Panel B — signal-level colocalization: percent of genes by the number of
# colocalized signals, with an inset of per-assay credible-set coloc rate.
# Inputs: enloc_cs, susie_eqtl, susie_sqtl, susie_pqtl
source("fig3/00_setup.R")

cs_level        <- coloc_cs_components()
genes_all_three_assays <- susie_finemapping()$genes_all_three_assays

signal_counts <- cs_level %>%
  filter(!is.na(coloc_component)) %>%       # the original used na.omit()
  distinct(Gene, coloc_component) %>%
  dplyr::count(Gene) %>%
  dplyr::count(n, name = "nn")

main <- signal_counts %>%
  ggplot(aes(x = factor(n), y = nn / sum(nn), fill = factor(n))) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = SIGNAL_FILL, guide = "none") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 0.55)) +
    labs(x = "Number of colocalized signals", y = "Percent of genes",
         title = "Signal-level Colocalization") +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = SZ_TITLE, face = "bold"),
          plot.title.position = "plot")

coloc_rate <- cs_level %>%
  filter(Gene %in% genes_all_three_assays) %>%
  group_by(type) %>%
  summarise(n_cs       = n_distinct(cs_id),
            n_coloc    = n_distinct(cs_id[has_coloc]),
            coloc_rate = n_coloc / n_cs,
            .groups = "drop")

inset <- coloc_rate %>%
  ggplot(aes(y = reorder(type, coloc_rate), x = coloc_rate, fill = type)) +
    geom_col(show.legend = FALSE) +
    paletteer::scale_fill_paletteer_d("ltc::minou") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = c(0, 0)) +
    labs(x = "% CS coloc", y = NULL) +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.y = element_line(color = "lightgrey"),
          axis.ticks.y = element_blank(),
          plot.tag = element_blank(),      # an inset is never tagged
          plot.margin = margin(0, 0, 0, 0))

pB <- main + inset_element(inset, left = 0.48, bottom = 0.60,
                           right = 0.98, top = 0.98, align_to = "panel")

saveRDS(pB, file.path(OBJ_DIR, "B.rds"))
message("Panel B built.")
```

- [ ] **Step 3: Build both panels against the fixture**

Run:
```bash
Rscript -e '
  source("tests/fixtures.R"); fx <- make_fig3_fixture(); use_fig3_fixture(fx)
  for (p in c("A","B")) source(sprintf("fig3/panels/panel%s.R", p))
  for (p in c("A","B")) {
    o <- readRDS(file.path(Sys.getenv("FIG3_OBJ"), paste0(p, ".rds")))
    stopifnot(inherits(o, "ggplot") || inherits(o, "patchwork"))
    cat("panel", p, "->", class(o)[1], "\n")
  }'
```
Expected: the computed `Panel A: ...%` message, `Panel A built.`, `Panel B built.`, then a class line for each

- [ ] **Step 4: Verify no panel sets its own tag**

Run: `grep -n "labs(tag" Fig3/panels/*.R; test $? -eq 1 && echo "OK: no panel sets a tag"`
Expected: `OK: no panel sets a tag`

- [ ] **Step 5: Commit**

```bash
git add Fig3/panels/panelA.R Fig3/panels/panelB.R
git commit -m "Add Fig3 panels A and B

Panel A computes the 'shared across >=1 molecular layer' percentage from the
data; the original hardcoded '86.2%' into the annotation while computing the
same value 130 lines earlier.

Panel B composes its coloc-rate inset itself, so the cached object is the
finished panel. Neither panel sets labs(tag = ...) - tagging is assemble.R's
job, which also removes the need for the inset's tag suppression."
```

---

## Task 9: Panels C, D and E

**Files:**
- Create: `Fig3/panels/panelC.R`
- Create: `Fig3/panels/panelD.R`
- Create: `Fig3/panels/panelE.R`

**Interfaces:**
- Consumes: `coloc_model_input()`, `SZ_ANNOT`, `blank_plot`, `OBJ_DIR`.
- Produces: `fig3/objects/C.rds`, `D.rds`, `E.rds`.

- [ ] **Step 1: Write panel C**

Create `Fig3/panels/panelC.R`:

```r
# Panel C — effect of -log(MAF) on cross-assay effect-size correlation:
# the interaction coefficient and 95% CI for each assay pair.
# Inputs: enloc_variant, annotations, susie_eqtl, susie_sqtl, susie_pqtl
source("fig3/00_setup.R")

# The original fit these models TWICE: once on the raw scale without an NA
# filter, then again with log1p() and NA filtering, discarding the first result.
# Only the second is kept.
#
# NOTE: the eQTL-pQTL model omits abs() on both posteriors, unlike the other two
# pairs. That asymmetry is preserved verbatim from the original; see NOTES.md,
# which flags it as possibly unintended (log1p() of a value below -1 is NaN).
FIT_FORMULA <- list(
  `eQTL-sQTL` = log1p(abs(posterior_mean.y)) ~ log1p(abs(posterior_mean.x)) * log_maf,
  `eQTL-pQTL` = log1p(posterior_mean.y)      ~ log1p(posterior_mean.x)      * log_maf,
  `sQTL-pQTL` = log1p(abs(posterior_mean.y)) ~ log1p(abs(posterior_mean.x)) * log_maf
)

fit_pair <- function(d, pair) {
  d <- filter(d, !is.na(posterior_mean.x), !is.na(posterior_mean.y))
  broom::tidy(lm(FIT_FORMULA[[pair]], data = d), conf.int = TRUE) %>%
    filter(!str_detect(term, "Int")) %>%
    mutate(type = pair)
}

interaction_res <- purrr::imap_dfr(coloc_model_input(), ~ fit_pair(.x, .y))

pC <- interaction_res %>%
  filter(str_detect(term, ":")) %>%
  ggplot(aes(x = estimate, y = type)) +
    geom_linerange(aes(xmin = conf.low, xmax = conf.high), color = "grey40") +
    geom_point(size = 4, color = "#364B63") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = "Effect of -log(MAF) on cross-ome correlation", y = NULL) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(l = 60, unit = "pt"))

saveRDS(pC, file.path(OBJ_DIR, "C.rds"))
message("Panel C built.")
```

- [ ] **Step 2: Write panel D**

Create `Fig3/panels/panelD.R`:

```r
# Panel D — eQTL vs pQTL posterior effect size, faceted by MAF bin.
# Inputs: enloc_variant, annotations, susie_eqtl, susie_pqtl
source("fig3/00_setup.R")

binned <- coloc_model_input()[["eQTL-pQTL"]] %>%
  filter(!is.na(max_af)) %>%
  mutate(maf_bin = case_when(max_af >= .05 ~ "Common",
                             max_af >= .01 ~ "Low-Freq",
                             TRUE          ~ "Rare"),
         posterior_mean.x = log1p(posterior_mean.x),
         posterior_mean.y = log1p(posterior_mean.y))

lims <- quantile(abs(c(binned$posterior_mean.x, binned$posterior_mean.y)),
                 0.97, na.rm = TRUE)

slopes <- binned %>%
  group_by(maf_bin) %>%
  summarise(n     = dplyr::n(),
            slope = coef(lm(posterior_mean.y ~ posterior_mean.x))[["posterior_mean.x"]],
            .groups = "drop") %>%
  mutate(label = paste0("β = ", round(slope, 2), "\nN = ", scales::comma(n)))

# Every row here is already Coloc == TRUE (coloc_model_input filters), so the
# original's redundant filter(Coloc == TRUE) inside geom_point is gone.
pD <- ggplot(binned, aes(x = posterior_mean.x, y = posterior_mean.y)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.35, color = "grey65") +
  geom_point(alpha = 0.6, size = 4, color = "#6EC5C3") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.7,
              color = "red", fill = "red", alpha = 0.12) +
  geom_label(data = slopes, aes(x = -Inf, y = Inf, label = label),
             inherit.aes = FALSE, hjust = -0.05, vjust = 1.15,
             size = SZ_ANNOT, label.size = 0, fill = "white", alpha = 0.85) +
  facet_wrap(~maf_bin, nrow = 1) +
  coord_cartesian(xlim = c(-lims, lims), ylim = c(-lims, lims)) +
  labs(x = "Expression Posterior", y = "Protein Posterior") +
  theme(panel.spacing.x = unit(0.8, "lines"),
        plot.margin = margin(5.5, 8, 5.5, 5.5))

saveRDS(pD, file.path(OBJ_DIR, "D.rds"))
message("Panel D built.")
```

- [ ] **Step 3: Write panel E**

Create `Fig3/panels/panelE.R`:

```r
# Panel E — reserved slot. The bottom row of Figure 3 leaves space for a panel
# that has not been decided yet (see fig3/NOTES.md). Rendering a blank keeps the
# layout stable and lets assemble.R stay unconditional.
source("fig3/00_setup.R")

pE <- blank_plot

saveRDS(pE, file.path(OBJ_DIR, "E.rds"))
message("Panel E built (reserved slot).")
```

- [ ] **Step 4: Build all three against the fixture**

Run:
```bash
Rscript -e '
  source("tests/fixtures.R"); fx <- make_fig3_fixture(); use_fig3_fixture(fx)
  for (p in c("C","D","E")) source(sprintf("fig3/panels/panel%s.R", p))
  for (p in c("C","D","E")) {
    o <- readRDS(file.path(Sys.getenv("FIG3_OBJ"), paste0(p, ".rds")))
    stopifnot(inherits(o, "ggplot") || inherits(o, "patchwork"))
    cat("panel", p, "->", class(o)[1], "\n")
  }'
```
Expected: `Panel C built.`, `Panel D built.`, `Panel E built (reserved slot).`, then a class line for each

> The fixture is deliberately tiny, so `lm()` may warn about a near-singular fit
> or an essentially perfect fit. Warnings are fine here; an **error** is not.

- [ ] **Step 5: Commit**

```bash
git add Fig3/panels/panelC.R Fig3/panels/panelD.R Fig3/panels/panelE.R
git commit -m "Add Fig3 panels C, D and the reserved slot E

Panel C keeps only the NA-filtered log1p models; the original fit every model
twice and discarded the first, unfiltered raw-scale result. The eQTL-pQTL
formula's missing abs() is preserved verbatim and flagged in NOTES.md rather
than silently unified with the other two pairs.

Panel D drops the redundant filter(Coloc == TRUE) inside geom_point, since
coloc_model_input already filters.

Panel E renders blank so the layout and assemble.R stay unconditional."
```

---

## Task 10: Assembly, runner and supplementals

**Files:**
- Create: `Fig3/assemble.R`
- Create: `Fig3/run_all.R`
- Create: `Fig3/supp_coloc_qc.R`
- Create: `Fig3/supp_prop_gene_coloc.R`

**Interfaces:**
- Consumes: `fig3/objects/A.rds` … `E.rds`, `coloc_gene_level()`, `coloc_cs_level()`, `coloc_variant_level()`, `SZ_TAG`, `DPI`, `OUTDIR`.
- Produces: `Figure3.pdf`, `Figure3.png`, `FigS_CrossQTLColocQC.pdf`, `FigS_CrossQTLPropGeneColoc.pdf`.

- [ ] **Step 1: Write assemble.R**

Create `Fig3/assemble.R`:

```r
# =============================================================================
# Figure 3 — assembly. Loads the per-panel ggplot objects (fig3/objects/*.rds)
# and composes the final figure. Panels must be built first (fig3/panels/*.R),
# or just run fig3/run_all.R.
#
# Panel E is a reserved, undecided slot (see NOTES.md); it renders blank.
# Tagging happens HERE and only here - no panel sets labs(tag = ...).
# =============================================================================
source("fig3/00_setup.R")

P <- setNames(lapply(LETTERS[1:5],
                     function(x) readRDS(file.path(OBJ_DIR, paste0(x, ".rds")))),
              LETTERS[1:5])

coloc_counts <- P$A + P$B + plot_layout(widths = c(1, 1))
top_row      <- (coloc_counts | P$C) + plot_layout(widths = c(2.5, 1))
bottom_row   <- (P$D | P$E) + plot_layout(widths = c(1.3, 1))

Fig3 <- (top_row / bottom_row) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = SZ_TAG, face = "bold"))

ggsave(file.path(OUTDIR, "Figure3.pdf"), Fig3,
       width = 18, height = 12, dpi = DPI, limitsize = FALSE)
# Raster preview: full physical size at lower DPI (showtext DPI matched).
showtext_opts(dpi = 100)
ggsave(file.path(OUTDIR, "Figure3.png"), Fig3,
       width = 18, height = 12, dpi = 100, limitsize = FALSE)
showtext_opts(dpi = DPI)

message("Assembled: ", file.path(OUTDIR, "Figure3.pdf"), " and Figure3.png")
```

- [ ] **Step 2: Write run_all.R**

Create `Fig3/run_all.R`:

```r
# =============================================================================
# Figure 3 — end-to-end runner. From the project root:  Rscript fig3/run_all.R
#   1. panels  : each stages the inputs it needs (fig3/inputs.tsv) and builds the
#                derived tables it needs (memoized in fig3/data/)
#                -> fig3/objects/<X>.rds
#   2. assemble: compose + save Figure3.pdf / Figure3.png
#
# To rebuild one panel, run just that panel script then fig3/assemble.R.
# The derived-data cache is keyed by NAME, not by content: clear fig3/data/
# after changing an input file.
#
# Supplementals are standalone:
#   Rscript fig3/supp_coloc_qc.R
#   Rscript fig3/supp_prop_gene_coloc.R
# =============================================================================
for (p in LETTERS[1:5]) source(sprintf("fig3/panels/panel%s.R", p))
source("fig3/assemble.R")
message("run_all complete.")
```

- [ ] **Step 3: Write supplemental S1**

Create `Fig3/supp_coloc_qc.R`:

```r
# =============================================================================
# Supplemental — cross-QTL colocalization QC. NOT part of Figure 3.
# Produces FigS_CrossQTLColocQC.pdf:
#   row 1: posterior distributions (GLCP / RCP / SCP) by assay pair
#   row 2: number of colocalizations discovered across posterior thresholds
#
# The original wrote this figure to CrossQTLPosteriorDistributions.pdf twice,
# with different content; the first save was silently overwritten.
#
# Run from the project root:  Rscript fig3/supp_coloc_qc.R
# =============================================================================
source("fig3/00_setup.R")

gene_lvl    <- coloc_gene_level()
cs_lvl      <- coloc_cs_level()
variant_lvl <- coloc_variant_level()

density_panel <- function(d, col, show_legend) {
  ggplot(d, aes(x = .data[[col]], fill = type)) +
    geom_density(alpha = 0.7, show.legend = show_legend) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
}

posterior_dists <- density_panel(gene_lvl,    "GLCP", FALSE) +
                   density_panel(cs_lvl,      "RCP",  FALSE) +
                   density_panel(variant_lvl, "SCP",  TRUE)

THRESHOLDS <- c(0.01, 0.1, 0.5, 0.9)

discovery <- function(d, col, label) {
  purrr::map_dfr(THRESHOLDS, function(t) {
    d %>%
      group_by(type) %>%
      summarise(n_tested = dplyr::n(),
                value    = sum(.data[[col]] > t, na.rm = TRUE),
                .groups  = "drop") %>%
      mutate(threshold = t, posterior = label)
  })
}

discovery_by_type <- bind_rows(
  discovery(gene_lvl,    "GLCP", "GLCP"),
  discovery(cs_lvl,      "RCP",  "RCP"),
  discovery(variant_lvl, "SCP",  "SCP")
)

discovery_plot <- discovery_by_type %>%
  ggplot(aes(x = threshold, y = log10(value), group = type, color = type)) +
    geom_point() +
    geom_line() +
    facet_wrap(~posterior, nrow = 1) +
    labs(x = "Posterior threshold", y = "log10(Number Colocs)")

SuppColocQC <- posterior_dists / discovery_plot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = SZ_TAG, face = "bold"))

ggsave(file.path(OUTDIR, "FigS_CrossQTLColocQC.pdf"), SuppColocQC,
       height = 7, width = 10, dpi = DPI)
message("Wrote ", file.path(OUTDIR, "FigS_CrossQTLColocQC.pdf"))
```

- [ ] **Step 4: Write supplemental S2**

Create `Fig3/supp_prop_gene_coloc.R`:

```r
# =============================================================================
# Supplemental — proportion of genes with a cross-QTL colocalization, by assay
# pair. NOT part of Figure 3. Produces FigS_CrossQTLPropGeneColoc.pdf.
#
# Run from the project root:  Rscript fig3/supp_prop_gene_coloc.R
# =============================================================================
source("fig3/00_setup.R")

ProportionGeneColocs <- coloc_gene_level() %>%
  group_by(type) %>%
  summarize(PropColoc = mean(Coloc), .groups = "drop") %>%
  ggplot(aes(x = type, y = PropColoc)) +
    geom_col(position = "dodge") +
    labs(x = NULL, y = "Proportion Colocalized") +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())

# Explicit plot = argument: the original relied on ggsave()'s last-plot default,
# which is correct in a notebook but fragile under Rscript.
ggsave(file.path(OUTDIR, "FigS_CrossQTLPropGeneColoc.pdf"), ProportionGeneColocs,
       height = 4, width = 6, dpi = DPI)
message("Wrote ", file.path(OUTDIR, "FigS_CrossQTLPropGeneColoc.pdf"))
```

- [ ] **Step 5: Run the whole pipeline end-to-end against the fixture**

Run:
```bash
Rscript -e '
  source("tests/fixtures.R"); fx <- make_fig3_fixture(); use_fig3_fixture(fx)
  source("fig3/run_all.R")
  source("fig3/supp_coloc_qc.R")
  source("fig3/supp_prop_gene_coloc.R")
  out <- Sys.getenv("FIG3_OUTDIR")
  for (f in c("Figure3.pdf","Figure3.png","FigS_CrossQTLColocQC.pdf",
              "FigS_CrossQTLPropGeneColoc.pdf")) {
    p <- file.path(out, f)
    stopifnot(file.exists(p), file.size(p) > 1000)
    cat(sprintf("%-34s %8d bytes\n", f, file.size(p)))
  }
  cat("end-to-end OK\n")'
```
Expected: `run_all complete.`, four size lines all > 1000 bytes, then `end-to-end OK`

- [ ] **Step 6: Verify every ggsave passes an explicit plot**

Run:
```bash
grep -rn "ggsave(" Fig3/*.R | grep -v "Fig3," | grep -v "ProportionGeneColocs" | grep -v "SuppColocQC" \
  && echo "FAIL: a ggsave relies on the last-plot default" || echo "OK: every ggsave passes an explicit plot"
```
Expected: `OK: every ggsave passes an explicit plot`

- [ ] **Step 7: Commit**

```bash
git add Fig3/assemble.R Fig3/run_all.R Fig3/supp_coloc_qc.R Fig3/supp_prop_gene_coloc.R
git commit -m "Add Fig3 assembly, runner and the two supplementals

The original built the assembled figure and never saved it; assemble.R now
writes Figure3.pdf plus a PNG preview at a matched showtext DPI.

Tagging lives here alone, via plot_annotation(tag_levels = 'A').

The QC supplemental replaces the original's hardcoded n_0.01/n_0.1/n_0.5/n_0.9
columns plus pivot_longer with a threshold loop, and is written once rather than
twice to the same filename with different content. Both supplementals pass an
explicit plot to ggsave instead of relying on the last-plot default.

The gsutil uploads of finished figures are dropped, matching Fig1; the README
documents the upload command."
```

---

## Task 11: Documentation, and retire Fig3.R

**Files:**
- Create: `Fig3/README.md`
- Create: `Fig3/NOTES.md`
- Delete: `Fig3/Fig3.R`

- [ ] **Step 1: Write NOTES.md**

Create `Fig3/NOTES.md`:

```markdown
# Figure 3 — deferred work and open questions

## Panel E — reserved slot (undecided)

The bottom row of Figure 3 reserves a slot beside panel D. Its content has not
been decided. `panels/panelE.R` renders a blank so the layout and `assemble.R`
stay unconditional.

Mechanics for whoever fills it: the bottom row is `(D | E)` with
`plot_layout(widths = c(1.3, 1))`, so E gets roughly 43% of an 18" canvas —
about 7.7 × 6 inches. Replace the body of `panels/panelE.R` and nothing else
needs to change.

## Colocalization covariate analysis (dropped from the pipeline)

`Fig3.R:512` built a per-gene table named `coloc_summary` that was never
plotted or written. It joined, per gene:

- credible-set colocalization rate (`PropColoc`, `NumColoc`, `Total`) from the
  component collapse,
- median expression (`median_tpms.tsv.gz` → `mean_TPM`),
- GeneBayes s_het selection estimates (`post_mean`, plus quintile `shet_q` and
  the four `shet_bin` categories) from the Nature Genetics supplement
  `41588_2024_1820_MOESM4_ESM.xlsx`, sheet 2,
- protein detectability (`LOD_proportions_per_proteins.tsv`), joined via the
  gene → UniProt map,

filtered to genes with finite `mean_TPM` and `post_mean`.

The evident question: **does cross-assay colocalization rate track expression
level, selection pressure, or protein detectability?** The code is dropped;
this note is the record. To resume, rebuild it on top of
`coloc_cs_components()` and `susie_finemapping()$gene_map`.

Manifest keys `median_tpm`, `genebayes_shet` and `protein_lod` are kept for
provenance but have **no consumer in the current pipeline**. Note that removing
this analysis also removed both `left_join(genebayes, ...)` calls (`Fig3.R:232`
and `:310`) — no retained panel reads the s_het columns they attached.

Also dropped as unused: `cs_count`, `n_eQTLs`, `n_sQTLs`, `n_pQTLs`,
`rare_gene_list`, and the exploratory `glm(Coloc ~ gvs_max_af + PIP_qtl +
PIP_gwas_marginal)` at `Fig3.R:719-724`.

## Open question — the eQTL-pQTL model formula

`panels/panelC.R` fits one model per assay pair. Two of them use `abs()`:

```
eQTL-sQTL:  log1p(abs(posterior_mean.y)) ~ log1p(abs(posterior_mean.x)) * log_maf
sQTL-pQTL:  log1p(abs(posterior_mean.y)) ~ log1p(abs(posterior_mean.x)) * log_maf
eQTL-pQTL:  log1p(posterior_mean.y)      ~ log1p(posterior_mean.x)      * log_maf
```

The eQTL-pQTL formula omits it. This is preserved **verbatim** from the
original rather than silently unified, because it changes the estimate: SuSiE
posterior means are signed, and `log1p()` of a value below −1 is `NaN`, so any
eQTL-pQTL variant with a posterior mean below −1 is dropped from that fit but
retained in the other two.

If the asymmetry was unintended, add `abs()` to the `eQTL-pQTL` entry of
`FIT_FORMULA` in `panels/panelC.R`.

## Inputs still missing a `gs://` origin

`inputs.tsv` records these as local-only; they must be pre-staged, or given a
URI so the pipeline can self-stage on Terra:

| key | files |
|---|---|
| `susie_eqtl`, `susie_sqtl`, `susie_pqtl` | the three cleaned COMB SuSiE tables |
| `enloc_gene`, `enloc_cs`, `enloc_variant` | the nine FastENLOC outputs |
| `median_tpm`, `genebayes_shet` | deferred analysis above |

In the original these were hardcoded to `/home/jupyter/multiomics_qtls/edit/FastEnlocInput/`
or to bare working-directory filenames.
```

- [ ] **Step 2: Write README.md**

Create `Fig3/README.md`:

````markdown
# Figure 3 — cross-QTL colocalization across molecular layers

Code to generate **Figure 3** (and two supplemental figures) for the manuscript.
The figure asks how often a genetic signal is **shared between molecular layers**
— expression (eQTL), protein (pQTL) and splicing (sQTL) — using FastENLOC
colocalization between each pair of assays at gene, credible-set (signal) and
variant level, and whether that sharing depends on allele frequency.

> **Data note.** The inputs derive from *All of Us* **controlled-access** genomic
> data and are **not** in this repository. See [Data & inputs](#data--inputs).
> Never commit staged inputs or figure outputs.

---

## The figure

**Main figure** (`Figure3.pdf` / `Figure3.png`), 18 × 12″:

| Panel | Content | Source |
|-------|---------|--------|
| **A** | Percent of genes by number of colocalized molecular layers (0–3), among genes fine-mapped in all three assays | gene-level FastENLOC + SuSiE |
| **B** | Percent of genes by number of colocalized signals, inset with per-assay credible-set colocalization rate | signal-level FastENLOC + SuSiE |
| **C** | Effect of −log(MAF) on cross-assay effect-size correlation, per assay pair | variant-level FastENLOC + annotations + SuSiE |
| **D** | eQTL vs pQTL posterior effect size, faceted by MAF bin | variant-level FastENLOC + annotations + SuSiE |
| **E** | *Reserved slot — content undecided (see [`NOTES.md`](NOTES.md))* | — |

**Supplemental figures:**

- `FigS_CrossQTLColocQC.pdf` — GLCP/RCP/SCP posterior distributions by assay
  pair, over colocalization counts across posterior thresholds.
- `FigS_CrossQTLPropGeneColoc.pdf` — proportion of genes colocalized, by pair.

---

## Repository layout

```
.
├── common/manifest.R            # shared input manifest + gsutil staging (all figures)
├── fig3/
│   ├── 00_setup.R               #   libraries, paths, theme, FastENLOC readers,
│   │                            #     FDR labelling, derived-data builders
│   ├── inputs.tsv               #   input manifest: key, type, local_subpath, gs_uri
│   ├── panels/panelA.R … panelE.R  #  each -> fig3/objects/<X>.rds
│   ├── assemble.R               #   composes panels -> Figure3.pdf / Figure3.png
│   ├── run_all.R                #   panels -> assemble (end-to-end)
│   ├── supp_coloc_qc.R          #   supplemental: posteriors + discovery thresholds
│   ├── supp_prop_gene_coloc.R   #   supplemental: proportion of genes colocalized
│   ├── NOTES.md                 #   deferred analyses, open questions
│   ├── data/                    #   memoized derived tables (.rds)  [git-ignored]
│   └── objects/                 #   per-panel ggplot objects (.rds) [git-ignored]
├── enloc/ finemapping/ annotated/ reference/   # inputs staged here [git-ignored]
└── Figure3.pdf / Figure3.png
```

### Two caches, and how they differ

**Staged inputs** (`common/manifest.R`) — `inputs.tsv` maps each logical input
(a `key`) to one or more files: a `local_subpath` relative to the staging root
(the cwd by default) and, where known, a `gs://` origin. A file is fetched with
`gsutil` **only when the local copy is missing**, so pre-staged files skip the
copy entirely. Panels reach inputs through `input_path(key)`,
`load_input(key)` and `input_files(key)`.

**Derived tables** (`fig3/data/`) — Fig3's panels share a heavy upstream: the
AoU annotation table and three SuSiE tables. Six builders in `00_setup.R`
memoize their results, so a panel still declares only what it needs while the
expensive work runs once:

| Builder | Used by |
|---|---|
| `susie_finemapping()` | A, B, C, D |
| `coloc_gene_level()` | A, S1, S2 |
| `coloc_cs_level()` | B, S1 |
| `coloc_variant_level()` | C, D, S1 |
| `coloc_cs_components()` | B |
| `coloc_model_input()` | C, D |

> **This cache is keyed by name, not by content.** After changing an input file
> or a builder, clear `fig3/data/` (or pass `refresh = TRUE`), or you will get a
> stale table back. The same applies to the `objects/` panel cache.

---

## Running

Run everything from the **project root**.

```bash
# Full pipeline: build every panel (staging inputs as needed) -> assemble
Rscript fig3/run_all.R

# Supplementals (standalone)
Rscript fig3/supp_coloc_qc.R
Rscript fig3/supp_prop_gene_coloc.R
```

Rebuild a single panel after editing it, then recompose:

```bash
Rscript fig3/panels/panelD.R
Rscript fig3/assemble.R
```

Force a full rebuild of the derived tables:

```bash
rm -rf fig3/data && Rscript fig3/run_all.R
```

Upload finished figures (the scripts no longer do this themselves):

```bash
gsutil cp Figure3.pdf gs://fc-secure-.../PaperFigures/
gsutil cp FigS_*.pdf  gs://fc-secure-.../PaperFigures/SupplementalFigures/CrossQTL/
```

### Tests

`tests/` holds plain `Rscript` checks that need no controlled-access data —
they run against mock inputs built by `tests/fixtures.R`:

```bash
Rscript tests/run_all.R
```

`tests/test_fig3_edges.R` is the regression test for the two defects that
reduced the credible-set colocalization graph to a single assay pair.

---

## Requirements

- **R** ≥ 4.3
- CRAN: `tidyverse`, `data.table`, `patchwork`, `igraph`, `broom`, `scales`,
  `showtext`, `sysfonts`, `paletteer`

```r
install.packages(c("tidyverse","data.table","patchwork","igraph","broom",
                   "scales","showtext","sysfonts","paletteer"))
```

Text uses **Inter**, fetched from Google Fonts at runtime via `showtext`, with
an automatic fallback to the default sans font. Output DPI is matched to the
`showtext` DPI so text scales correctly.

---

## Data & inputs

Every input is declared in **`fig3/inputs.tsv`**, which is also the record of
each file's origin. Columns: `key`, `type` (assay pair, `-` where not
applicable), `local_subpath` (relative to the staging root, overridable via
`FIG3_INPUT`), `gs_uri` (`-` = local-only).

| Manifest key | Files | Used by |
|---|---|---|
| `susie_eqtl` / `susie_sqtl` / `susie_pqtl` | `finemapping/coloc/{e,s,p}QTL_susie.COMB.annotated.cleaned.tsv.gz` | A, B, C, D |
| `enloc_gene` | `enloc/*.enloc.gene.out` (3 pairs) | A, S1, S2 |
| `enloc_cs` | `enloc/*.enloc.sig.out` (3 pairs) | B, S1 |
| `enloc_variant` | `enloc/*.enloc.snp.out` (3 pairs) | C, D, S1 |
| `annotations` | `annotated/AoU.COMB.AC10.AN95.biallelic.v9.annotations.tsv.bgz` | C, D |
| `protein_lod` / `median_tpm` / `genebayes_shet` | `reference/…` | *(deferred — see NOTES.md)* |

**Only `annotations` and `protein_lod` currently carry a `gs://` URI.** The rest
are local-only and must be pre-staged; a missing one fails with a message naming
the file and the key. See [`NOTES.md`](NOTES.md) for the full list.

The nine FastENLOC outputs are treated as **inputs**, not pipeline products —
this figure does not run FastENLOC.

To add or re-point an input, edit `fig3/inputs.tsv`; no script changes needed.

Access to the underlying *All of Us* data is via the *All of Us* Researcher
Workbench under its Data Use Agreement; it is not shared through this repository.

---

## Provenance

This directory replaces a single 858-line notebook export (`Fig3.R`). The
restructuring and the 21 defects it fixed are recorded in
[`docs/superpowers/specs/2026-07-23-fig3-cleanup-design.md`](../docs/superpowers/specs/2026-07-23-fig3-cleanup-design.md).

**Two of those defects changed results.** A pasted terminal escape sequence in a
regex silently emptied the eQTL-pQTL credible-set frame, and a duplicated
`bind_rows()` argument dropped eQTL-sQTL entirely — so the colocalization graph
behind **panel B and its inset** was built from sQTL-pQTL edges alone. Panel B
is expected to differ from the pre-cleanup version; that difference is the fix,
not a regression.
````

- [ ] **Step 3: Delete the superseded script**

Run:
```bash
/bin/rm -f Fig3/Fig3.R && test ! -e Fig3/Fig3.R && echo "Fig3.R removed"
```
Expected: `Fig3.R removed`

> `rm` is aliased to `rm -i` in this shell; `/bin/rm -f` bypasses the prompt.
> Verify the file is actually gone — a silent prompt has left it in place before.

- [ ] **Step 4: Run the full test suite and a final end-to-end build**

Run:
```bash
Rscript tests/run_all.R
for f in common/manifest.R Fig3/00_setup.R Fig3/panels/*.R Fig3/assemble.R \
         Fig3/run_all.R Fig3/supp_coloc_qc.R Fig3/supp_prop_gene_coloc.R; do
  Rscript -e "invisible(parse('$f'))" || echo "PARSE FAIL: $f"
done
Rscript -e '
  source("tests/fixtures.R"); fx <- make_fig3_fixture(); use_fig3_fixture(fx)
  source("fig3/run_all.R"); cat("final end-to-end OK\n")'
```
Expected: `All test files passed.`, no `PARSE FAIL` lines, `final end-to-end OK`

- [ ] **Step 5: Confirm nothing sensitive or generated is staged**

Run:
```bash
git status --short
git check-ignore Fig3/objects/A.rds Fig3/data/susie_finemapping.rds && echo "caches ignored"
```
Expected: only the intended source files listed; `caches ignored`

- [ ] **Step 6: Commit**

```bash
git add Fig3/README.md Fig3/NOTES.md
git rm --cached Fig3/Fig3.R 2>/dev/null || true
git add -u Fig3
git commit -m "Document Fig3 and retire the notebook export

README covers the layout, the two distinct caches (staged inputs vs memoized
derived tables) and their invalidation rule, how to run and rebuild, and which
inputs still need a gs:// URI. It flags that panel B is expected to differ from
the pre-cleanup version because of the two data-corrupting defects.

NOTES.md records the dropped colocalization covariate analysis, the undecided
panel E slot and its layout mechanics, and the eQTL-pQTL model formula's missing
abs() - preserved verbatim rather than silently unified, since it changes the
estimate.

Deletes Fig3/Fig3.R, now fully superseded."
```

---

## Final verification

- [ ] **Push and open the PR**

```bash
git push -u origin fig3-cleanup
gh pr create --base main --head fig3-cleanup \
  --title "Clean up Figure 3: modular coloc pipeline with shared manifest loaders" \
  --body "$(cat <<'BODY'
Replaces the 858-line notebook export `Fig3/Fig3.R` with a Fig2-style modular
pipeline, and extracts the manifest/gsutil staging layer into `common/manifest.R`
so Fig1, Fig3 and later figures share one implementation.

Spec: `docs/superpowers/specs/2026-07-23-fig3-cleanup-design.md`

## Two defects changed results

- `Fig3.R:462` held a pasted terminal escape sequence inside a regex
  (`"^EN[<43;65;42MSG[0-9]+"`), so `gene_1` was `NA` for every eQTL-pQTL row and
  the following `filter(gene_1 == gene_2)` dropped all of them.
- `Fig3.R:477` listed `cleaned_eqtl_pqtl` twice in `bind_rows()` and never bound
  `cleaned_eqtl_sqtl`.

Together the credible-set graph behind **panel B and its inset** was built from
sQTL-pQTL edges alone. **Panel B is expected to differ from the previous
version** — please eyeball it. `tests/test_fig3_edges.R` regression-tests both.

Also fixed: a hardcoded `"86.2%"` annotation that is now computed; the assembled
figure was never saved; two `ggsave()` calls relied on the last-plot default; one
supplemental was written twice to the same filename with different content;
`showtext` DPI 400 vs `ggsave`'s 300 default. Full list in the spec (21 items).

## Structure

- `common/manifest.R` — shared loaders; metadata attach generalised so it is
  figure-agnostic. Fig1 migrated (two token changes; its loaded frames are
  unchanged).
- `Fig3/00_setup.R` — paths, theme, FastENLOC readers, FDR labelling, and six
  memoized derived-data builders cached to `fig3/data/`, so panels declare only
  what they need without re-running the shared upstream three times.
- `Fig3/panels/panelA.R … panelE.R`, `assemble.R`, `run_all.R`, two supplementals.
- `tests/` — plain Rscript checks against mock inputs; no controlled-access data
  required. `Rscript tests/run_all.R`.

## Not verified here

The real inputs are *All of Us* controlled-access and unavailable outside the
Researcher Workbench, so validation is parse checks, unit tests and a full
end-to-end build against mock fixtures. The real test is `Rscript fig3/run_all.R`
from the project root on Terra.

Panel E is a reserved, undecided slot; `NOTES.md` records its layout mechanics.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_013yRNqXjLZrciXpcd2btdSP
BODY
)"
```

---

## Traceability — spec section → task

| Spec § | Requirement | Task |
|---|---|---|
| §3 | Target directory structure | 2, 8, 9, 10, 11 |
| §4 | Extract `common/manifest.R` | 1 |
| §4.1 | Generalised metadata attach; Fig1's two token changes | 1 |
| §5 | `build_or_cache` + six builders | 4, 5, 6, 7 |
| §6 | Panel decomposition and layout | 8, 9, 10 |
| §7.1.1 | Escape-sequence regex | 6 |
| §7.1.2 | Duplicated `bind_rows` argument | 6 |
| §7.1.3 | Hardcoded 86.2% | 8 |
| §7.1.4 | Duplicate model fits | 9 |
| §7.1.5 | Dead `rare` assignment | 7 |
| §7.1.6 | `QTL_union_gene_list` → `genes_all_three_assays` | 4 |
| §7.2.7 | Save the assembled figure | 10 |
| §7.2.8 | Explicit `plot=` in `ggsave` | 10 |
| §7.2.9 | Supplemental written twice | 10 |
| §7.2.10 | Matched `showtext`/`ggsave` DPI | 2, 10 |
| §7.2.11 | Remove `options(repr.plot.*)` | 2 (never carried over) |
| §7.2.12 | Drop figure uploads; document in README | 10, 11 |
| §7.3.13 | Delete `LoadEnlocGeneLevel` | 3 |
| §7.3.14 | Unify the two enloc readers | 3 |
| §7.3.15 | Explicit return | 3 |
| §7.3.16 | `components(g)` computed once | 6 |
| §7.3.17 | `bind_cols` → `mutate` | 7 |
| §7.3.18 | Report the FDR threshold | 3 |
| §7.3.19 | No-ops and redundant filters | 4, 6, 9 |
| §7.3.20 | Typos | 8, 9, 10 |
| §7.3.21 | Double tagging | 8, 9, 10 |
| §8 | Deferred analyses, panel E | 9, 11 |
| §9 | Input manifest | 2 |
