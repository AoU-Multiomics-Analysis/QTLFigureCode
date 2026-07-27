# Supplemental Figure 1 Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the flat 316-line `SuppFig1/SuppFig1.R` (AoU-vs-GTEx expression comparison) with a modular, manifest-driven pipeline matching Fig3/Fig4, extend the shared loader to stage `https://` inputs, fix the TF double-log, and generate the figure's notebook.

**Architecture:** `SuppFig1/00_setup.R` sources the shared `common/manifest.R` (now https-capable), then adds paths, theme, folded helpers, and memoized builders cached to `suppfig1/data/`. Four panel scripts build the main figure; four `extra/` scripts reproduce the retained per-gene-type comparisons; `assemble.R` composes the panels; `run_all.R` runs everything.

**Tech Stack:** R 4.3, tidyverse, data.table, rtracklayer, enrichR, hexbin, ggrepel, rcartocolor, broom, patchwork, showtext/sysfonts, paletteer.

**Spec:** `docs/superpowers/specs/2026-07-24-suppfig1-cleanup-design.md`

## Global Constraints

- **Branch:** `suppfig1-cleanup` (off `fig4-cleanup`; spec committed as `7ca5ab2`).
- **Run every script from the project root**, referring to it as lowercase `suppfig1/...` (macOS is case-insensitive; matches the other figures).
- **Never commit inputs or figure outputs** — `aou_tpms` and `go_annotations` are *All of Us* controlled-access.
- **`showtext` DPI must equal every `ggsave()` DPI.** SuppFig1 uses `DPI <- 300`.
- **`element_text()` sizes are pt; `geom_*`/`annotate()` sizes are mm.**
- **Tagging is `assemble.R`'s job.** No panel script may call `labs(tag = ...)`.
- **The derived-data cache is keyed by name only.**
- Commit after every task. Never `git add -A` (controlled-access files may be untracked).

---

## File Structure

| File | Responsibility |
|---|---|
| `common/manifest.R` | **Modify.** Extend the staging loop to handle `https`/`http`/`file` URIs, not only `gs://`. |
| `tests/test_manifest.R` | **Modify.** Add an https-branch (via `file://`) staging test. |
| `SuppFig1/00_setup.R` | **Create.** Libs, `SUPPFIG1_*` paths, theme, folded helpers, `build_or_cache`, four builders. |
| `SuppFig1/inputs.tsv` | **Create.** 5 keys. |
| `SuppFig1/panels/panelA.R`–`panelD.R` | **Create.** |
| `SuppFig1/extra/*.R` | **Create.** Four retained per-gene-type scripts. |
| `SuppFig1/assemble.R`, `run_all.R` | **Create.** |
| `SuppFig1/README.md`, `NOTES.md`, `.gitignore` | **Create.** |
| `SuppFig1/SuppFig1.R` | **Commit unchanged (Task 2), then `git rm` (Task 10).** |
| `.gitignore` | **Modify.** Add `/expression/`. |
| `tools/build_notebooks.py` | **Modify.** Add `SuppFig1` config (Task 9). |
| `tests/fixtures.R` | **Modify.** Add `make_suppfig1_fixture()`. |
| `tests/test_suppfig1_helpers.R` | **Create.** |
| `tests/test_suppfig1_assembly.R` | **Create.** |

---

## Task 1: Extend `common/manifest.R` to stage https/http/file URIs

**Files:**
- Modify: `common/manifest.R` (the `input_files` staging loop)
- Modify: `tests/test_manifest.R`

**Interfaces:**
- Produces: staging now dispatches on the `gs_uri` scheme — `gs://` via `gsutil cp`, `http(s)://` and `file://` via `download.file`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_manifest.R`, before its final `cat(...)`:

```r
# --- https/file staging branch (via a local file:// URI, no network) --------
src <- file.path(tmp, "sub", "remote_src.tsv")
wtsv(data.frame(v = 42), src)
writeLines(c(
  "key\ttype\tlocal_subpath\tgs_uri",
  paste0("remote\t-\tsub/fetched.tsv\tfile://", src)
), file.path(tmp, "inputs.tsv"))
.manifest_env$m <- NULL                     # force a re-read of the rewritten manifest
f <- input_files("remote")
ok("file:// URI stages a missing local input",
   file.exists(file.path(tmp, "sub/fetched.tsv")))
ok("staged file has the source contents",
   identical(readLines(f$path)[1], "v"))
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript tests/test_manifest.R`
Expected: FAIL — the `remote` key errors (`gsutil cp failed` / unsupported), because staging only knows `gs://`.

- [ ] **Step 3: Write the implementation**

In `common/manifest.R`, replace the two staging lines inside `input_files`:

```r
    message("  staging ", gs, " -> ", lp)
    status <- system2("gsutil", c("cp", gs, lp))
    if (status != 0L || !file.exists(lp)) stop("gsutil cp failed for ", gs)
```

with a call to a new helper:

```r
    message("  staging ", gs, " -> ", lp)
    .stage_uri(gs, lp)
```

and add the helper just above `input_files` (after `manifest_rows`):

```r
# Stage a single remote URI to a local destination, dispatching on scheme:
#   gs://            -> gsutil cp   (needs the Google Cloud SDK)
#   http(s):// file:// -> download.file
# Errors (naming the URI) if the scheme is unsupported or the fetch fails.
.stage_uri <- function(uri, dest) {
  if (grepl("^gs://", uri)) {
    status <- system2("gsutil", c("cp", uri, dest))
    if (status != 0L || !file.exists(dest)) stop("gsutil cp failed for ", uri)
  } else if (grepl("^(https?|file)://", uri)) {
    ok <- tryCatch({ utils::download.file(uri, dest, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (!isTRUE(ok) || !file.exists(dest)) stop("download.file failed for ", uri)
  } else {
    stop("Unsupported URI scheme for staging (key's gs_uri = '", uri, "').")
  }
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/test_manifest.R`
Expected: PASS — both new `ok - ...` lines, plus the existing assertions.

- [ ] **Step 5: Run the full suite (no regressions in the other figures' loaders)**

Run: `Rscript tests/run_all.R`
Expected: `All test files passed.`

- [ ] **Step 6: Commit**

```bash
git add common/manifest.R tests/test_manifest.R
git commit -m "Stage https/file URIs in the shared manifest loader

The loader only knew gsutil cp; dispatch on the gs_uri scheme so http(s):// and
file:// inputs stage via download.file. Lets public inputs (GENCODE, GTEx, the TF
list) self-stage. Tested with a local file:// URI so no network is needed."
```

---

## Task 2: Preserve original, SuppFig1 setup skeleton, manifest, gitignore

**Files:**
- Commit then later remove: `SuppFig1/SuppFig1.R`
- Create: `SuppFig1/00_setup.R` (paths/libs/theme/`build_or_cache` only)
- Create: `SuppFig1/inputs.tsv`, `SuppFig1/.gitignore`
- Modify: `.gitignore`

**Interfaces:**
- Consumes: `common/manifest.R`.
- Produces: globals `INPUT`, `DATA_DIR`, `OBJ_DIR`, `OUTDIR`, `MANIFEST`, `FONT`, `DPI`, `SZ_*`; `build_or_cache(name, expr, refresh = FALSE)`.

- [ ] **Step 1: Preserve the original in history**

```bash
git add SuppFig1/SuppFig1.R
git commit -m "Add the original SuppFig1.R before cleanup

SuppFig1/SuppFig1.R was never committed. Recording it unchanged so the original
analysis stays recoverable, and so the restructuring is reviewable as a diff."
```

- [ ] **Step 2: Create the manifest**

Create `SuppFig1/inputs.tsv` — **tab-separated**, 1 header + 5 rows:

```
key	local_subpath	gs_uri
aou_tpms	expression/aou_rnaseq_121624.gene_tpm.fixed.gct.gz	gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/QTL_5500_subset/expression/aou_rnaseq_121624.gene_tpm.fixed.gct.gz
gtex_tpms	external/gene_tpm_v10_whole_blood.gct.gz	https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_whole_blood.gct.gz
gencode_gtf	reference/gencode.v45.annotation.gtf.gz	https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
tf_list	reference/TF_names_v_1.01.txt	https://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt
go_annotations	reference/GO_annotations.csv	gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/annotations/GO_annotations.csv
```

- [ ] **Step 3: Create the setup skeleton**

Create `SuppFig1/00_setup.R`:

```r
# =============================================================================
# Supplemental Figure 1 — shared setup (sourced by every panel / extra / assembly
# script). Run all scripts from the PROJECT ROOT:
#   source("suppfig1/00_setup.R")
#
# SuppFig1 compares All of Us bulk RNA-seq expression to GTEx whole blood across
# gene types (coding correlation + detection, GO enrichment, TF expression,
# lncRNA), with additional per-gene-type comparisons in suppfig1/extra/.
#
# See docs/superpowers/specs/2026-07-24-suppfig1-cleanup-design.md
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(rtracklayer)
  library(enrichR)
  library(hexbin)
  library(ggrepel)
  library(rcartocolor)
  library(broom)
  library(patchwork)
  library(showtext)
  library(sysfonts)
  library(paletteer)
})

# --- Paths (override via env vars if needed) --------------------------------
INPUT    <- Sys.getenv("SUPPFIG1_INPUT",    ".")                 # staging root (cwd)
DATA_DIR <- Sys.getenv("SUPPFIG1_DATA",     "suppfig1/data")     # memoized tables
OBJ_DIR  <- Sys.getenv("SUPPFIG1_OBJ",      "suppfig1/objects")  # per-panel objects
OUTDIR   <- Sys.getenv("SUPPFIG1_OUTDIR",   ".")                 # final outputs
MANIFEST <- Sys.getenv("SUPPFIG1_MANIFEST", "suppfig1/inputs.tsv")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OBJ_DIR,  showWarnings = FALSE, recursive = TRUE)

source("common/manifest.R")

# --- Font: Inter via Google if reachable, else default sans -----------------
FONT <- "sans"
if (!"inter" %in% sysfonts::font_families()) {
  if (isTRUE(tryCatch({ font_add_google("Inter", "inter"); TRUE }, error = function(e) FALSE)))
    FONT <- "inter"
} else FONT <- "inter"

DPI <- 300                     # showtext DPI; must match every ggsave() DPI
showtext_opts(dpi = DPI)
showtext_auto()

# --- Design system ----------------------------------------------------------
# element_text() sizes are in pt; geom_text()/annotate() are in mm. This figure
# uses bordered panels (folded in from the original's theme_update).
SZ_AXIS_TITLE <- 20; SZ_AXIS_TEXT <- 15; SZ_STRIP <- 20
SZ_LEG_TEXT   <- 12; SZ_LEG_TITLE <- 14; SZ_TAG   <- 25

theme_set(theme_classic())
theme_update(
  text             = element_text(family = FONT),
  panel.border     = element_rect(fill = NA, linewidth = 0.9),
  axis.line        = element_blank(),
  axis.title       = element_text(size = SZ_AXIS_TITLE),
  axis.text        = element_text(size = SZ_AXIS_TEXT),
  strip.text       = element_text(size = SZ_STRIP),
  strip.background = element_blank(),
  legend.text      = element_text(size = SZ_LEG_TEXT),
  legend.title     = element_text(size = SZ_LEG_TITLE)
)

# =============================================================================
# Derived-data cache (identical contract to Fig3/Fig4): lazy `expr`, cache to
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

- [ ] **Step 4: Create SuppFig1/.gitignore**

```
# Per-panel ggplot object cache (regenerated by the pipeline)
objects/

# Memoized derived tables (regenerated by the builders in 00_setup.R)
data/

# macOS
.DS_Store
```

- [ ] **Step 5: Extend the root .gitignore**

Add `/expression/` to the controlled-access inputs block of `.gitignore`
(after `/external/`). `/external/` and `/reference/` are already present. Add the
outputs to the "Figure outputs" block:

```
/SupplementalFig1.pdf
/SupplementalFig1.png
```

- [ ] **Step 6: Verify setup sources and the manifest resolves**

Run:
```bash
Rscript -e '
  source("suppfig1/00_setup.R")
  stopifnot(nrow(input_manifest()) == 5L)
  stopifnot(length(unique(input_manifest()$key)) == 5L)
  stopifnot(is.function(build_or_cache), DPI == 300)
  cat("SuppFig1 setup OK: 5 rows, 5 keys\n")'
```
Expected: `SuppFig1 setup OK: 5 rows, 5 keys`

- [ ] **Step 7: Verify gitignore coverage**

Run:
```bash
git check-ignore -v expression/x.gct.gz external/x.gz reference/x.txt \
  SupplementalFig1.pdf SuppFig1/objects/A.rds SuppFig1/data/merged_tpms.rds
git check-ignore Fig2/objects/A.rds && echo "REGRESSION" || echo "Fig2 objects still tracked - OK"
```
Expected: a match per path; `Fig2 objects still tracked - OK`

- [ ] **Step 8: Commit**

```bash
git add SuppFig1/00_setup.R SuppFig1/inputs.tsv SuppFig1/.gitignore .gitignore
git commit -m "Add SuppFig1 setup skeleton and input manifest

Paths, font/DPI, bordered-panel theme and the derived-data cache; the manifest
declares 5 inputs. Only aou_tpms and go_annotations are controlled-access; the
GTEx TPMs, GENCODE GTF and TF list are public https inputs that self-stage."
```

---

## Task 3: Folded helpers and the `gene_annotations` / `merged_tpms` builders

**Files:**
- Modify: `SuppFig1/00_setup.R` (append)
- Modify: `tests/fixtures.R`
- Create: `tests/test_suppfig1_helpers.R`

**Interfaces:**
- Consumes: `build_or_cache`, `input_path`.
- Produces: `compute_mean_TPM(TPM_df, gene_annotations)`, `plot_gene_comparison(TPM_df, point_size = 3, detection_threshold = 1, gene_type = "", hex = FALSE)`, `join_mean_tpms(aou_means, gtex_means)` (pure — takes two `compute_mean_TPM` outputs), `gene_annotations(refresh = FALSE)`, `merged_tpms(refresh = FALSE)`. `merged_tpms()` / `join_mean_tpms()` return a tibble with `Description`, `gene_id`, `gene_type`, `mean_TPM_AoU`, `mean_TPM_GTEx`.

- [ ] **Step 1: Add the fixture**

Append to `tests/fixtures.R`:

```r
# --- SuppFig1 mock project ---------------------------------------------------
# Mock AoU / GTEx GCT TPM matrices (Name, Description, then sample columns) and a
# small GENCODE gene-type table, enough to exercise the pure transforms. The GTF
# import and enrichR calls are not run here (rtracklayer/enrichR/network).
make_suppfig1_fixture <- function(root = file.path(tempdir(), "suppfig1_fixture")) {
  unlink(root, recursive = TRUE)
  for (d in c("expression", "external", "reference"))
    dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
  wtsv <- function(df, p) write.table(df, p, sep = "\t", row.names = FALSE, quote = FALSE)

  # Six genes spanning the gene types the panels filter on.
  ids  <- c("ENSG001.1","ENSG002.1","ENSG003.1","ENSG004.1","ENSG005.1","ENSG006.1")
  desc <- c("GENEA","TFGENE","LNC1","SNR1","HBB","RRNA1")
  types <- c("protein_coding","protein_coding","lncRNA","snRNA","protein_coding","rRNA")

  aou <- data.frame(Name = ids, Description = desc,
                    s1 = c(10, 8, 0.5, 0.3, 40, 0.05),
                    s2 = c(12, 9, 0.4, 0.2, 44, 0.05))
  wtsv(aou, file.path(root, "expression/aou.gct.tsv"))   # plain tsv stand-in

  gtex <- data.frame(Name = ids, Description = desc,
                     s1 = c(11, 2, 0.02, 0.25, 38, 0.06),
                     s2 = c(9,  3, 0.03, 0.30, 42, 0.05))
  wtsv(gtex, file.path(root, "external/gtex.gct.tsv"))

  ann <- data.frame(gene_id = ids, gene_type = types)
  saveRDS(ann, file.path(root, "reference/gene_annotations.rds"))

  list(root = root, aou_tpms = file.path(root, "expression/aou.gct.tsv"),
       gtex_tpms = file.path(root, "external/gtex.gct.tsv"),
       annotations = ann)
}
```

- [ ] **Step 2: Write the failing test**

Create `tests/test_suppfig1_helpers.R`:

```r
# Tests for SuppFig1's pure helper transforms. Run from the project root:
#   Rscript tests/test_suppfig1_helpers.R
source("tests/fixtures.R")
fx <- make_suppfig1_fixture()
source("suppfig1/00_setup.R")
source("tests/helpers.R")

# --- compute_mean_TPM ------------------------------------------------------
aou <- as_tibble(data.table::fread(fx$aou_tpms))
means <- compute_mean_TPM(aou, fx$annotations)
ok("compute_mean_TPM returns one row per gene", nrow(means) == 6L)
ok("mean_TPM is the row mean across samples",
   abs(means$mean_TPM[means$Description == "GENEA"] - 11) < 1e-9)
ok("gene_type is joined on", "gene_type" %in% names(means))
ok("no SE column is emitted", !"SE" %in% names(means))

# --- join_mean_tpms (pure) -------------------------------------------------
gtex <- as_tibble(data.table::fread(fx$gtex_tpms))
aou_means  <- compute_mean_TPM(aou,  fx$annotations)
gtex_means <- compute_mean_TPM(gtex, fx$annotations)
merged <- join_mean_tpms(aou_means, gtex_means)
ok("join_mean_tpms inner-joins on version-stripped gene id",
   all(c("mean_TPM_AoU", "mean_TPM_GTEx", "gene_id", "gene_type") %in% names(merged)))
ok("join keeps the shared genes", nrow(merged) == 6L)
ok("GENEA carries both cohorts' means",
   merged$mean_TPM_AoU[merged$Description == "GENEA"] == 11 &&
   merged$mean_TPM_GTEx[merged$Description == "GENEA"] == 10)

# --- plot_gene_comparison returns a patchwork of two plots -----------------
p <- plot_gene_comparison(merged, gene_type = "coding", hex = FALSE)
ok("plot_gene_comparison returns a 2-plot patchwork",
   inherits(p, "patchwork") && length(p$patches$plots) >= 1)

cat("test_suppfig1_helpers.R (Task 3): all assertions passed\n")
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `Rscript tests/test_suppfig1_helpers.R`
Expected: FAIL — `could not find function "compute_mean_TPM"`

- [ ] **Step 4: Write the implementation**

Append to `SuppFig1/00_setup.R`:

```r
# =============================================================================
# Folded helpers
# =============================================================================
# Mean TPM per gene from a GCT-style matrix (Name, Description, then samples),
# joined to gene-type annotation. The original also emitted SE = sd/n, which is
# neither a standard error (sd/sqrt(n)) nor plotted anywhere; it is dropped.
compute_mean_TPM <- function(TPM_df, gene_annotations) {
  annotation_dat <- gene_annotations %>% distinct()
  TPM_df %>%
    pivot_longer(!c(Name, Description)) %>%
    group_by(Description, Name) %>%
    summarize(mean_TPM = mean(value), .groups = "drop") %>%
    left_join(annotation_dat, by = c("Name" = "gene_id"))
}

# Correlation scatter (AoU vs GTEx log2 TPM) beside cohort-specific detection
# counts. `hex = TRUE` uses geom_hex; after_stat() replaces the deprecated stat().
plot_gene_comparison <- function(TPM_df, point_size = 3, detection_threshold = 1,
                                 gene_type = "", hex = FALSE) {
  axis_label <- paste0("Unique ", gene_type, " genes")
  cor_values <- cor(log2(TPM_df$mean_TPM_AoU + 1), log2(TPM_df$mean_TPM_GTEx + 1),
                    use = "complete.obs")
  plot_text <- paste0("Cor = ", round(cor_values, digits = 2))

  cor_plot <- TPM_df %>%
    ggplot(aes(x = log2(mean_TPM_AoU + 1), y = log2(mean_TPM_GTEx + 1))) +
    xlab("log2(AoU mean TPM)") + ylab("log2(GTEx mean TPM)") +
    annotate("text", x = 1.8, y = 19, label = plot_text, size = 7)
  if (!hex) {
    cor_plot <- cor_plot +
      geom_point(shape = 21, size = point_size, fill = "black", alpha = 0.7) +
      geom_abline(slope = 1, color = "red")
  } else {
    cor_plot <- cor_plot +
      geom_hex(aes(fill = after_stat(log10(count)))) +
      scale_fill_viridis_c() +
      geom_abline(slope = 1, color = "red")
  }

  detection_plot <- TPM_df %>%
    filter((mean_TPM_AoU > detection_threshold & mean_TPM_GTEx < detection_threshold) |
           (mean_TPM_AoU < detection_threshold & mean_TPM_GTEx > detection_threshold)) %>%
    mutate(type = case_when(mean_TPM_AoU > mean_TPM_GTEx ~ "AoU", TRUE ~ "GTEx")) %>%
    dplyr::count(type) %>%
    ggplot(aes(x = type, y = n, fill = type)) +
    geom_col(show.legend = FALSE) +
    ylab(axis_label) + xlab("Cohort") +
    scale_fill_carto_d(palette = "Safe")

  cor_plot + detection_plot
}

# =============================================================================
# Builders
# =============================================================================
gene_annotations <- function(refresh = FALSE)
  build_or_cache("gene_annotations", refresh = refresh, {
    rtracklayer::import(input_path("gencode_gtf")) %>%
      data.frame() %>% select(gene_id, gene_type) %>% distinct()
  })

# Pure join of two compute_mean_TPM() outputs: version-strip the gene id,
# inner-join AoU to GTEx, set missing means to 0. Factored out so it is testable
# without the GTF import / GCT reads that merged_tpms() performs.
join_mean_tpms <- function(aou_means, gtex_means) {
  strip <- function(d) d %>% distinct() %>% mutate(gene_id = str_remove(Name, "\\..*"))
  strip(aou_means) %>%
    select(-Name) %>%
    inner_join(strip(gtex_means) %>% select(-Description, -gene_type, -Name),
               by = "gene_id", suffix = c("_AoU", "_GTEx")) %>%
    mutate(mean_TPM_GTEx = coalesce(mean_TPM_GTEx, 0),
           mean_TPM_AoU  = coalesce(mean_TPM_AoU, 0))
}

# AoU + GTEx per-gene mean TPM, joined and gene-type-attached. The GCT files have
# a 2-line header, hence skip = 2.
merged_tpms <- function(refresh = FALSE)
  build_or_cache("merged_tpms", refresh = refresh, {
    ann  <- gene_annotations()
    aou  <- as_tibble(data.table::fread(input_path("aou_tpms"),  skip = 2)) %>% compute_mean_TPM(ann)
    gtex <- as_tibble(data.table::fread(input_path("gtex_tpms"), skip = 2)) %>% compute_mean_TPM(ann)
    join_mean_tpms(aou, gtex)
  })
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `Rscript tests/test_suppfig1_helpers.R`
Expected: PASS — 5 `ok - ...` lines then the completion message.

- [ ] **Step 6: Commit**

```bash
git add SuppFig1/00_setup.R tests/fixtures.R tests/test_suppfig1_helpers.R
git commit -m "Add SuppFig1 TPM helpers and gene_annotations/merged_tpms builders

Folds compute_mean_TPM (dropping the unused mislabeled SE column) and
plot_gene_comparison (after_stat replaces the deprecated stat()) into 00_setup.R,
and adds the gene_annotations + merged_tpms builders. tests/fixtures.R gains a
mock AoU/GTEx GCT project."
```

---

## Task 4: The `tf_list` and `go_enrichment` builders

**Files:**
- Modify: `SuppFig1/00_setup.R` (append)
- Modify: `tests/test_suppfig1_helpers.R` (append)

**Interfaces:**
- Consumes: `build_or_cache`, `input_path`, `merged_tpms`.
- Produces: `tf_list(refresh = FALSE)` (character vector of TF symbols); `cohort_specific_genes(merged, cohort)` (pure — the protein-coding genes detected in one cohort only); `go_enrichment(refresh = FALSE)` (a tibble of enrichR results with a `type` ∈ {AoU, GTEx}); `GO_DBS`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_suppfig1_helpers.R`, before the final `cat(...)`:

```r
# --- cohort_specific_genes (pure) ------------------------------------------
spec_mock <- tibble::tibble(
  Description = c("A", "B", "C", "D"),
  gene_type = "protein_coding",
  mean_TPM_AoU  = c(5,   0.2, 5,   0.2),
  mean_TPM_GTEx = c(0.2, 5,   5,   0.2))
ok("AoU-specific = coding, AoU>1 & GTEx<1",
   identical(cohort_specific_genes(spec_mock, "AoU"), "A"))
ok("GTEx-specific = coding, AoU<1 & GTEx>1",
   identical(cohort_specific_genes(spec_mock, "GTEx"), "B"))
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript tests/test_suppfig1_helpers.R`
Expected: FAIL — `could not find function "cohort_specific_genes"`

- [ ] **Step 3: Write the implementation**

Append to `SuppFig1/00_setup.R`:

```r
tf_list <- function(refresh = FALSE)
  build_or_cache("tf_list", refresh = refresh, {
    data.table::fread(input_path("tf_list"), header = FALSE)$V1
  })

GO_DBS <- c("GO_Molecular_Function_2025", "GO_Cellular_Component_2025",
            "GO_Biological_Process_2025")

# Protein-coding genes detected (mean TPM > 1) in exactly one cohort.
cohort_specific_genes <- function(merged, cohort) {
  d <- merged %>% filter(gene_type == "protein_coding")
  if (cohort == "AoU")
    d %>% filter(mean_TPM_AoU > 1,  mean_TPM_GTEx < 1) %>% pull(Description)
  else
    d %>% filter(mean_TPM_AoU < 1,  mean_TPM_GTEx > 1) %>% pull(Description)
}

# enrichR GO enrichment for AoU- vs GTEx-specific genes, row-bound with a `type`
# label. Cached so the enrichR HTTP call runs once.
go_enrichment <- function(refresh = FALSE)
  build_or_cache("go_enrichment", refresh = refresh, {
    merged <- merged_tpms()
    aou  <- enrichR::enrichr(cohort_specific_genes(merged, "AoU"),  GO_DBS) %>% bind_rows()
    gtex <- enrichR::enrichr(cohort_specific_genes(merged, "GTEx"), GO_DBS) %>% bind_rows()
    bind_rows(aou %>% mutate(type = "AoU"), gtex %>% mutate(type = "GTEx"))
  })
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/test_suppfig1_helpers.R`
Expected: PASS — the two new `cohort_specific_genes` assertions plus the earlier ones.

- [ ] **Step 5: Commit**

```bash
git add SuppFig1/00_setup.R tests/test_suppfig1_helpers.R
git commit -m "Add SuppFig1 tf_list and go_enrichment builders

tf_list reads the human TF names; go_enrichment computes cohort-specific genes
(factored into a pure cohort_specific_genes, tested directly) and caches the
enrichR GO call. GO_DBS lists the three ontologies queried."
```

---

## Task 5: Panels A and B

**Files:**
- Create: `SuppFig1/panels/panelA.R`, `SuppFig1/panels/panelB.R`
- Modify: `tests/test_suppfig1_helpers.R` (append — Panel B's plot-prep)

**Interfaces:**
- Consumes: `merged_tpms`, `plot_gene_comparison`, `go_enrichment`, `format_go_plot_data`, `OBJ_DIR`.
- Produces: `SuppFig1/objects/A.rds`, `B.rds`; `format_go_plot_data(go_df, top_n = 10)`.

- [ ] **Step 1: Write the failing test**

Append to `tests/test_suppfig1_helpers.R`, before the final `cat(...)`:

```r
# --- format_go_plot_data (Panel B prep, pure) ------------------------------
go_mock <- tibble::tibble(
  Term = c(paste0("T", 1:12, " (GO:000", 1:12, ")")),
  Adjusted.P.value = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08),
  type = rep(c("AoU", "GTEx"), 6))
gp <- format_go_plot_data(go_mock, top_n = 3)
ok("format_go_plot_data keeps <= top_n per cohort",
   all(gp %>% count(type) %>% pull(n) <= 3))
ok("the GO id suffix is stripped from Term", !any(grepl("\\(GO", gp$Term)))
ok("neglog10_p is -log10(adjusted p)",
   abs(gp$neglog10_p[1] - -log10(gp$Adjusted.P.value[1])) < 1e-9)
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript tests/test_suppfig1_helpers.R`
Expected: FAIL — `could not find function "format_go_plot_data"`

- [ ] **Step 3: Implement the Panel B prep helper**

Append to `SuppFig1/00_setup.R`:

```r
# Panel B prep: keep the top_n GO terms (by adjusted p) per cohort, strip the
# "(GO:...)" id from the term, and add neglog10_p.
format_go_plot_data <- function(go_df, top_n = 10) {
  go_df %>%
    group_by(type) %>%
    arrange(Adjusted.P.value, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(Term = str_remove(Term, "\\(GO.*"),
           neglog10_p = -log10(Adjusted.P.value))
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/test_suppfig1_helpers.R`
Expected: PASS — the three `format_go_plot_data` assertions plus the earlier ones.

- [ ] **Step 5: Write panel A**

Create `SuppFig1/panels/panelA.R`:

```r
# Panel A — protein-coding AoU-vs-GTEx TPM correlation (hex) + detection counts
# Inputs: aou_tpms, gtex_tpms, gencode_gtf
source("suppfig1/00_setup.R")

protein_coding_genes <- merged_tpms() %>%
  filter(gene_type == "protein_coding") %>%
  filter(mean_TPM_AoU > 1 | mean_TPM_GTEx > 1)

protein_coding_comparison <- protein_coding_genes %>%
  plot_gene_comparison(gene_type = "coding", hex = TRUE)

saveRDS(protein_coding_comparison, file.path(OBJ_DIR, "A.rds"))
message("Panel A built.")
```

- [ ] **Step 6: Write panel B**

Create `SuppFig1/panels/panelB.R`:

```r
# Panel B — GO enrichment of cohort-specific genes (AoU vs GTEx)
# Inputs: aou_tpms, gtex_tpms, gencode_gtf (via go_enrichment)
source("suppfig1/00_setup.R")

AoU_GTEx_enrichment_plot <- go_enrichment() %>%
  format_go_plot_data(top_n = 10) %>%
  ggplot(aes(x = neglog10_p, y = reorder(Term, neglog10_p))) +
  geom_col() +
  facet_wrap(~type, scales = "free_y") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
  xlab("-log10(FDR)") + ylab("") +
  geom_vline(xintercept = 1.3) +
  theme(panel.spacing = unit(1, "cm"), strip.text = element_text(size = 20))

saveRDS(AoU_GTEx_enrichment_plot, file.path(OBJ_DIR, "B.rds"))
message("Panel B built.")
```

- [ ] **Step 7: Verify no panel sets its own tag**

Run: `grep -n "labs(tag" SuppFig1/panels/*.R; test $? -eq 1 && echo "OK: no panel sets a tag"`
Expected: `OK: no panel sets a tag`

- [ ] **Step 8: Commit**

```bash
git add SuppFig1/panels/panelA.R SuppFig1/panels/panelB.R SuppFig1/00_setup.R tests/test_suppfig1_helpers.R
git commit -m "Add SuppFig1 panels A and B

Panel A (protein-coding AoU-vs-GTEx correlation + detection) via merged_tpms +
plot_gene_comparison; Panel B (GO enrichment) via go_enrichment + the pure
format_go_plot_data prep, tested directly."
```

---

## Task 6: Panels C and D

**Files:**
- Create: `SuppFig1/panels/panelC.R`, `SuppFig1/panels/panelD.R`
- Modify: `tests/test_suppfig1_helpers.R` (append — the TF single-log check)

**Interfaces:**
- Consumes: `merged_tpms`, `tf_list`, `plot_gene_comparison`, `tf_expression_table`, `OBJ_DIR`.
- Produces: `SuppFig1/objects/C.rds`, `D.rds`; `tf_expression_table(merged, tf_symbols)` returning the TF genes with raw `mean_TPM_AoU`/`mean_TPM_GTEx` (NOT pre-logged — the fix).

- [ ] **Step 1: Write the failing test**

Append to `tests/test_suppfig1_helpers.R`, before the final `cat(...)`:

```r
# --- tf_expression_table does NOT pre-log the TPMs (the double-log fix) -----
tf_mock <- tibble::tibble(
  Description = c("TF1", "TF2", "OTHER"),
  gene_type = "protein_coding",
  mean_TPM_AoU  = c(8, 16, 10),
  mean_TPM_GTEx = c(4,  4,  4))
tf_tab <- tf_expression_table(tf_mock, c("TF1", "TF2"))
ok("tf_expression_table keeps only TF-list genes above threshold",
   setequal(tf_tab$Description, c("TF1", "TF2")))
ok("TPMs are raw, not pre-logged (TF2 AoU stays 16, not log2(16+1))",
   tf_tab$mean_TPM_AoU[tf_tab$Description == "TF2"] == 16)

cat("test_suppfig1_helpers.R (Task 3): all assertions passed\n")
```

(Delete the old final `cat(...)` line that this block now precedes — keep only one.)

- [ ] **Step 2: Run the test to verify it fails**

Run: `Rscript tests/test_suppfig1_helpers.R`
Expected: FAIL — `could not find function "tf_expression_table"`

- [ ] **Step 3: Implement the TF table helper**

Append to `SuppFig1/00_setup.R`:

```r
# TF genes present in the TF list and detected (TPM > 1) in either cohort. TPMs
# are returned RAW; the plot / lm / t-test each apply log2(x+1) exactly once. The
# original pre-logged here and then logged again downstream (see NOTES.md).
tf_expression_table <- function(merged, tf_symbols) {
  merged %>%
    filter(gene_type == "protein_coding") %>%
    filter(Description %in% tf_symbols) %>%
    filter(mean_TPM_GTEx > 1 | mean_TPM_AoU > 1)
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/test_suppfig1_helpers.R`
Expected: PASS — the two TF assertions plus everything earlier; final completion line.

- [ ] **Step 5: Write panel C**

Create `SuppFig1/panels/panelC.R`:

```r
# Panel C — TF expression boxplot (AoU vs GTEx), log2(TPM+1) once
# Inputs: aou_tpms, gtex_tpms, gencode_gtf, tf_list
source("suppfig1/00_setup.R")

TF_genes <- tf_expression_table(merged_tpms(), tf_list())

# Reported statistics (log exactly once, matching the axis label).
tf_lm <- TF_genes %>%
  do(broom::tidy(lm(log2(mean_TPM_GTEx + 1) ~ log2(mean_TPM_AoU + 1), data = .)))
tf_ttest <- t.test(log2(TF_genes$mean_TPM_AoU + 1),
                   log2(TF_genes$mean_TPM_GTEx + 1), paired = TRUE)
message("Panel C TF lm:"); print(tf_lm)
message("Panel C TF paired t-test p = ", signif(tf_ttest$p.value, 3))

TF_gene_comparison <- TF_genes %>%
  select(Description, contains("TPM")) %>%
  pivot_longer(!Description) %>%
  mutate(name = str_remove(name, "mean_TPM_")) %>%
  ggplot(aes(x = name, y = log2(value + 1), fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  ylab("TF Expression \n log2(TPM + 1)") +
  scale_fill_carto_d(palette = "Safe")

saveRDS(TF_gene_comparison, file.path(OBJ_DIR, "C.rds"))
message("Panel C built.")
```

- [ ] **Step 6: Write panel D**

Create `SuppFig1/panels/panelD.R`:

```r
# Panel D — lncRNA AoU-vs-GTEx TPM correlation (hex) + detection counts
# Inputs: aou_tpms, gtex_tpms, gencode_gtf
source("suppfig1/00_setup.R")

lncRNA_comparison <- merged_tpms() %>%
  filter(gene_type == "lncRNA") %>%
  filter(mean_TPM_AoU > 0.1 | mean_TPM_GTEx > 0.1) %>%
  plot_gene_comparison(gene_type = "lncRNA", hex = TRUE)

saveRDS(lncRNA_comparison, file.path(OBJ_DIR, "D.rds"))
message("Panel D built.")
```

- [ ] **Step 7: Verify no panel sets its own tag**

Run: `grep -n "labs(tag" SuppFig1/panels/*.R; test $? -eq 1 && echo "OK: no panel sets a tag"`
Expected: `OK: no panel sets a tag`

- [ ] **Step 8: Commit**

```bash
git add SuppFig1/panels/panelC.R SuppFig1/panels/panelD.R SuppFig1/00_setup.R tests/test_suppfig1_helpers.R
git commit -m "Add SuppFig1 panels C and D

Panel C (TF expression boxplot + lm/t-test) via tf_expression_table, which
returns RAW TPMs so the plot/lm/t-test log exactly once (fixing the original's
double log). Panel D (lncRNA correlation + detection) via plot_gene_comparison."
```

---

## Task 7: The retained `extra/` scripts

**Files:**
- Create: `SuppFig1/extra/gene_type_comparisons.R`, `rrna_expression.R`, `hemoglobin_genes.R`, `gene_type_detection.R`

**Interfaces:**
- Consumes: `merged_tpms`, `plot_gene_comparison`, `OUTDIR`.
- Produces: standalone PDFs in `OUTDIR` (not part of the main assembly).

- [ ] **Step 1: Write the four extra scripts**

Create `SuppFig1/extra/gene_type_comparisons.R`:

```r
# Extra — snRNA and miRNA AoU-vs-GTEx comparisons (not in the main figure).
# Inputs: aou_tpms, gtex_tpms, gencode_gtf
source("suppfig1/00_setup.R")

merged <- merged_tpms()

for (gt in c("snRNA", "miRNA")) {
  p <- merged %>%
    filter(gene_type == gt) %>%
    filter(mean_TPM_AoU > 0.1 | mean_TPM_GTEx > 0.1) %>%
    plot_gene_comparison(gene_type = gt, hex = FALSE)
  ggsave(file.path(OUTDIR, paste0("SuppFig1_extra_", gt, ".pdf")),
         p, width = 10, height = 5, dpi = DPI)
  message("wrote SuppFig1_extra_", gt, ".pdf")
}
```

Create `SuppFig1/extra/rrna_expression.R`:

```r
# Extra — rRNA mean-TPM histogram (AoU vs GTEx). Not in the main figure.
# Inputs: aou_tpms, gtex_tpms, gencode_gtf
source("suppfig1/00_setup.R")

rRNA_hist <- merged_tpms() %>%
  filter(gene_type == "rRNA") %>%
  select(Description, mean_TPM_AoU, mean_TPM_GTEx) %>%
  pivot_longer(!Description) %>%
  mutate(name = str_remove(name, "mean_TPM_")) %>%
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(position = "identity", alpha = 0.7) +
  xlab("Mean rRNA TPM") + ylab("Number genes") +
  scale_fill_carto_d(palette = "Safe")

ggsave(file.path(OUTDIR, "SuppFig1_extra_rRNA.pdf"), rRNA_hist,
       width = 7, height = 5, dpi = DPI)
message("wrote SuppFig1_extra_rRNA.pdf")
```

Create `SuppFig1/extra/hemoglobin_genes.R`:

```r
# Extra — hemoglobin (HB*) protein-coding genes, AoU-vs-GTEx scatter. Not in the
# main figure. Inputs: aou_tpms, gtex_tpms, gencode_gtf
source("suppfig1/00_setup.R")

hb_plot <- merged_tpms() %>%
  filter(str_detect(Description, "^HB") & gene_type == "protein_coding") %>%
  ggplot(aes(x = log2(mean_TPM_AoU), y = log2(mean_TPM_GTEx))) +
  geom_point(size = 6, shape = 21, fill = "black", alpha = 0.6) +
  geom_abline(slope = 1, color = "red") +
  ylim(-5, 20) + xlim(-5, 20) +
  xlab("log2(Mean AoU TPM)") + ylab("log2(Mean GTEx TPM)") +
  ggrepel::geom_text_repel(aes(label = Description))

ggsave(file.path(OUTDIR, "SuppFig1_extra_hemoglobin.pdf"), hb_plot,
       width = 7, height = 7, dpi = DPI)
message("wrote SuppFig1_extra_hemoglobin.pdf")
```

Create `SuppFig1/extra/gene_type_detection.R`:

```r
# Extra — number of genes detected (TPM > 0.1) per gene type in AoU. Not in the
# main figure. Inputs: aou_tpms, gencode_gtf
source("suppfig1/00_setup.R")

small_rnas_list <- c("misc_RNA", "sRNA", "vault_RNA", "scaRNA", "rRNA", "ribozyme")

# AoU per-gene means (one side of merged_tpms) with the gene type.
aou_means <- merged_tpms() %>% transmute(Description, gene_type, mean_TPM = mean_TPM_AoU)

gene_detection_summary <- aou_means %>%
  filter(mean_TPM > 0.1) %>%
  dplyr::count(gene_type) %>%
  mutate(gene_type = case_when(gene_type == "protein_coding" ~ "Coding", TRUE ~ gene_type)) %>%
  filter(!str_detect(gene_type, "IG|TR|pseudo|artifact|Mt|TEC|rRNA"))

other_small_rnas <- gene_detection_summary %>%
  filter(gene_type %in% small_rnas_list) %>%
  summarize(gene_type = "other small RNAs", n = sum(n))

detection_plot <- gene_detection_summary %>%
  bind_rows(other_small_rnas) %>%
  filter(!gene_type %in% small_rnas_list) %>%
  ggplot(aes(x = gene_type, y = n, fill = gene_type, color = gene_type)) +
  geom_col(show.legend = FALSE, alpha = 0.7) +
  ylab("Number detected \n (TPM > 0.1)") +
  scale_y_log10() +
  scale_fill_carto_d(palette = "Safe") + scale_color_carto_d(palette = "Safe") +
  geom_text(aes(label = n), vjust = -0.25, color = "black") +
  xlab("Type of gene")

ggsave(file.path(OUTDIR, "SuppFig1_extra_gene_detection.pdf"), detection_plot,
       width = 12, height = 7, dpi = DPI)
message("wrote SuppFig1_extra_gene_detection.pdf")
```

- [ ] **Step 2: Parse-check the four scripts**

Run: `for f in SuppFig1/extra/*.R; do Rscript -e "invisible(parse('$f'))" && echo "parse OK: $f"; done`
Expected: `parse OK` for each.

- [ ] **Step 3: Commit**

```bash
git add SuppFig1/extra/
git commit -m "Add SuppFig1 extra per-gene-type comparison scripts

Retained (but not in the assembled figure): snRNA/miRNA comparisons, the rRNA
mean-TPM histogram, the hemoglobin-gene scatter, and the per-gene-type detection
bar. Each reuses merged_tpms/plot_gene_comparison and saves its own PDF. The
broken isoform histogram (undefined input) is not carried over -- see NOTES.md."
```

---

## Task 8: Assembly, runner, and the rendered-tag test

**Files:**
- Create: `SuppFig1/assemble.R`, `SuppFig1/run_all.R`, `tests/test_suppfig1_assembly.R`

**Interfaces:**
- Consumes: `SuppFig1/objects/A.rds`–`D.rds`, `SZ_TAG`, `DPI`, `OUTDIR`.
- Produces: `SupplementalFig1.pdf`, `SupplementalFig1.png`.

- [ ] **Step 1: Write assemble.R**

Create `SuppFig1/assemble.R`:

```r
# =============================================================================
# Supplemental Figure 1 — assembly. Loads the per-panel objects
# (suppfig1/objects/*.rds) and composes the final figure.
#
# Panels A and D are each a 2-plot patchwork, so plot_annotation(tag_levels="A")
# letters SIX leaf plots (A,B on coding; C on GO; D on TF; E,F on lncRNA). That
# six-tag layout is intentional (see the design spec) and pinned by
# tests/test_suppfig1_assembly.R. Tagging happens HERE and only here.
# =============================================================================
source("suppfig1/00_setup.R")

P <- setNames(lapply(LETTERS[1:4],
                     function(x) readRDS(file.path(OBJ_DIR, paste0(x, ".rds")))),
              LETTERS[1:4])

top_row    <- P$A | (P$B + theme(plot.margin = margin(l = 2, unit = "cm")))
bottom_row <- P$C | P$D

SuppFig1 <- top_row / bottom_row +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = SZ_TAG, face = "bold"))

ggsave(file.path(OUTDIR, "SupplementalFig1.pdf"), SuppFig1,
       width = 35, height = 23, dpi = DPI, limitsize = FALSE)
showtext_opts(dpi = 100)
ggsave(file.path(OUTDIR, "SupplementalFig1.png"), SuppFig1,
       width = 35, height = 23, dpi = 100, limitsize = FALSE)
showtext_opts(dpi = DPI)

message("Assembled: ", file.path(OUTDIR, "SupplementalFig1.pdf"), " and .png")
```

- [ ] **Step 2: Write run_all.R**

Create `SuppFig1/run_all.R`:

```r
# =============================================================================
# Supplemental Figure 1 — end-to-end runner. From the project root:
#   Rscript suppfig1/run_all.R
#   1. panels  : each stages the inputs it needs (suppfig1/inputs.tsv) and builds
#                the derived tables it needs (memoized in suppfig1/data/)
#   2. assemble: compose + save SupplementalFig1.pdf / .png
# The extra/ scripts are standalone: Rscript suppfig1/extra/<name>.R
# The derived-data cache is keyed by NAME: clear suppfig1/data/ after changing an input.
# =============================================================================
for (p in LETTERS[1:4]) source(sprintf("suppfig1/panels/panel%s.R", p))
source("suppfig1/assemble.R")
message("run_all complete.")
```

- [ ] **Step 3: Write the assembly tag test**

Create `tests/test_suppfig1_assembly.R` (complete and self-contained — it builds
the four panel objects from the fixture and a stand-in Panel B, since Panel B's
real builder makes a live enrichR HTTP call the test host may lack; the contract
is **exactly six tags A–F**):

```r
# Asserts the assembled Supplemental Figure 1 carries exactly six tags A-F, by
# rendering and inspecting the tag grobs. Panels A and D are 2-plot patchworks,
# so six leaves are lettered -- intended (see the design spec). Panel B's real
# builder needs a live enrichR call, so a stand-in B object is used here; the
# real Panel B is exercised on the Workbench. Run from the project root.
source("tests/fixtures.R")
source("tests/helpers.R")
fx <- make_suppfig1_fixture()
Sys.setenv(SUPPFIG1_INPUT = fx$root,
           SUPPFIG1_DATA   = file.path(fx$root, "data"),
           SUPPFIG1_OBJ    = file.path(fx$root, "objects"),
           SUPPFIG1_OUTDIR = fx$root)
source("suppfig1/00_setup.R")

# Reuse the proven extractor from tests/test_fig3_assembly.R.
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

# Build merged TPMs from the fixture via the pure join (no GTF/GCT read needed).
ann    <- fx$annotations
means  <- function(p) compute_mean_TPM(as_tibble(data.table::fread(p)), ann)
merged <- join_mean_tpms(means(fx$aou_tpms), means(fx$gtex_tpms))

# Panels A, C, D are the real panel plots; B is a stand-in (real B needs enrichR).
saveRDS(merged %>% filter(gene_type == "protein_coding") %>%
          plot_gene_comparison(gene_type = "coding", hex = TRUE),
        file.path(OBJ_DIR, "A.rds"))
saveRDS(ggplot() + theme_void(base_family = FONT) +
          labs(title = "GO enrichment (stand-in)"),
        file.path(OBJ_DIR, "B.rds"))
saveRDS(tf_expression_table(merged, c("TFGENE")) %>%
          select(Description, contains("TPM")) %>% pivot_longer(!Description) %>%
          ggplot(aes(name, log2(value + 1))) + geom_boxplot(),
        file.path(OBJ_DIR, "C.rds"))
saveRDS(merged %>% filter(gene_type == "lncRNA") %>%
          plot_gene_comparison(gene_type = "lncRNA", hex = FALSE),
        file.path(OBJ_DIR, "D.rds"))

source("suppfig1/assemble.R")

tags <- extract_tag_labels(SuppFig1)
ok("Supplemental Figure 1 carries exactly 6 tags", length(tags) == 6L)
ok("its tags are exactly A, B, C, D, E, F",
   identical(sort(unique(tags)), c("A", "B", "C", "D", "E", "F")))
ok("the figure was written", file.exists(file.path(fx$root, "SupplementalFig1.pdf")))

cat("test_suppfig1_assembly.R: all assertions passed\n")
```

- [ ] **Step 4: Run the assembly test (this is the end-to-end build)**

Run: `Rscript tests/test_suppfig1_assembly.R`
Expected: `Assembled: …`, three `ok - ...` lines, `test_suppfig1_assembly.R: all assertions passed`.
(The real end-to-end run with the live Panel B happens on the Workbench.)

- [ ] **Step 5: Run the test suite**

Run: `Rscript tests/run_all.R`
Expected: all test files pass (assembly test may print its skip line if enrichR is absent).

- [ ] **Step 6: Commit**

```bash
git add SuppFig1/assemble.R SuppFig1/run_all.R tests/test_suppfig1_assembly.R
git commit -m "Add SuppFig1 assembly, runner and rendered-tag test

assemble.R composes the four panels into the 2x2 figure; because A and D are
2-plot patchworks the tag_levels='A' pass letters six leaves (A-F), which the
assembly test pins against the rendered grobs."
```

---

## Task 9: Notebook

**Files:**
- Modify: `tools/build_notebooks.py`
- Create: `SuppFig1/SuppFig1.ipynb` (generated)

**Interfaces:**
- Consumes: the SuppFig1 pipeline scripts.

- [ ] **Step 1: Add the SuppFig1 config to the generator**

In `tools/build_notebooks.py`, add a `"SuppFig1"` entry to the `FIGURES` dict
(after `"Fig4"`), and its panel notes to `PANEL_NOTES`:

```python
 "SuppFig1": {
   "num": 0, "cap": "SuppFig1", "lower": "suppfig1", "obj_letters": list("ABCD"),
   "title": "Supplemental Figure 1 — All of Us vs GTEx expression comparison",
   "overview": """
This notebook reproduces **Supplemental Figure 1**, comparing *All of Us* bulk
RNA-seq to GTEx whole blood across gene types: coding-gene TPM correlation and
cohort-specific detection, GO enrichment of cohort-specific genes, transcription-
factor expression, and lncRNA expression. Additional per-gene-type comparisons
are reproduced in the `extra/` sections at the end.

`00_setup.R` provides memoized builders (`gene_annotations()`, `merged_tpms()`,
`tf_list()`, `go_enrichment()`) cached to `SuppFig1/data/`.
""",
   "data_note": """
### Inputs

> **Controlled-access.** `aou_tpms` (and `go_annotations`) derive from *All of Us*
> controlled-access data. The GTEx TPMs, GENCODE GTF and TF list are public and
> self-stage from their `https://` URIs.

Inputs are manifest-driven (`SuppFig1/inputs.tsv`); builders cache to
`SuppFig1/data/` **by name**, so clear it after changing an input.

> **Results differ from the pre-cleanup script.** Panel C's TF statistics are
> logged once, not twice. See `SuppFig1/NOTES.md`.
""",
   "env": {"SUPPFIG1_MANIFEST": "inputs.tsv", "SUPPFIG1_OBJ": "objects",
           "SUPPFIG1_DATA": "data"},
   "steps": [
     ("panel","SuppFig1/panels/panelA.R","A","protein_coding_comparison"),
     ("panel","SuppFig1/panels/panelB.R","B","AoU_GTEx_enrichment_plot"),
     ("panel","SuppFig1/panels/panelC.R","C","TF_gene_comparison"),
     ("panel","SuppFig1/panels/panelD.R","D","lncRNA_comparison"),
     ("assemble","SuppFig1/assemble.R",None,"SuppFig1"),
     ("supp","SuppFig1/extra/gene_type_comparisons.R",None,None,
        "**Extra — snRNA / miRNA comparisons.** Retained per-gene-type AoU-vs-GTEx "
        "comparisons; not part of the assembled figure."),
     ("supp","SuppFig1/extra/rrna_expression.R",None,None,
        "**Extra — rRNA expression histogram.**"),
     ("supp","SuppFig1/extra/hemoglobin_genes.R",None,None,
        "**Extra — hemoglobin-gene scatter.**"),
     ("supp","SuppFig1/extra/gene_type_detection.R",None,None,
        "**Extra — genes detected per gene type.**"),
   ],
 },
```

Add to `PANEL_NOTES`:

```python
 ("SuppFig1","A"): "Protein-coding genes detected in either cohort: AoU-vs-GTEx log2(TPM+1) correlation (hex) beside counts of genes detected in only one cohort. Sub-plots are tagged A and B.",
 ("SuppFig1","B"): "Top GO terms (by adjusted p) enriched among genes specific to each cohort, faceted AoU vs GTEx, with a -log10(FDR)=1.3 reference line. Tagged C.",
 ("SuppFig1","C"): "Transcription-factor expression, AoU vs GTEx, as log2(TPM+1) boxplots, with an lm and a paired t-test reported in the run log. The TPMs are logged exactly once (the original logged twice). Tagged D.",
 ("SuppFig1","D"): "lncRNA AoU-vs-GTEx log2(TPM+1) correlation (hex) beside cohort-specific detection counts. Sub-plots are tagged E and F.",
```

- [ ] **Step 2: Generate the notebook**

Run: `python3 tools/build_notebooks.py SuppFig1`
Expected: `wrote SuppFig1/SuppFig1.ipynb  (N cells)`

- [ ] **Step 3: Validate the notebook**

Run:
```bash
python3 - <<'PY'
import json, subprocess, tempfile, os
nb=json.load(open("SuppFig1/SuppFig1.ipynb"))
assert nb["nbformat"]==4 and nb["metadata"]["kernelspec"]["name"]=="ir"
codes=[c for c in nb["cells"] if c["cell_type"]=="code"]
bad=[]
for i,c in enumerate(codes):
    with tempfile.NamedTemporaryFile("w",suffix=".R",delete=False) as f:
        f.write("".join(c["source"])); tmp=f.name
    r=subprocess.run(["Rscript","-e",f"invisible(parse('{tmp}'))"],capture_output=True,text=True)
    os.unlink(tmp)
    if r.returncode: bad.append((i,r.stderr.strip().splitlines()[-1]))
print(f"SuppFig1.ipynb: {len(nb['cells'])} cells | JSON OK | parse-fail: {len(bad)}")
for i,e in bad: print("  cell",i,e)
PY
```
Expected: `parse-fail: 0`.

- [ ] **Step 4: Commit**

```bash
git add tools/build_notebooks.py SuppFig1/SuppFig1.ipynb
git commit -m "Add SuppFig1 notebook and generator config

Extends tools/build_notebooks.py with SuppFig1 and commits SuppFig1/SuppFig1.ipynb
(the standard per-figure cleanup deliverable), covering setup, the four panels,
assembly, and the extra/ scripts."
```

---

## Task 10: Documentation, and retire SuppFig1.R

**Files:**
- Create: `SuppFig1/README.md`, `SuppFig1/NOTES.md`
- Remove: `SuppFig1/SuppFig1.R`

- [ ] **Step 1: Write NOTES.md**

Create `SuppFig1/NOTES.md`:

```markdown
# Supplemental Figure 1 — notes and deferred work

## Result-affecting change: the TF double log (Panel C)

`SuppFig1.R` applied `log2(x+1)` to the TF TPM columns in `TF_genes`
(`SuppFig1.R:173`), then applied `log2(·+1)` **again** in the boxplot
(`:184`), the `lm` (`:179`) and the paired `t.test` (`:192`) — so Panel C plotted
and tested `log2(log2(TPM+1)+1)` while the axis read `log2(TPM + 1)`. The cleanup
returns raw TPMs from `tf_expression_table()` and logs exactly once, so Panel C's
values and the reported lm / t-test change.

## Dropped: the isoform histogram

`SuppFig1.R:291–299` read `isoform_results`, a variable that is never defined (its
source path is commented out at `:22`), so the block errors. It is not carried
over. To reinstate it, add an isoform-results table to `inputs.tsv` and a small
`extra/` script that reads it and plots isoforms-per-gene.

## GO_annotations.csv is staged but not yet consumed

`go_annotations` is declared in `inputs.tsv` and can be read with
`fread(input_path("go_annotations"))`, but nothing in the current pipeline uses it
(the original only `head()`-ed it). It is kept as a hook; wire it into an analysis
when needed.

## `compute_mean_TPM` no longer emits `SE`

The original emitted `SE = sd(value)/n`, which is `sd/n` (not the standard error
`sd/sqrt(n)`) and was never plotted. The column is dropped.

## Retained extras (not in the assembled figure)

`extra/gene_type_comparisons.R` (snRNA, miRNA), `extra/rrna_expression.R`,
`extra/hemoglobin_genes.R`, and `extra/gene_type_detection.R` reproduce the
per-gene-type comparisons the original computed but did not assemble. Each is
standalone and saves its own `SuppFig1_extra_*.pdf`.
```

- [ ] **Step 2: Write README.md**

Create `SuppFig1/README.md`:

```markdown
# Supplemental Figure 1 — All of Us vs GTEx expression

Code to generate **Supplemental Figure 1**: a comparison of *All of Us* bulk
RNA-seq expression to GTEx whole blood across gene types.

> **Data note.** `aou_tpms` (and `go_annotations`) derive from *All of Us*
> **controlled-access** data and are not in this repository. The GTEx TPMs,
> GENCODE GTF and TF list are public and self-stage from `https://` URIs.

## The figure

`SupplementalFig1.pdf` — a 2×2 of composites, lettered A–F (panels A and D are
each a 2-plot pair):

| Tags | Content |
|------|---------|
| **A, B** | protein-coding AoU-vs-GTEx TPM correlation (hex) + cohort-specific detection counts |
| **C** | GO enrichment of cohort-specific genes (AoU vs GTEx) |
| **D** | transcription-factor expression (AoU vs GTEx), log2(TPM+1) |
| **E, F** | lncRNA AoU-vs-GTEx TPM correlation (hex) + detection counts |

Additional per-gene-type comparisons (snRNA, miRNA, rRNA, hemoglobin genes, and
per-type detection counts) are in `extra/`, each saving its own PDF.

## Running

From the **project root**:

    Rscript suppfig1/run_all.R              # panels -> assemble -> SupplementalFig1.pdf/.png
    Rscript suppfig1/extra/rrna_expression.R  # a standalone extra
    rm -rf suppfig1/data && Rscript suppfig1/run_all.R   # force-rebuild the cache
    Rscript tests/run_all.R                 # unit tests (no controlled-access data)

The panels share upstream through memoized builders in `00_setup.R`
(`gene_annotations`, `merged_tpms`, `tf_list`, `go_enrichment`), cached to
`suppfig1/data/` **by name** — clear it after changing an input. `go_enrichment`
makes a live enrichR HTTP call, cached so it runs once.

## Requirements

R ≥ 4.3, and CRAN/Bioconductor packages `tidyverse`, `data.table`, `rtracklayer`,
`enrichR`, `hexbin`, `ggrepel`, `rcartocolor`, `broom`, `patchwork`, `showtext`,
`sysfonts`, `paletteer`.

## Data & inputs

Every input is declared in `suppfig1/inputs.tsv`; the three public inputs
self-stage over https. See [`NOTES.md`](NOTES.md), in particular the
result-affecting TF double-log fix (Panel C).

## Provenance

Replaces a single 316-line notebook export (`SuppFig1.R`, recoverable from git
history). The restructuring and its fixes are recorded in
[`docs/superpowers/specs/2026-07-24-suppfig1-cleanup-design.md`](../docs/superpowers/specs/2026-07-24-suppfig1-cleanup-design.md).
```

- [ ] **Step 3: Delete the superseded script**

Run:
```bash
git rm SuppFig1/SuppFig1.R
test ! -e SuppFig1/SuppFig1.R && echo "SuppFig1.R removed"
```
Expected: `SuppFig1.R removed`

- [ ] **Step 4: Full suite + parse check + end-to-end**

Run:
```bash
Rscript tests/run_all.R
for f in SuppFig1/00_setup.R SuppFig1/panels/*.R SuppFig1/extra/*.R SuppFig1/assemble.R SuppFig1/run_all.R; do
  Rscript -e "invisible(parse('$f'))" || echo "PARSE FAIL: $f"
done
echo "(no PARSE FAIL above = clean)"
```
Expected: `All test files passed.`, no `PARSE FAIL`.

- [ ] **Step 5: Confirm nothing sensitive is staged**

Run:
```bash
git status --short
git check-ignore SuppFig1/objects/A.rds SuppFig1/data/merged_tpms.rds && echo "caches ignored"
```
Expected: only the intended docs; `caches ignored`.

- [ ] **Step 6: Commit**

```bash
git add SuppFig1/README.md SuppFig1/NOTES.md
git add -u SuppFig1
git commit -m "Document SuppFig1 and retire the notebook export

README covers the six-tag layout, the extras, the https-staged public inputs and
the requirements. NOTES records the TF double-log fix, the dropped isoform
histogram, the GO_annotations hook and the removed SE column. Deletes SuppFig1.R."
```

---

## Final verification

- [ ] **Push and open the PR**

```bash
git push -u origin suppfig1-cleanup
gh pr create --base main --head suppfig1-cleanup \
  --title "Clean up Supplemental Figure 1: AoU-vs-GTEx expression pipeline" \
  --body "$(cat <<'BODY'
Replaces the 316-line `SuppFig1/SuppFig1.R` with a modular, manifest-driven
pipeline, extends `common/manifest.R` to stage https inputs, and adds
`SuppFig1/SuppFig1.ipynb`.

Spec: `docs/superpowers/specs/2026-07-24-suppfig1-cleanup-design.md`
Original preserved in history (committed before removal).

## Results change — see `SuppFig1/NOTES.md`
- Panel C's TF statistics were logged **twice** in the original; now logged once.

## Structure
Four memoized builders; four panels (six tags A–F, pinned by a rendered-tag test);
four retained `extra/` scripts; a shared-loader https extension.

## Merge order
Stacked on `fig4-cleanup`. Merge `figure-notebooks` (#4) and `fig4-cleanup` (#5)
first; both touch `common/manifest.R`/`tools/build_notebooks.py` — take this
branch's superset if git flags either.

## Verified here
All tests pass; the pipeline assembles end-to-end against a mock fixture (Panel B's
live enrichR call is exercised on the Workbench). 

🤖 Generated with [Claude Code](https://claude.com/claude-code)
BODY
)"
```

---

## Traceability — spec section → task

| Spec § | Requirement | Task |
|---|---|---|
| §3 | Target structure | 2, 5–8, 10 |
| §4 | https loader extension | 1 |
| §5 | Folded helpers + four builders | 3, 4 |
| §6 | Panels + six-tag assembly + extras | 5, 6, 7, 8 |
| §7.1 | TF double-log fix | 6 |
| §7.2 | Cleanup (dead staging, dup def, after_stat, repr options, libs, SE) | 1, 2, 3 |
| §7.3 | Isoform dropped; GO hook | 7, 10 |
| §8 | Input manifest | 2 |
| §9 | Tests | 1, 3–8 |
| §10 | Notebook | 9 |
