#!/usr/bin/env python3
"""
Generate an R-kernel Jupyter notebook per figure that reproduces the figure's
results, with authored markdown explaining each analysis step.

The code cells are pulled from the actual pipeline scripts at build time, so the
notebooks stay faithful to the code and can be regenerated after any change:

    python3 tools/build_notebooks.py          # all figures
    python3 tools/build_notebooks.py Fig4     # just one figure

Generating a figure's notebook is a **standard part of cleaning that figure up**:
when a figure gains a modular pipeline, add its config to FIGURES below and commit
`FigN/FigN.ipynb` alongside the pipeline.

Design
------
* One notebook per cleaned figure: Fig1/Fig1.ipynb, Fig2/Fig2.ipynb,
  Fig3/Fig3.ipynb.
* A single setup cell sources <Fig>/00_setup.R (case-robustly, so it also works
  on Linux/Terra where 'Fig3' != 'fig3'); every later code cell is the *body* of
  a panel / prep / assemble / supplemental script with its `source(...)` line
  removed and a trailing line that renders the object inline.
* Panels still `saveRDS()` their object to <Fig>/objects, exactly as the scripts
  do, because assemble.R reads those .rds files back.

The notebooks are meant to run on the *All of Us* Workbench, where the
controlled-access inputs live; they are not executed by this generator.
"""
import json, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# --- helpers -----------------------------------------------------------------

def md(text):
    return {"cell_type": "markdown", "metadata": {},
            "source": _lines(text.strip("\n") + "\n")}

def code(text):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": _lines(text.rstrip("\n") + "\n")}

def _lines(text):
    return text.splitlines(keepends=True)

def read_script(relpath, fig_lower, strip_header=True, display=None):
    """Return a script's body ready to drop into a code cell:
       - drop the leading contiguous comment/blank header block (moved to md),
       - drop the `source("figN/00_setup.R")` line (setup is done once up top),
       - append `display` on its own line so the object renders inline.
    """
    lines = (ROOT / relpath).read_text().splitlines()
    src_re = re.compile(r'^\s*source\("' + fig_lower + r'/00_setup\.R"\)\s*$')

    # drop leading header (comments + blanks) until first real code line
    i = 0
    if strip_header:
        while i < len(lines) and (lines[i].lstrip().startswith("#") or not lines[i].strip()):
            i += 1
    body = [ln for ln in lines[i:] if not src_re.match(ln)]
    text = "\n".join(body).strip("\n")
    if display:
        text += "\n\n" + display
    return text

def script_header(relpath):
    """First contiguous comment lines of a script, cleaned of leading '# '."""
    out = []
    for ln in (ROOT / relpath).read_text().splitlines():
        if ln.lstrip().startswith("#"):
            out.append(re.sub(r"^\s*#\s?", "", ln))
        elif ln.strip() == "":
            if out:
                break
        else:
            break
    # drop banner rules like ==========
    return [l for l in out if not set(l.strip()) <= {"=", "-"} and l.strip()]

def setup_cell_code(fig_lower, fig_cap, env):
    envs = ",\n           ".join(f'{k} = file.path(fig_dir, "{v}")' for k, v in env.items())
    setenv = f"Sys.setenv({envs})\n" if env else ""
    return f'''## Locate the figure directory. Terra / Linux is case-sensitive, so match
## either "{fig_cap}" or "{fig_lower}" and point the pipeline's env-overridable
## paths at it before sourcing the shared setup.
fig_dir <- list.files(".", pattern = "^({fig_cap}|{fig_lower})$", include.dirs = TRUE)
stopifnot("Run this notebook from the repository root (the folder holding {fig_cap}/)." =
          length(fig_dir) == 1)
{setenv}source(file.path(fig_dir, "00_setup.R"))
cat("Setup loaded from", fig_dir, "\\n")'''

# --- per-figure authored content --------------------------------------------

PANEL_NOTES = {
 ("Fig1","A"): "Panel A is a reserved slot for the hand-drawn study-design schematic. The code emits a blank placeholder so the assembled layout keeps A's position; the diagram is dropped in externally.",
 ("Fig1","B"): "Counts participants per genetic-ancestry group by joining the ancestry predictions to the combined QTL sample list, then plots the per-group totals — the cohort composition the rest of the figure conditions on.",
 ("Fig1","C"): "For every ancestry × assay, computes the proportion of tested genes with at least one significant QTL — FDR-controlled **within** each permutation file — and plots it against the group's sample size, so discovery is shown scaling with N. Expression sits beside splicing + protein.",
 ("Fig1","D"): "Counts eGenes / sGenes discovered here that are **novel** relative to GTEx, INTERVAL and TOPMed, broken down by gene type — what the *All of Us* cohort adds over prior consortia.",
 ("Fig1","E"): "Builds the cumulative distribution of the number of SuSiE credible sets per eGene (fraction of eGenes with ≥ *k* credible sets), summarizing the allelic complexity of expression regulation.",
 ("Fig1","F"): "Compares credible-set sizes across pQTL ancestry groups, with a UKBB reference point — a readout of fine-mapping resolution for proteins.",
 ("Fig1","G"): "Per assay, the proportion of credible sets resolved to a **single** causal variant — a direct comparison of fine-mapping resolution across expression, protein and splicing.",
 ("Fig1","H"): "Functional enrichment of fine-mapped causal variants across regulatory / coding annotation groups, per assay, testing where causal variants concentrate.",

 ("Fig2","A"): "Cumulative minor-allele-frequency distribution of fine-mapped causal variants (fraction with MAF ≤ *k*), pooled across expression and splicing (protein where available) — the allele-frequency spectrum of causal regulatory variation.",
 ("Fig2","B"): "For every fine-mapped variant (PIP > 0.9, one row per variant) assigns a primary population and plots the distribution of population enrichment, split rare vs. common, with counts of ancestry-enriched variants — how causal variants partition across ancestries.",
 ("Fig2","C"): "Effect size (|log2 allelic fold-change|) by variant molecular consequence, rare (dark) vs. common (light), annotated with the rare−common effect shift per consequence — do rare causal variants carry larger effects?",
 ("Fig2","D"): "Cumulatively combines a gene's independent causal variants to estimate expected heterozygous expression, contrasting deletion-like and duplication-like directions — a dosage model built from the fine-mapped effects.",
 ("Fig2","E"): "Expression-outlier enrichment across effect-size bins (under- vs. over-expression), from the QTL-burden summary — linking causal-variant burden to outlier expression.",
 ("Fig2","F"): "Median number of dosage-affected protein-coding genes per individual by percent-change bin, with the deletion-/duplication-like bins highlighted.",
 ("Fig2","G"): "Counts of validated dosage-like events, deletion- vs. duplication-like, drawn as a dumbbell (segment + two end points, no extra dependency).",

 ("Fig3","A"): "Among genes fine-mapped in **all three** assays, the percent of genes by the number of molecular layers that colocalize (0–3), from gene-level FastENLOC — how often a signal is shared across expression, protein and splicing. The '≥1 layer' figure in the annotation is computed from the data, not hard-coded.",
 ("Fig3","B"): "Collapses colocalized credible sets into connected components (a signal shared across layers, via an undirected graph) and plots the percent of genes by number of colocalized signals, with an inset giving the per-assay-pair credible-set colocalization rate.",
 ("Fig3","C"): "Fits, per assay pair, the interaction between −log(MAF) and cross-assay effect-size correlation, and plots each coefficient with its 95% CI — does the strength of a shared signal track allele frequency? (The three per-pair model formulas are deliberately not identical; see `NOTES.md`.)",
 ("Fig3","D"): "eQTL vs. pQTL fine-mapped effect sizes, faceted by MAF bin (Common / Low-Freq / Rare), with a per-bin slope and N — how expression and protein effects relate across the frequency spectrum.",
 ("Fig3","E"): "A reserved slot in the bottom row; its content is undecided, so the code renders a blank to keep the layout stable (see `Fig3/NOTES.md`).",

 ("Fig4","A"): "For each trait category, the number of colocalizing genes (intensity-shaded bars) beside the credible-set colocalization rate (lollipop, with a dashed median reference). The rate plotted is the **lenient** union rate (`pct_lenient_union`): `coloc_by_trait()` calls `build_plot_data(rate_type = \"lenient\")`, even though that function's own default is `\"stringent\"` — a caller reusing it gets the stringent rate unless it asks otherwise. The panel keeps **traits** with more than five colocalizing genes, then the top ten traits *within each category*, and drops the `neuropsychiatric` and `renal_urological` categories.",
 ("Fig4","B"): "Histogram of the R² increase (%) that rare variants add to TWAS, computed as the difference in adjusted R² between allele-frequency thresholds 0 and 0.01 for genes with a significant fit — how much predictive signal rare variants contribute. Note the axis is truncated with `xlim(-20, 40)`, so genes whose R² increase falls outside that window are dropped from the histogram rather than shown at the edge.",
 ("Fig4","C"): "Common and rare conditional-TWAS QQ plots over nine blood traits, comparing colocalized genes against **MatchIt propensity-matched** non-colocalized genes (\"Null\") — an inflation test that controls for expression, selection, rare-variant count and LD.",
 ("Fig4","D"): "Odds-ratio forest (logistic regression) for the enrichment of rare vs. common conditional-TWAS signal across four gene-set outcomes — Mendelian (any), Mendelian (immune/hematological), druggable, and approved drug target. Both continuous and binary predictor sets are fit; the panel plots the continuous rare/common terms (see `Fig4/NOTES.md`).",

 ("SuppFig1","A"): "Protein-coding genes detected in either cohort: AoU-vs-GTEx log2(TPM+1) correlation (hex) beside counts of genes detected in only one cohort. Sub-plots are tagged A and B.",
 ("SuppFig1","B"): "Top GO terms (by adjusted p) enriched among genes specific to each cohort, faceted AoU vs GTEx, with a -log10(FDR)=1.3 reference line. Tagged C.",
 ("SuppFig1","C"): "Transcription-factor expression, AoU vs GTEx, as log2(TPM+1) boxplots, with an lm and a paired t-test reported in the run log. The TPMs are logged exactly once (the original logged twice). Tagged D.",
 ("SuppFig1","D"): "lncRNA AoU-vs-GTEx log2(TPM+1) correlation (hex) beside cohort-specific detection counts. Sub-plots are tagged E and F.",

 ("Fig5","A"): "Fine-mapped rare ClinGen variants by genetic ancestry, faceted AD vs AR (lighter = novel, darker = ClinVar-annotated), over a per-inheritance ClinVar-classification proportion bar. A 2-plot patchwork; tagged A.",
 ("Fig5","B"): "|Average QTL effect| by variant-effect category (Control / AD / AR, split Aggregate / Decreasing / Increasing) with 95% CIs. Plots |posterior_mean| (the SuSiE effect), not log2 aFC. A gene curated under BOTH AD and AR is counted in both classes (preserved fan-out; see Fig5/NOTES.md). Tagged B.",
 ("Fig5","C"): "log2 aFC of the strongest decreasing variant per gene, top 7 per AD/AR (HEXA and MLC1 forced in), optionally refined by the private no-homozygote list. Tagged C.",
 ("Fig5","D"): "ACVRL1 eQTL fine-mapping: regional PIP with a zoom inset and an aligned gene-model track, ENCODE cCREs shaded. A patchwork wrapped as one tag leaf; tagged D.",
 ("Fig5","E"): "ACVRL1 causal promoter variants — PromoterAI base-resolution effect on a signed-square axis, one window per causal variant; the focal chr12:51,912,316 C>T sits in an SP1 GC-box. Tagged E.",
 ("Fig5","F"): "SPR eQTL fine-mapping: regional PIP + zoom inset + gene-model track. A patchwork wrapped as one tag leaf; tagged F.",
 ("Fig5","G"): "SPR Borzoi multi-track allelic browser (reference vs alternate across 16 marks, broad→narrow zoom) with the CTCF base-resolution attribution logo. A cowplot object; tagged G.",
}

FIGURES = {
 "Fig1": {
   "num": 1, "cap": "Fig1", "lower": "fig1", "obj_letters": list("ABCDEFGH"),
   "title": "Figure 1 — multi-omic QTL discovery across ancestries in *All of Us*",
   "overview": """
This notebook reproduces **Figure 1** and its supplemental. The figure summarizes
cis-QTL discovery across three molecular assays — expression (eQTL), protein
(pQTL) and splicing (sQTL) — in the ancestrally diverse *All of Us* cohort:
cohort composition, the proportion of genes with significant QTLs by sample size,
novelty relative to prior consortia, and the SuSiE fine-mapping architecture of
causal variants.

It runs the modular pipeline in `Fig1/` step by step: a shared **setup**, then
one cell per **panel** (each loads only the inputs it needs and caches its ggplot
object), then **assembly** into the final figure.
""",
   "data_note": """
### Inputs

> **Controlled-access.** Every input derives from *All of Us* controlled-access
> data and is not in this repository. This notebook is meant to run on the
> Researcher Workbench, where the data lives.

Inputs are **manifest-driven** (`Fig1/inputs.tsv`). Each panel names the logical
inputs it needs by *key*; the loaders (`input_path`, `load_input`, `input_files`)
resolve a key to file(s) and stage them from Google Cloud with `gsutil` **only if
the local copy is missing**. Pre-staged files are used as-is.
""",
   "env": {"FIG1_MANIFEST": "inputs.tsv", "FIG1_OBJ": "objects"},
   "steps": [
     ("panel","Fig1/panels/panelA.R","A","pA"),
     ("panel","Fig1/panels/panelB.R","B","pB"),
     ("panel","Fig1/panels/panelC.R","C","pC"),
     ("panel","Fig1/panels/panelD.R","D","pD"),
     ("panel","Fig1/panels/panelE.R","E","pE"),
     ("panel","Fig1/panels/panelF.R","F","pF"),
     ("panel","Fig1/panels/panelG.R","G","pG"),
     ("panel","Fig1/panels/panelH.R","H","pH"),
     ("assemble","Fig1/assemble.R",None,"Fig1"),
     ("supp","Fig1/supp_metasoft.R",None,None,
        "**Supplemental (MetaSoft).** Cross-ancestry meta-analysis M-value summary. "
        "This step needs the external **MetaSoft** Java tool (`METASOFT_JAR`) and is "
        "skipped if its outputs are absent; it is not part of the main figure."),
   ],
 },
 "Fig2": {
   "num": 2, "cap": "Fig2", "lower": "fig2", "obj_letters": list("ABCDEFG"),
   "title": "Figure 2 — causal-variant genetics across expression & splicing QTLs",
   "overview": """
This notebook reproduces **Figure 2** and its supplemental. The figure summarizes
properties of statistically fine-mapped causal variants (SuSiE, PIP > 0.9) across
molecular QTL assays in *All of Us*: their allele-frequency spectrum, ancestry
enrichment, molecular effect sizes, and their relationship to expression-outlier
and dosage-like phenotypes.

Figure 2 uses a **prep-then-panels** pattern: `01_prep_data.R` reads the raw
inputs once and derives every per-panel data frame into `Fig2/data/*.rds`; each
panel then reads its own prepared frame, so iterating on one panel never
re-reads the large SuSiE tables. Run the prep cell before the panel cells.
""",
   "data_note": """
### Inputs

> **Controlled-access.** The raw inputs derive from *All of Us* controlled-access
> data and are not in this repository. The prep step reads them **by manifest key**
> (`Fig2/inputs.tsv`) from the staging root `FIG2_INPUT` (default the cwd; the files
> live under `input_data/`). The derived per-panel `Fig2/data/*.rds` are
> participant-derived and git-ignored, so run the prep cell before the panels.
""",
   "env": {"FIG2_OBJ": "objects", "FIG2_DATA": "data"},
   "steps": [
     ("prep","Fig2/01_prep_data.R",None,None,
        "**Data prep.** Reads the raw SuSiE / QTL-burden inputs once and writes one "
        "prepared data frame per panel to `Fig2/data/<panel>.rds`. Run this before the "
        "panels (or rely on the committed `Fig2/data/`)."),
     ("panel","Fig2/panels/panelA.R","A","pA"),
     ("panel","Fig2/panels/panelB.R","B","pB"),
     ("panel","Fig2/panels/panelC.R","C","pC"),
     ("panel","Fig2/panels/panelD.R","D","pD"),
     ("panel","Fig2/panels/panelE.R","E","pE"),
     ("panel","Fig2/panels/panelF.R","F","pF"),
     ("panel","Fig2/panels/panelG.R","G","pG"),
     ("assemble","Fig2/assemble.R",None,"Fig2"),
     ("supp","Fig2/supp_ancestry_stacked.R",None,"stacked_histogram",
        "**Supplemental (ancestry composition).** Composition of ancestry-enriched "
        "fine-mapped expression variants across AoU allele frequency, stacked by "
        "enriched ancestry, with the proportion shared across all continents overlaid."),
   ],
 },
 "Fig3": {
   "num": 3, "cap": "Fig3", "lower": "fig3", "obj_letters": list("ABCDE"),
   "title": "Figure 3 — cross-QTL colocalization across molecular layers",
   "overview": """
This notebook reproduces **Figure 3** and its two supplementals. The figure asks
how often a genetic signal is **shared between molecular layers** — expression,
protein and splicing — using FastENLOC colocalization between each assay pair at
gene, credible-set (signal) and variant level, and whether that sharing depends on
allele frequency.

Figure 3's panels share a heavy upstream, so `00_setup.R` provides **memoized
derived-data builders** (`susie_finemapping()`, `coloc_gene_level()`,
`coloc_cs_level()`, `coloc_variant_level()`, `coloc_cs_components()`,
`coloc_model_input()`) that each build once and cache to `Fig3/data/`. Panels call
only the builders they need.
""",
   "data_note": """
### Inputs

> **Controlled-access.** Every input derives from *All of Us* controlled-access
> data and is not in this repository; run this on the Workbench.

Inputs are **manifest-driven** (`Fig3/inputs.tsv`): the AoU annotation table
auto-stages from its `gs://` URI, while the SuSiE and FastENLOC tables are marked
local-only and must already be staged. The derived-data builders cache to
`Fig3/data/` **by name**, so clear that directory (or pass `refresh = TRUE`) after
changing an input.

> **Results differ from the pre-cleanup figure.** Several fixes change output —
> most visibly Panel B, whose colocalization graph was previously built from one
> assay pair. See `Fig3/NOTES.md`.
""",
   "env": {"FIG3_MANIFEST": "inputs.tsv", "FIG3_OBJ": "objects", "FIG3_DATA": "data"},
   "steps": [
     ("panel","Fig3/panels/panelA.R","A","pA"),
     ("panel","Fig3/panels/panelB.R","B","pB"),
     ("panel","Fig3/panels/panelC.R","C","pC"),
     ("panel","Fig3/panels/panelD.R","D","pD"),
     ("panel","Fig3/panels/panelE.R","E","pE"),
     ("assemble","Fig3/assemble.R",None,"Fig3"),
     ("supp","Fig3/supp_coloc_qc.R",None,"SuppColocQC",
        "**Supplemental (colocalization QC).** GLCP / RCP / SCP posterior "
        "distributions by assay pair, over colocalization counts across posterior "
        "thresholds."),
     ("supp","Fig3/supp_prop_gene_coloc.R",None,"ProportionGeneColocs",
        "**Supplemental (gene colocalization rate).** Proportion of genes "
        "colocalized, by assay pair."),
   ],
 },
 "Fig4": {
   "num": 4, "cap": "Fig4", "lower": "fig4", "obj_letters": list("ABCD"),
   "title": "Figure 4 — QTL–GWAS colocalization and rare-variant TWAS",
   "overview": """
This notebook reproduces **Figure 4**. The figure connects molecular QTL
colocalization with GWAS traits to downstream rare-variant signal: how many genes
colocalize per trait category and at what lenient credible-set rate, the R² that
rare variants add to TWAS, whether colocalized genes carry inflated conditional-TWAS
signal versus matched controls, and whether that signal enriches in clinically
important gene sets.

Like Figure 3, `00_setup.R` provides **memoized derived-data builders**
(`read_conditional_twas()`, `gene_sets()`, `coloc_by_trait()`,
`qq_match_source()`, `gene_twas_data()`) that each build once and cache to
`Fig4/data/`. Panels call only the builders they need.
""",
   "data_note": """
### Inputs

> **Controlled-access.** Every input derives from *All of Us* controlled-access
> data and is not in this repository; run this on the Workbench.

Inputs are **manifest-driven** (`Fig4/inputs.tsv`). The two
`high_level_coloc_outputs/` tables are extracted from a single archive
(`gsutil cp gs://…/high_level_coloc_outputs.tar.gz . && tar -xf …`); the ClinGen
parquet, drug-interaction table, R² and conditional-TWAS tables are local-only and
must be pre-staged. The derived-data builders cache to `Fig4/data/` **by name**,
so clear that directory (or pass `refresh = TRUE`) after changing an input.

> **Results differ from the pre-cleanup script.** The duplicated
> conditional-TWAS table (Panel D) is deduplicated to its clean version, and the
> unused GeneBayes read is dropped. See `Fig4/NOTES.md`.
""",
   "env": {"FIG4_MANIFEST": "inputs.tsv", "FIG4_OBJ": "objects", "FIG4_DATA": "data"},
   "steps": [
     ("panel","Fig4/panels/panelA.R","A","coloc_plot"),
     ("panel","Fig4/panels/panelB.R","B","R2_dist"),
     ("panel","Fig4/panels/panelC.R","C","p_blood_qq"),
     ("panel","Fig4/panels/panelD.R","D","p"),
     ("assemble","Fig4/assemble.R",None,"Fig4"),
   ],
 },
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
 "Fig5": {
   "num": 5, "cap": "Fig5", "lower": "fig5", "obj_letters": list("ABCDEFG"),
   "title": "Figure 5 — rare-variant eQTL constraint & sequence-model variant-effect dissection",
   "overview": """
This notebook reproduces **Figure 5**. The figure pairs a population / constraint
view of fine-mapped **rare** eQTLs at ClinGen disease genes (A–C) with two worked
sequence-model dissections of causal regulatory variants: the **ACVRL1** promoter
(PromoterAI, D/E) and the **SPR** regulatory element (Borzoi allelic browser, F/G).

It runs the modular pipeline in `Fig5/` step by step: a shared **setup** with
memoized loaders (`get_reference()`, `get_eqtl()`, `get_variant_freq()`), one cell
per **panel** (each builds a single *untagged* object), then **assembly**, where
the panel tags A–G are applied. Additional AlphaGenome analyses and supplements
live in `Fig5/extra/` and are not run here.
""",
   "data_note": """
### Inputs

> **Controlled-access.** The SuSiE eQTL table and the no-homozygote list derive
> from *All of Us* controlled-access data and are not in this repository; the SPR
> Borzoi bigWigs are large model-derived files. Public / derived reference tables
> (ClinGen, PromoterAI, ENCODE cCREs, GENCODE gene models) are committed under
> `Fig5/reference/`.

Inputs are **manifest-driven** (`Fig5/inputs.tsv`). Panel C degrades gracefully
when the private no-homozygote list is absent.

> **Panel B/C note.** A gene curated under **both** AD and AR ClinGen inheritance
> is counted under both classes (preserved from the original). See `Fig5/NOTES.md`.
""",
   "env": {"FIG5_MANIFEST": "inputs.tsv", "FIG5_OBJ": "objects", "FIG5_DATA": "data"},
   "steps": [
     ("panel","Fig5/panels/panelA.R","A","A"),
     ("panel","Fig5/panels/panelB.R","B","B"),
     ("panel","Fig5/panels/panelC.R","C","C"),
     ("panel","Fig5/panels/panelD.R","D","D"),
     ("panel","Fig5/panels/panelE.R","E","E"),
     ("panel","Fig5/panels/panelF.R","F","F_panel"),
     ("panel","Fig5/panels/panelG.R","G","G"),
     ("assemble","Fig5/assemble.R",None,"Fig5"),
   ],
 },
}

# --- build -------------------------------------------------------------------

def build(figkey):
    F = FIGURES[figkey]
    lower, cap = F["lower"], F["cap"]
    cells = []

    cells.append(md(f"# {F['title']}\n\n{F['overview']}"))
    cells.append(md(
        "## How to run\n\n"
        "* **Kernel:** R (IRkernel). **Working directory:** the repository root "
        f"(the folder that contains `{cap}/`).\n"
        "* Run the cells top to bottom. Each panel caches its ggplot object to "
        f"`{cap}/objects/`, and the assembly cell reads those back.\n"
        "* Re-running a single panel cell rebuilds just that panel; re-run the "
        "assembly cell to recompose.\n"
    ))
    cells.append(md(F["data_note"]))

    cells.append(md(
        "## Setup\n\n"
        "`00_setup.R` is the shared foundation every step builds on: it loads the "
        "libraries, sets the **Inter** font with a matched output **DPI**, defines "
        "the `theme_set`/`SZ_*` design system and the env-overridable paths, and "
        + ("brings in the manifest loaders and the folded-in analysis / plot "
           "functions" if figkey != "Fig2" else "brings in the shared helper "
           "functions")
        + (" plus the memoized derived-data builders" if figkey in ("Fig3", "Fig4", "SuppFig1") else "")
        + " — so the figure is self-contained."
    ))
    cells.append(code(setup_cell_code(lower, cap, F["env"])))

    if figkey in ("Fig1", "Fig3", "Fig4", "SuppFig1"):
        cells.append(md(
            "### Panels\n\nEach panel below **loads only the inputs it needs** "
            "(staging them on demand), does its wrangling, builds one ggplot object, "
            f"caches it to `{cap}/objects/`, and renders it inline."
        ))

    for step in F["steps"]:
        kind, path = step[0], step[1]
        display = step[3]
        extra_md = step[4] if len(step) > 4 else None

        if kind == "panel":
            letter = step[2]
            hdr = script_header(path)
            inputs_line = next((h for h in hdr if h.lower().startswith("inputs")), None)
            # the title sentence may wrap across several comment lines before "Inputs:"
            title_lines = hdr[:hdr.index(inputs_line)] if inputs_line in hdr else hdr
            title = " ".join(title_lines) if title_lines else f"Panel {letter}"
            note = PANEL_NOTES.get((figkey, letter), "")
            body = f"### {title}\n\n{note}"
            if inputs_line:
                body += f"\n\n*{inputs_line}*"
            cells.append(md(body))
            cells.append(code(read_script(path, lower, display=display)))

        elif kind == "prep":
            cells.append(md("### Data prep\n\n" + extra_md))
            cells.append(code(read_script(path, lower, display=display)))

        elif kind == "assemble":
            cells.append(md(
                "## Assemble the figure\n\n"
                f"Loads the cached panel objects from `{cap}/objects/` and composes "
                "them with patchwork into the final figure. Panel **tagging happens "
                "here and only here**. Writes the PDF plus a PNG preview, and renders "
                "the assembled figure inline."
            ))
            cells.append(code(read_script(path, lower, display=display)))

        elif kind == "supp":
            cells.append(md("## " + extra_md))
            cells.append(code(read_script(path, lower, display=display)))

    nb = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"name": "ir", "display_name": "R", "language": "R"},
            "language_info": {"name": "R", "codemirror_mode": "r",
                              "file_extension": ".r", "mimetype": "text/x-r-source",
                              "pygments_lexer": "r", "version": "4.3"},
        },
        "nbformat": 4, "nbformat_minor": 5,
    }
    out = ROOT / cap / f"{cap}.ipynb"
    out.write_text(json.dumps(nb, indent=1))
    return out, len(cells)


if __name__ == "__main__":
    # No args -> build every figure's notebook. Pass figure keys (e.g. "Fig4")
    # to build only those -- the per-figure regeneration used during a cleanup.
    import sys
    targets = sys.argv[1:] or list(FIGURES.keys())
    unknown = [t for t in targets if t not in FIGURES]
    if unknown:
        sys.exit(f"unknown figure(s): {', '.join(unknown)}; "
                 f"choose from {', '.join(FIGURES)}")
    for k in targets:
        out, n = build(k)
        print(f"wrote {out.relative_to(ROOT)}  ({n} cells)")
