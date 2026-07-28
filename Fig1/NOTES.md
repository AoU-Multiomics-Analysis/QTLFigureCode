# Figure 1 — deferred analyses & notes

## Resolved: novel pGene comparison (Panel D)

Panel D shows **novel eGene + sGene + pGene** counts relative to external
consortia. The pGene bar was previously absent: the original notebook referenced
`UKBB_pqtl_path` but never assigned it, so the third bar silently disappeared --
a latent bug rather than a deliberate omission.

It is now wired up. The UKBB-PPP reference is the Nature 2023 supplementary
workbook, declared as `ukbb_pqtl` in `fig1/inputs.tsv` and read by
`LoadUKBBpQTL()` from **sheet 10 (`ST9`), skipping 4 header rows**, filtered to
`cis/trans == "cis"`.

**This changed Panel D.** It gained a third bar; the eGene and sGene bars are
unaffected. Counts on the staged data:

| gene class | eGene | sGene | pGene |
|------------|------:|------:|------:|
| Protein-Coding | 1,207 | 2,035 | 1,783 |
| lncRNA         | 8,696 |   415 |     0 |
| non-coding RNA |   246 |     0 |     0 |

`novel_pQTLs` drops lncRNA explicitly: the proteomics assay cannot measure it,
so a non-zero "novel lncRNA pGene" count would be an artefact of the gene-type
join rather than a finding. The two zeros in the pGene column are therefore
expected, not missing data.

Note this is a **different** UKBB reference from the one on Panel F, which
hard-codes `mean_cs = 9.7` as its credible-set comparison point
(`panels/panelF.R`). That constant is unrelated to this workbook and is still
not derived from data.

## Panel H: annotation source, VEP expansion, and a casing bug

Three compounding changes. **Panel H's result changed** — it went from 3
annotations to 12.

### 1. `annot` now reads the top-level SuSiE tables

The key pointed at per-assay `AnnotateSusie` outputs, one of which
(`annotated/pQTL_susie_annotated.tsv`) had no recorded origin anywhere, so the
panel was permanently unbuildable. It now reads the same three tables Fig2 and
Fig3 use:

    COMB.{expression,splicing,proteomics}_SusieMerged.annotated.ancestry_skew.tsv.gz

All three figures therefore share one SuSiE provenance chain. The `gs_uri` on
these rows was cleared: the old bucket paths point at the superseded
`AnnotateSusie` outputs and would re-stage the wrong files on a clean checkout.

### 2. VEP consequences are expanded before the enrichment

`ComputeFinemappingAnnotationEnrichment()` picks annotation columns up via
`where(is.logical)`. These tables carry consequences as ONE comma-separated
string (`"intron_variant,NMD_transcript_variant"`), ~20% of rows holding more
than one term, so the entire coding / splice / UTR / noncoding-RNA half of
`variant_functional_groups` was invisible: the panel rendered with three cCRE
annotations (PLS, pELS, dELS) and no error.

`.expand_vep_consequences()` in `00_setup.R` one-hots each term into its own
logical column, applied after the `pip` filter so it runs on the reduced set. It
builds the matrix directly rather than `separate_rows()` + `pivot_wider()` —
these tables have ~80 columns and pivot_wider would treat every one as an id
column. It drops VEP's `.` placeholder and refuses to overwrite an existing
column (the cCRE flags are already one-hot).

### 3. Casing bug: both UTR bars were being dropped

`PlotFinemappingEnrichment()` filters against a 12-annotation whitelist by
equality. The whitelist and the tribble spelled the UTR terms
`3_prime_utr_variant` / `5_prime_utr_variant`; VEP emits `3_prime_UTR_variant` /
`5_prime_UTR_variant`.

The two UTR enrichments were computed correctly (3,025 and 1,696 variants,
both significant) and then **silently discarded at the join** — which is exactly
why only 10 of the 12 whitelisted annotations appeared. Names now match VEP.

The whitelist in `annotation_labels` is exactly these twelve:
`splice_donor_variant`, `splice_acceptor_variant`, `splice_region_variant`,
`stop_gained`, `missense_variant`, `inframe_insertion`, `inframe_deletion`,
`5_prime_UTR_variant`, `3_prime_UTR_variant`, `PLS`, `pELS`, `dELS`.

Note there is **no microRNA (or any other noncoding-RNA) entry**, so the
"Noncoding RNA" group defined in `variant_functional_groups`
(`non_coding_transcript_exon_variant`, `non_coding_transcript_variant`,
`mature_miRNA_variant`) is unreachable: those rows are dropped by the
`name %in% names(annotation_labels)` filter and the group can never render.
Widening the panel to noncoding RNA means adding labels to `annotation_labels`,
not just to the tribble.

This bug most likely predates the change: the old `AnnotateSusie` files probably
carried lowercase one-hot columns, so the whitelist matched. Repointing to the
`ancestry_skew` exports, which preserve VEP's native spelling, exposed it.

Both failures share a shape worth remembering: **an annotation whose name does
not match is dropped, never reported.** The panel built cleanly at 3 annotations
and again at 10.

## Panel D: the splicing and protein colours were swapped (rendered tag **E**)

Three panels colour by assay, and each lists its categories in a different
order: panel D as eGene / sGene / pGene, panel G as Expression / Protein /
Splicing, panel H as eQTL / pQTL / sQTL. All three took their colours
positionally off `ltc::minou`, so panel D handed **sGene the protein colour
(red) and pGene the splicing colour (yellow)** while every other panel used red
for protein and yellow for splicing. Read across the figure, the same assay
changed colour between panels.

`AssayColors` in `00_setup.R` now maps every alias of an assay
(`Expression`/`eQTL`/`eGene`, `Protein`/`pQTL`/`pGene`,
`Splicing`/`sQTL`/`sGene`) to one colour by name, and panels D, G and H use it
via `scale_fill_manual()`. **Panel D's bars changed colour** — sGene is now
yellow and pGene red; bar order, legend order and every count are unchanged.
Panels G and H are unchanged, they already agreed with the palette.

`PlotNovelQTLCount()` also carried a `scale_color_*` alongside its fill scale
even though it maps no colour aesthetic; it is dropped.

## Log axes are labelled as powers of ten (rendered tags **G** and **I**)

`PlotCSSizes()` and the enrichment panel now pass `scales::label_log()` to their
log10 scales, so ticks read 10^0, 10^1, 10^2 … instead of 1, 10, 100. Breaks and
data are unchanged. The novel-gene panel (rendered **E**) is also log10 and was
deliberately left with plain numeric labels, because its bars are annotated with
exact counts and the two would compete.

## Panel E y axis now follows the data (rendered tag **F**)

`panels/panelE.R` pinned the y axis to `limits = c(0, 1)`, but the cumulative
curve only ever runs from 100% (k = 1) down to ~49% (k = 5+), so half the panel
was empty. The hard limits are gone; the axis takes its range from the data with
`expansion(mult = 0.1)` — enough padding that the size-10 points at 100% and at
the bottom of the curve are not clipped — and breaks are labelled every 10 points
from 50% to 100%.

No value changed; only the visible range did. Anyone reading percentages off
this panel should note the axis no longer starts at zero.

## SUPERSEDED — Panel H: infinite odds ratios drawn, not clipped (tag **I**)

> Kept because it is why the panel changed method. Panel H no longer uses the
> Fisher exact test at all; see *Panel H is now a PIP logistic regression*
> below. `PlotFinemappingEnrichment()` / `RenderFinemappingEnrichment()` still
> exist in `00_setup.R` and still behave as described here, but nothing in the
> figure calls them.


Three pQTL annotations — **inframe insertion**, **stop gained** and **splice
acceptor** — have an *infinite* Fisher odds ratio and an infinite upper
confidence bound, because a cell of their 2x2 table is empty (no PIP < 0.01
background variant carries the annotation). `ggplot` placed those points at
`+Inf`, i.e. exactly on the panel edge, where they were drawn as clipped half
circles with intervals running off the top.

They are now drawn at a finite ceiling, `1.5 x` the largest bounded value on the
panel (`max` over all finite `odds_ratio` and `conf.high`, ~5,503 → ceiling
~8,255), as an **upward-pointing triangle** rather than a circle, with the
confidence interval drawn from its finite lower bound up to the same ceiling.
The triangle marks an estimate that is unbounded above; its height carries no
quantitative meaning. The y scale gained upper expansion (1.5% → 6%) so the
markers sit inside the panel.

Two smaller consequences:

- **The x order inside the Coding facet changed.** Annotations are ordered by
  mean odds ratio; with three `Inf` values the mean was `Inf` for every
  annotation containing one, so ties fell back to input order. `Inframe
  insertion` and `Stop gained` swap places. No other facet moved.
- The alternating grey bands now use `ymax = Inf` instead of a computed
  `max(conf.high) * 1.5`, which was `Inf` in practice anyway. The lower edge is
  still computed, because it is what keeps the `y = 1` reference line inside the
  axis range.

The underlying enrichment statistics are untouched — no continuity correction
was applied, and the three odds ratios are still genuinely infinite. If the
manuscript needs finite estimates for these annotations, that is a change to
`ComputeFinemappingAnnotationEnrichment()` (e.g. a Haldane–Anscombe correction),
not to the plot.

The rendering half of `PlotFinemappingEnrichment()` was split out as
`RenderFinemappingEnrichment(plot_data)`, a pure function of the prepared
enrichment table, so the panel can be redrawn from an already-computed table
without the controlled-access inputs.

## Panel H is now a PIP logistic regression (rendered tag **I**)

**Panel H's result changed.** The enrichment is no longer a Fisher exact test.

| | before | now |
|---|---|---|
| model | Fisher exact test on (PIP > 0.9) × annotated | `glm(annotated ~ PIP, binomial)` |
| variants used | only PIP > 0.9 or PIP < 0.01 | all of them |
| effect | odds ratio | slope, shown exponentiated |
| axis | log10 fold enrichment | log10, reference line at 1 |
| infinite estimates | 3, drawn as capped triangles | none |

The effect is the change in the log odds of carrying the annotation as PIP goes
from 0 to 1. It is plotted `exp()`-transformed on a log10 axis, so it reads as
an odds ratio between PIP = 1 and PIP = 0 and the axis matches the figure's
other log panels. The y axis is labelled *Fine-mapping Enrichment (Log Odds)*.

The two methods agree closely: across the 34 annotation × assay cells both
report, **the sign agrees 34/34** and Spearman ρ between the Fisher log odds
ratio and the regression slope is **0.965** over the 31 cells with a finite odds
ratio.

The reason for switching is the other three cells. The pQTL **stop gained**,
**inframe insertion** and **splice acceptor** enrichments had an infinite Fisher
odds ratio (see the superseded note above, where they are drawn as triangles at
a capped ceiling). The regression gives all three a finite, tightly determined
slope — 4.80, 4.14 and 3.74 — because it uses the whole PIP range instead of a
2×2 table whose background cell happens to be empty. Separation does not vanish
entirely: `stop_lost` (pQTL) and `start_retained_variant` (sQTL) have no finite
maximum-likelihood slope. Those are marked `separated = TRUE` in the results and
dropped from the panel with a message rather than drawn with an axis-destroying
interval. Neither is on the panel's whitelist, so today nothing is dropped.

The panel also gained two cells the Fisher version could not show at all
(34 → 36), because nothing is discarded for having an empty background cell.

### Which annotations the panel shows

`PANEL_ANNOTATION_LABELS` in `00_setup.R` is the whitelist — membership is what
puts an annotation on the panel. The coding slots were **inframe insertion** and
**inframe deletion**; they are now **frameshift** and **start lost**, which are
significant in all three assays, whereas inframe insertion was not significant
for sQTL. The other nine annotations are unchanged, and all 36 remaining
annotation × assay cells are FDR < 0.05.

Points that are not significant are drawn faded, and the significance key only
appears in the legend when something is actually faded — with everything
significant it would be a legend entry that never varies.

### Presentation

Following Fig2 panel B, the x axis line and ticks are gone and each annotation
gets a thin grey rule beneath the span its three dodged assay points occupy. The
rule sits below the odds-ratio-of-1 reference line: taking its height from the
data alone put it on top of that line whenever no annotation was depleted, which
is the usual case.

The alternating grey background bands the Fisher panel used are gone — the rules
delimit the categories on their own, as they do in Fig2 panel B, and the two
together read as clutter. The y axis title is set on two lines; at this figure's
axis-title size it otherwise runs nearly the full panel height.

One trap worth recording even though the code no longer hits it, because it
fails quietly and the bands may come back. Full-height `geom_rect` bands on this
panel must use `ymin = 0`, not `ymin = -Inf`: the axis is log10-transformed,
`log10(-Inf)` is `NaN`, and every band is then dropped with nothing but a
"Removed 5 rows containing missing values" warning. `log10(0)` is `-Inf`, which
is the panel edge, which is what was wanted.

### Where the fits live

`FitAnnotationEnrichmentGlm()` memoizes to `fig1/objects/enrichment_glm.tsv`,
shared between panel H and the standalone `enrichment_glm.R` figure, so
whichever runs first pays for the fit (a few minutes, dominated by the eQTL
table). Delete that file or pass `refresh = TRUE` after changing the model or
the inputs. `enrichment_glm.R` draws all 31 annotations rather than the panel's
twelve, on the same scale.

### `where(is.logical)` does not find the annotation columns

Worth recording because it cost a debugging cycle and the same idiom is still
live in `ComputeFinemappingAnnotationEnrichment()`. In the **splicing** table
`tss`, `gene_name`, `distTSS`, `Constrained` and `lof.pLI` are `NA` for *every
row* — not merely sparse — so `fread` types all five as logical and
`where(is.logical)` collects them as annotations. One is then all-`NA` after the
pivot, which leaves `glm()` with an empty response and aborts the run outright.
The expression and proteomics tables have real values in those columns, so they
are typed integer/character/numeric there and the problem does not appear.

`enrichment_glm_fit.R` therefore names its seven flag columns
(`CCRE_FLAG_COLS`) and errors if any is absent or non-logical. The idiom is
still live in `ComputeFinemappingAnnotationEnrichment()`, which now feeds
nothing in the figure; it survived the same hazard only because its 12-name
whitelist discarded the junk columns downstream — not robust, just lucky.

**The splicing SuSiE table having no gene-level annotation at all** (`tss`,
`gene_name`, `distTSS` entirely `NA`) may matter beyond this script; anything
joining splicing variants to genes through those columns will silently get
nothing.

## MetaSoft supplemental

`fig1/supp_metasoft.R` reproduces the cross-ancestry meta-analysis summary
(`xGeneMetasoftSummary.pdf`). It is **not** part of Figure 1 and is not run by
`run_all.R`. It requires the MetaSoft Java tool (`METASOFT_JAR`); see the script
header.

## Original source

Reconciled from `QTLAnalysisUtils/notebooks/QTL_paper_Fig1.ipynb` (the authoritative
working notebook). The pre-refactor `Fig1.R` was a linearized export of that
notebook and contained several bugs (undefined assembly variables, a `ptions()`
typo, "casual"→"causal" misspellings, and duplicated read logic), all fixed here.
