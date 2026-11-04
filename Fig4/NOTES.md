# Figure 4 — notes and deferred work

## Panel A was rebuilt from a newer coloc run (2026-07-28)

The `high_level_coloc_outputs` copy under the usual staging root was known bad,
so Panel A's two inputs (`coloc_gene`, `coloc_rate`) were restaged from
`submissions_757aa6ab-…_call-CollectHighLevelOutputs_high_level_coloc_outputs (1).tar`.
Panels B/C/D are untouched — their inputs did not change. To reproduce, extract
the two tables into a staging root and point `FIG4_INPUT` at it:

    tar -xf '…high_level_coloc_outputs (1).tar' -C "$STAGE/Fig4Data" \
      high_level_coloc_outputs/harmonized/harmonized_coloc.gene.tsv.gz \
      high_level_coloc_outputs/final_summaries/coloc_rate_by_trait.tsv
    rm -f Fig4/data/coloc_by_trait.rds   # the derived cache is keyed by NAME
    FIG4_INPUT=$STAGE Rscript fig4/run_all.R

**The staging root still needs fixing up properly** — this run used a throwaway
directory that symlinked the unchanged B/C/D inputs back to the old root.

What changed in the data, beyond the counts (monocyte 1010 → 1256, height
543 → 708, creatinine 216 → 354): every trait name is now lower-case, `autoimmune
disease`, `systemic lupus erythematosus` and `multiple sclerosis` are new,
`gastroesophageal reflux disease` drops below the `count > 5` cutoff, and the
duplicated `Asthma`/`asthma` and `body height`/`height` entries are gone. The
lower-casing broke every `TRAIT_SHORTHAND` lookup — see the labelling section
below.

## Appearance-affecting change: Panel C's axis titles

`plot_qq()`'s axes now read "Expected TWAS −log₁₀(P)" / "Observed TWAS
−log₁₀(P)", naming the test the *P*-values come from. `LEGEND.md` is updated to
match. Labels only; nothing plotted changes.

## Appearance-affecting change: Panel A's bold column headers

`build_plot()` sets `face = "bold"` on the **column** strips only ("Colocalizing
Genes" / "CS coloc rate "). The trait-category strips down the side are left at
normal weight: they are already set apart by size and by the grey rule beside
them, and bolding those too made the panel read as all heading. `assemble.R`
re-specifies `strip.text.y.left` to resize it to 13pt, so that override has to
stay unbolded as well or it puts the bold back. No numbers change; the plotted
data is identical.

## New: the labelled Panel A variant, and a second output file

`build_plot(label_top = n)` annotates the n highest-count traits per category
with a short label, to use the dead space to the right of the bars; `Inf` labels
every trait, which is what the variant uses. `label_top = 0` is the default and
reproduces the previous panel exactly — no labels are drawn, nothing extends the
axis and the breaks are ggplot's defaults, so `Figure4.pdf` is unaffected.

In the labelled variant the count axis runs past the data to make room for the
label column, and its **breaks are restricted to the data range**. Letting them
run to the extended limit both implies values that do not exist out there and put
a "1500" hard against the rate facet's "0%" across the 8pt panel gap, which read
as "15000%". Its x-axis breaks therefore differ from the plain panel's.

All the labels are right-aligned at one x, in a column reserved by **extending
the count axis past the longest bar** (to ~1500 with the current data, whose
longest bar is 1256). Every label therefore shares an x and none of them touches
a bar.

Extending the axis is what makes the alignment possible, and the earlier attempts
are worth recording because they both fail on real data. Drawing each label just
past the end of its own bar destroys the column. Keeping the column at the data
maximum and recolouring the labels that bars reach — white on the bar, grey off
it — fails whenever a bar stops *part-way through* a label: the label then
straddles bar and background and no single colour is legible on both. That is
exactly what `eosinophil` (1108 of 1256) did; anchoring it on its own bar end
made it legible but left it out of line with every other label. With the column
reserved there is no straddle case, so the white-text rule is gone and every
label is plain grey.

How much room the column needs is estimated per label from its character count,
since a 15-character label needs far more than a 3-character one:
`label_char_frac` (0.0165) is the width of one character and `label_gap` (0.03)
the clearance between the longest bar and the column, both as fractions of the
count-axis range. **They are estimates of how wide the text renders**, so they
need retuning if `label_size` or the panel geometry changes; the symptom is a
label touching its bar.

`panels/panelA_labeled.R` builds it to `objects/A_labeled.rds`, and `assemble.R`
now writes **two** files from one layout: `Figure4.pdf`/`.png` from the plain
panel and `Figure4_labeled.pdf`/`.png` from the labelled one. Both are covered by
`.gitignore`.

Label text comes from `trait_label()`: the `TRAIT_SHORTHAND` override if the
trait has one, else the upstream name tidied (underscores to spaces; first letter
capitalised only when the name is entirely lower-case). `TRAIT_SHORTHAND` is an
override map for the long, clinical and acronym names, **not** a whitelist, so a
new upstream trait still gets a sensible label.

**`TRAIT_SHORTHAND` is keyed in lower case and matched case-insensitively**, and
that is load-bearing. The coloc table changed the case of its trait names between
releases (`Asthma` → `asthma`, `AST` → `ast`); with the exact-match lookup the
map originally used, every override silently stopped matching, and the acronyms
then sentence-cased to `Ast`, `Bun`, `Mchc`, `Rdw`, `Hba1c`. Hence the explicit
lower-case acronym entries and the two regression tests in
`tests/test_fig4_helpers.R`.

`build_plot()` also reports, via `message()`, any label that appears twice within
one category. The previous coloc release needed it: it carried `Asthma` and
`asthma` as separate `allergic_respiratory` traits, which collapsed to the same
text on two bars. The current release has a single `asthma`, so nothing is
reported — but deduplicating near-duplicate trait strings is an upstream data
decision, not a plotting one, so the check stays.

`Fig4.ipynb` gained a cell for the labelled variant and was regenerated and
re-executed with the panels rendered. That cell renders taller than the other
panels (`repr.plot.height = 9.5`, restored afterwards): the panel is 70 rows, and
at the notebook's default 7.5 × 5 in preview the rows are shorter than the label
text, so the labels collide there even though they do not in the 15 × 10 figure.

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

## Open issue: Panel C's "Null" series falls off the colour scale

`panels/panelC.R` renames the matched control group to `"Null"`
(`case_when(str_detect(group, "Matched") ~ "Null", ...)`) *before* calling
`plot_qq()`, whose `scale_colour_manual` only maps
`"Matched non-colocalized genes"` → `#2B2B2B` and `"Colocalized genes"` →
`#275D8E`. `"Null"` therefore matches no level and renders in ggplot2's default
grey `na.value` instead of the intended `#2B2B2B`. The legend is switched off in
this panel, so the mismatch is not otherwise visible. Left as-is for now —
fixing it changes the rendered colour of the null series.

## Hand-curated heme/immune keyword list

The heme/immune gene set is defined by a regex of disease-name keywords
(`HEME_IMMUNE_KW` in `00_setup.R`) matched against ClinGen `diseaseFromSource`.
It is a curated list, not an ontology mapping.

## Local input requirements

`inputs.tsv` records these local-only. The two `high_level_coloc_outputs/`
tables are extracted from the private `high_level_coloc_outputs.tar.gz` archive
(README documents extraction); `R2_TWAS_combined.tsv`,
`ConditionalTWAS.cleaned.tsv.gz`, `reference/interactions.tsv` and the ClinGen
parquet must be supplied at their listed local paths.
