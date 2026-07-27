# Fig5 cleanup — implementation plan

*2026-07-27. Executes `docs/superpowers/specs/2026-07-27-fig5-cleanup-design.md`.*

## Global constraints

- **Run everything from the project root** `/Users/evinmpadhi/QTLFigureCode`.
- Scripts self-reference **lowercase** `fig5/…` (dir is `Fig5/`; case-insensitive macOS).
- Single `DPI <- 400` feeds `showtext_opts()` + every `ggsave()`. `element_text` sizes in
  **pt**; `geom/annotate/geom_label` in **mm**.
- **Tagging only in `assemble.R`** (explicit `labs(tag=)`); panels return untagged objects.
- **Never `git add -A`** — controlled-access files sit untracked in the tree; stage named
  paths only. Never commit `susie`, `nohom`, `SPR_bigwig`, `data/`, `objects/`, or outputs.
- No hard-coded AlphaGenome key anywhere; provenance scripts read env `AG_KEY`.
- Faithful source for verbatim bodies: the **validated ad-hoc refactor** at
  `$REF = /private/tmp/claude-501/-Users-evinmpadhi-Desktop-susie-files/122327f8-4348-4577-bd1a-0eefdd0625ec/scratchpad/fig5_adhoc_reference`
  (builds the exact `Figure5.pdf`, 1,136,396 B). Original flat scripts: `~/Fig5_QTL_paper/`.
- Each task ends in its own commit (no amends). Run `Rscript tests/run_all.R` at tasks 3,4,6,7,10.

## File structure

| Path | Responsibility | Action |
|---|---|---|
| `Fig5/_original/{PlotFigure5,acvrl1_panels,spr_browser_panel}.R` | recoverable original | Commit (T1) → `git rm` (T10) |
| `Fig5/.gitignore` `README.md` `NOTES.md` | scaffolding/docs | Create |
| `Fig5/inputs.tsv` | manifest (12 keys) | Create |
| `Fig5/reference/*.{tsv,csv}` | committed small tables | Create |
| `Fig5/00_setup.R` | libs, paths, manifest, font/DPI, theme, palettes, utils, build_or_cache, folded loaders + A–C/ACVRL1/SPR helpers | Create |
| `Fig5/panels/panel{A..G}.R` | one untagged object each → `objects/X.rds` | Create |
| `Fig5/assemble.R` `run_all.R` | tag + compose; run | Create |
| `Fig5/extra/*.R` + `Fig5/extra/alphagenome/*` | standalone analyses + provenance | Create |
| `tests/fixtures.R` | `make_fig5_fixture()` | Modify |
| `tests/test_fig5_helpers.R` `tests/test_fig5_assembly.R` | pure-seam + rendered-tag | Create |
| `tools/build_notebooks.py` | add Fig5 config + PANEL_NOTES | Modify |
| `Fig5/Fig5.ipynb` | generated notebook | Create |

---

## Task 1 — Scaffold + preserve original

- [ ] `mkdir -p Fig5/{panels,reference,extra,_original}`.
- [ ] Copy `~/Fig5_QTL_paper/{PlotFigure5.R,acvrl1_panels.R,spr_browser_panel.R}` →
      `Fig5/_original/`. **Scrub**: `grep -rnE "AIza|api_key|API_KEY" Fig5/_original` must be empty.
- [ ] `Fig5/.gitignore`:
      ```
      objects/
      data/
      inputs/
      ```
      (plus `.DS_Store`). `inputs/` holds staged controlled-access files → ignored.
- [ ] Stub `Fig5/README.md`, `Fig5/NOTES.md` (filled in T10).
- [ ] Commit: `Fig5 cleanup: preserve original flat scripts + scaffolding`.

## Task 2 — 00_setup.R core + inputs.tsv + reference tables

- [ ] Copy `$REF/reference/*` → `Fig5/reference/` (6 files: `clingen_gene_validity.csv`,
      `acvrl1_promoterAI.tsv`, `acvrl1_ccre_region.tsv`, `spr_ccre_region.tsv`,
      `gene_exons.tsv`, `gene_spans.tsv`). Verify none are participant-derived (they are
      public/ENCODE/GENCODE/derived-annotation) → safe to commit.
- [ ] `Fig5/inputs.tsv` = `$REF/inputs.tsv` verbatim (12 keys; header `key	local_subpath	gs_uri`).
- [ ] `Fig5/00_setup.R` **head** (through `build_or_cache`) = `$REF/00_setup.R:1–122` with the
      single change **`"Fig5/…"` → `"fig5/…"`** in the five `Sys.getenv` path defaults and the
      three `source("Fig5/panels/_helpers_*.R")` lines removed (helpers fold in below, T3–T6).
      Keep: libs (rtracklayer namespaced comment), env paths, `source("common/manifest.R")`,
      Inter/sans font, `DPI <- 400`, `SZ_*`, theme, palettes, `signed_square_trans`, `fmt_af`,
      `make_freq_label`, `normalize_variant_id`, `make_variant_id`, `build_or_cache`.
- [ ] Parse-check: `Rscript -e 'source("fig5/00_setup.R")'` (will error only if a later
      helper is referenced — none yet). Commit: `Fig5 cleanup: 00_setup core + manifest + reference tables`.

## Task 3 — Folded loaders + pure seams + fixture + first tests (TDD)

- [ ] **Refactor loaders into pure seams** when folding `$REF/00_setup.R:124–178` into
      `00_setup.R`:
  - `derive_clingen(clingen_raw, susie_gene_map)` — the `filter/distinct/transmute/left_join/rename`
    body of `load_reference`'s clingen block, returning the clingen frame. `load_reference()`
    calls it after the two `fread`/`read_excel` reads.
  - `annotate_eqtl(susie, genebayes, clingen_genes)` — the `mutate/left_join/mutate` body of
    `load_eqtl_susie`. `load_eqtl_susie(reference)` extracts locals then calls it.
  - `build_variant_freq(eqtl)` folds verbatim (already pure).
  - keep `load_nohomozygotes`, `get_reference/get_eqtl/get_variant_freq` getters verbatim.
- [ ] **`tests/fixtures.R`**: add `make_fig5_fixture(root=…)` writing mock `susie.tsv`
      (cols `gene_id, gene_name, variant, gvs_max_af, gvs_max_subpop, log2_aFC, posterior_mean,
      pip, ref, clinvar_classification, position, cs_id, type`; ~8 rows spanning
      ClinGen/control, rare/common, AD/AR, a gene present as **both AD and AR**), mock
      `clingen.csv` (6 filler lines then header-less rows matching `fread(skip=6, select=1:7)`),
      mock `genebayes` as an `.rds` with cols `ensg, post_mean`, and a small SPR-track tibble
      for `pct_at_element`. Return `list(root, susie, clingen, genebayes, spr_tracks)`.
- [ ] **`tests/test_fig5_helpers.R`** (test-first — write, watch fail, then the fold makes it
      pass): source fixtures + `fig5/00_setup.R` + `tests/helpers.R`; `ok()` on:
  `derive_clingen` (keeps Def/Strong AD/AR/SD, drops others), `annotate_eqtl`
  (`ClinGen`/`rare`/`post_mean` correct), `build_variant_freq` (label formats),
  `make_variant_id` (delim normalization).
- [ ] `Rscript tests/run_all.R`. Commit: `Fig5 cleanup: folded loaders as pure seams + fixture + helper tests`.

## Task 4 — A–C helpers folded + panel scripts + fan-out test

- [ ] Fold `$REF/panels/_helpers_abc.R` into `00_setup.R` (section "Panels A–C"),
      **removing every `tag = "A"/"B"/"C"` and `labs(..., tag=…)`** from `build_panelA/B/C`
      (tagging moves to assemble). Keep `variant_effect_by_inheritance` verbatim (preserves
      the AD+AR fan-out per spec §2).
- [ ] Extend `test_fig5_helpers.R`: `variant_effect_by_inheritance` on the fixture →
      assert Control/AD/AR assignment **and** that the AD+AR gene yields a row in *both* AD
      and AR (pins the preserved double-count).
- [ ] `panels/panelA.R`, `panelB.R`, `panelC.R`: each `source("fig5/00_setup.R")`, build via
      `build_panel{A,B,C}(get_eqtl(), get_reference()$clingen[, nohom for C])`,
      `saveRDS(obj, file.path(OBJ_DIR, "X.rds"))`, `message("Panel X built.")`.
- [ ] `Rscript tests/run_all.R`. Commit: `Fig5 cleanup: panels A-C (untagged) + fan-out pinned`.

## Task 5 — Fine-mapping helpers + ACVRL1 panels D/E

- [ ] Fold the fine-mapping block `$REF/00_setup.R:180–337` (`.gene_exons_all`, `.gene_spans_all`,
      `ccre_bands`, `gene_model_track`, `plot_finemapping_with_inset`) into `00_setup.R`
      **removing `tag = panel_tag`** from `plot_finemapping_with_inset` (drop the `panel_tag`
      arg + its `labs(tag=)`; keep `panel_title`). Fold `$REF/panels/_helpers_acvrl1.R`
      (ACVRL1 fine-map + PromoterAI builders), removing self-tags.
- [ ] `panels/panelD.R` (ACVRL1 fine-mapping patchwork), `panelE.R` (PromoterAI), untagged,
      saveRDS + message.
- [ ] Parse-check `Rscript -e 'source("fig5/panels/panelD.R")'` needs inputs; instead
      `Rscript -e 'source("fig5/00_setup.R")'` clean. Commit: `Fig5 cleanup: fine-map helpers + ACVRL1 panels D/E`.

## Task 6 — SPR helpers + panels F/G + pct_at_element test

- [ ] Fold `$REF/panels/_helpers_spr.R` (SPR metadata/loaders, `pct_at_element`,
      `plot_browser_grouped`, `motif_logo`, `make_spr_browser`) into `00_setup.R`, self-tags
      removed. `rtracklayer::import` stays namespaced.
- [ ] Extend `test_fig5_helpers.R`: `pct_at_element(fixture$spr_tracks)` → % change sign/scale.
- [ ] `panels/panelF.R` (SPR fine-map patchwork), `panelG.R` (`make_spr_browser()`), untagged.
- [ ] `Rscript tests/run_all.R`. Commit: `Fig5 cleanup: SPR helpers + panels F/G + pct_at_element test`.

## Task 7 — assemble.R (explicit tags) + run_all.R + rendered-tag test

- [ ] `assemble.R`: read `objects/{A..G}.rds`; **apply tags explicitly** — ggplots B/C/E via
      `+ labs(tag=…)`; patchworks A/D/F via `+ plot_annotation(...) & theme(plot.tag=…)` scoped
      to the sub-object, or wrap + `labs`; cowplot G via `wrap_elements(P$G) + labs(tag="G")`.
      Compose the locked layout (`$REF/assemble.R:20–29`), `ggsave` pdf(DPI)+png(100 toggle).
- [ ] `run_all.R`: `for (p in LETTERS[1:7]) source(sprintf("fig5/panels/panel%s.R", p)); source("fig5/assemble.R")`.
- [ ] `tests/test_fig5_assembly.R`: reuse `extract_tag_labels`; build 7 lightweight stand-in
      objects (patchwork for A/D/F, ggplot for B/C/E, `ggdraw()` cowplot for G), source
      `assemble.R` against them, assert **exactly `{A,B,C,D,E,F,G}` (7 tags)** + pdf written.
- [ ] `Rscript tests/run_all.R`. Commit: `Fig5 cleanup: assemble with explicit tags + run_all + rendered-tag test`.

## Task 8 — extra/ analyses + AlphaGenome provenance

- [ ] Copy `$REF/analysis/*.R` + `$REF/supplements/*` → `Fig5/extra/`, adjusting `source()`
      paths to `fig5/00_setup.R` and self-tag/`panels/` references. Copy the AlphaGenome
      Python scripts → `Fig5/extra/alphagenome/`; **scrub keys** (`grep -rnE "AIza|api_key"`
      empty; ensure each reads `os.environ["AG_KEY"]`). Add `extra/alphagenome/README.md`
      noting the venv + `AG_KEY` requirement.
- [ ] Parse-check each R supplement sources cleanly. Commit: `Fig5 cleanup: extra analyses + AlphaGenome provenance`.

## Task 9 — Notebook

- [ ] Add `Fig5` to `tools/build_notebooks.py` `FIGURES` (`num=5, cap="Fig5", lower="fig5",
      obj_letters=list("ABCDEFG")`, controlled-access `data_note` + `NOTES.md` pointer, `env`
      FIG5_MANIFEST/OBJ/DATA, `steps` = 7 panels + assemble + `extra/*` supps) and `PANEL_NOTES`
      A–G. Setup cell regex `^(Fig5|fig5)$`.
- [ ] `python3 tools/build_notebooks.py Fig5`; confirm round-trip (regenerate → no diff).
      Commit: `Fig5 cleanup: notebook config + Fig5.ipynb`.

## Task 10 — Docs, remove original, full run, whole-branch review

- [ ] Fill `Fig5/README.md` (how to run, inputs + which need URIs, controlled-access note) and
      `Fig5/NOTES.md` (every result-affecting item from spec §7.1 with original line refs:
      preserved AD+AR fan-out, `ref`-column masking, Panel B `|posterior_mean|` choice).
- [ ] `git rm -r Fig5/_original`.
- [ ] Full pipeline from root: `Rscript fig5/run_all.R` → `Figure5.pdf` builds (compare size
      ≈1,136,396 B). `Rscript tests/run_all.R` all pass.
- [ ] Whole-branch review (most capable model): cross-file coherence, no stale artifacts, no
      committed controlled-access data (`git ls-files Fig5 | grep -vE '…'`), tags correct.
- [ ] Commit: `Fig5 cleanup: docs, remove preserved original, final verification`. Open PR.

## Self-review checklist (playbook step 3)

- [ ] Fixture `susie.tsv` columns cover **every** field the folded code reads (esp. `ref`,
      `clinvar_classification`, `post_mean` via join, `cs_id`, `type`, `position`) — a missing
      column silently NA-fills and the join "matches nothing".
- [ ] `derive_clingen`/`annotate_eqtl` signatures match both `load_*` callers and the tests.
- [ ] Removing `panel_tag` from `plot_finemapping_with_inset` doesn't strip `panel_title`.
- [ ] assemble tag logic verified by the **rendered** test, not by reading — patchwork A/D/F
      must contribute exactly one tag each (insets `ignore_tag`), G exactly one.
- [ ] No `Fig5/` (uppercase) self-ref survives: `grep -rn '"Fig5/' Fig5` empty.
