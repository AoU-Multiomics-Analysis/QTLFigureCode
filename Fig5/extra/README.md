# Fig5 `extra/` — standalone analyses & supplements

These are retained provenance and supplementary analyses. **None are run by
`fig5/run_all.R`** and none are part of the assembled Figure 5. All read
controlled-access or large model-derived inputs, so their rendered outputs are
`.gitignore`d and never committed.

- **`analysis/`** — AlphaGenome cell-type deep-dives (CAGE / DNase / RNA
  deconvolution, disease-context screens). Self-contained: each `plot_*.R`
  sources the local `analysis/_setup.R` and reads AlphaGenome result parquets
  directly (via `FIG5_AG`, default `~/Fig5_QTL_paper/alphagenome`). Run from
  inside `fig5/extra/analysis/`. Figures are written to `analysis/figures/`.

- **`supplements/`** — Borzoi/AlphaGenome supplementary figures (the ADA
  allelic-series supplement and per-gene browsers). `supp_browser.R` and
  `plot_ADA_supp.R` reuse the SPR browser machinery folded into
  `fig5/00_setup.R`, so **run them from the project root**
  (`Rscript fig5/extra/supplements/supp_browser.R`); the local `_setup.R` sources
  the folded setup when run from root. The `ada_scifig/` sub-figure runs from its
  own directory. Outputs go to `supplements/output/` (override with
  `FIG5_SUPP_OUT`).

- **`alphagenome/`** — the Python AlphaGenome scoring scripts that generated the
  SPR / ACVRL1 model inputs (bigWigs, attribution tables, screen parquets). See
  `alphagenome/README.md`.
