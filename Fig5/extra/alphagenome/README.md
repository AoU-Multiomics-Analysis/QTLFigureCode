# AlphaGenome scoring — provenance

The Python scripts here generated the sequence-model inputs used by Figure 5
(SPR / ACVRL1 Borzoi ref/alt bigWigs, CTCF attribution tables, and the
rare-variant bridge screens). They are **provenance only** — not part of the R
figure pipeline and not run by `fig5/run_all.R`.

**Requirements**
- A Python environment with the `alphagenome` client installed (the project used
  a dedicated venv; see the AlphaGenome docs).
- An AlphaGenome API key exported as `AG_KEY`:

  ```sh
  export AG_KEY=...            # never hard-code the key
  ```

`config.require_key()` reads `AG_KEY` from the environment and errors if it is
unset. Output locations default to `FIG5_AG` (`~/Fig5_QTL_paper/alphagenome`) and
can be overridden by env var. No keys or participant-derived data are committed.
