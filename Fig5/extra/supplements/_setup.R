# =============================================================================
# Shared setup for the Figure 5 supplements (Fig5/supplements/*.R).
# Self-contained; run each script from this directory (Fig5/supplements/).
# Provides fonts, the ANC palette, borzoi/AG/output paths, and the SPR browser
# machinery (Fig5/panels/_helpers_spr.R). Override paths via env vars.
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(cowplot); library(ggseqlogo)
  library(patchwork); library(rcartocolor); library(ggrepel); library(scales)
  library(showtext); library(sysfonts)
})

setup_fonts <- function() {
  if (!"inter" %in% sysfonts::font_families())
    try(sysfonts::font_add_google("Inter", "inter"), silent = TRUE)
  showtext::showtext_opts(dpi = 400); showtext::showtext_auto(); invisible(TRUE)
}
setup_fonts()

# AlphaGenome ancestry palette (Okabe-Ito; common = grey)
ANC <- c(AFR = "#E69F00", EAS = "#56B4E9", EUR = "#009E73",
         SAS = "#CC79A7", common = "#999999")

DIR_BORZOI <- Sys.getenv("FIG5_BORZOI", "/Users/evinmpadhi/borzoi_run")
DIR_AG     <- Sys.getenv("FIG5_AG",     "/Users/evinmpadhi/Fig5_QTL_paper/alphagenome")
DIR_OUT    <- Sys.getenv("FIG5_SUPP_OUT", "output")
if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
PATH_SUSIE <- Sys.getenv("FIG5_SUSIE",
  "/Users/evinmpadhi/Desktop/susie_files/COMB.expression_SusieMerged.annotated.ancestry_skew.tsv.gz")

# The SPR / Borzoi browser machinery (plot_browser_grouped, .spr_*, motif_logo)
# now lives in the folded fig5/00_setup.R. Supplement scripts that need it
# (supp_browser.R, plot_ADA_supp.R) are run from the PROJECT ROOT, where this
# sources it; scripts run from a nested dir (ada_scifig/) skip it — they only
# build ADA panels and don't use the browser.
if (file.exists("fig5/00_setup.R")) source("fig5/00_setup.R")
