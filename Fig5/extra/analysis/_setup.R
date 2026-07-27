# =============================================================================
# Shared setup for the AlphaGenome analysis panels (Fig5/analysis/*.R).
# Self-contained (these read AlphaGenome result parquets directly, not via the
# figure manifest). Source at the top of each plot_*.R:  source("_setup.R")
# — run each script from this directory (Fig5/analysis/).
# Provides: AG (results dir), FIGDIR (output), ANC palette, compact theme.
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(arrow); library(showtext); library(sysfonts)
})
if (!"inter" %in% sysfonts::font_families())
  try(font_add_google("Inter", "inter"), silent = TRUE)
showtext_opts(dpi = 400); showtext_auto()

theme_set(theme_classic(base_size = 8, base_family = "inter") +
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "bold", hjust = 0),
        axis.text   = element_text(colour = "grey20"),
        plot.margin = margin(4, 6, 4, 6)))

# AlphaGenome ancestry palette (Okabe-Ito; common = grey)
ANC <- c(AFR = "#E69F00", EAS = "#56B4E9", EUR = "#009E73",
         SAS = "#CC79A7", common = "#999999")

AG     <- Sys.getenv("FIG5_AG", "/Users/evinmpadhi/Fig5_QTL_paper/alphagenome")  # results (parquets)
FIGDIR <- Sys.getenv("FIG5_AG_FIG", "figures")                                   # figures written here
if (!dir.exists(FIGDIR)) dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
