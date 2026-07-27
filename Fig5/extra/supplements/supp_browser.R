## ============================================================================
## Supplement Borzoi browsers for the huge-effect / FAT1 candidates.
## Reuses the SPR browser machinery (fig5/00_setup.R) generalized per gene.
## Run from the PROJECT ROOT:  Rscript fig5/extra/supplements/supp_browser.R
## ============================================================================
source("fig5/extra/supplements/_setup.R")   # fonts, DIR_BORZOI/DIR_OUT + folds in SPR machinery

EXPORTS <- file.path(DIR_BORZOI, "exports")
GENES <- c("FAT1","SLC1A2","SLC6A19","SCN11A","BSND")

make_browser <- function(gene) {
  bw_dir <- file.path(EXPORTS, paste0(gene, "_bigwig"))
  meta   <- fread(file.path(bw_dir, "meta.tsv"))
  vpos   <- meta$pos[1]
  md     <- .spr_build_metadata(bw_dir)
  tracks <- .spr_load_all(md, variant_pos = vpos)
  ## generic element highlight: 200 bp around the variant
  br <- tryCatch(plot_browser_grouped(tracks, element = c(-100, 100), min_pct = 8),
                 error = function(e) plot_browser_grouped(tracks, element = c(-100, 100), min_pct = 0))
  title <- sprintf("%s  %s:%s %s>%s   (aFC in caption)", gene, meta$chrom[1],
                   format(vpos, big.mark=","), meta$ref[1], meta$alt[1])
  cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label(title, x=0.01, hjust=0, size=12, fontface="bold"),
                     br, ncol=1, rel_heights=c(0.04, 1))
}

for (g in GENES) {
  bw_dir <- file.path(EXPORTS, paste0(g, "_bigwig"))
  if (!file.exists(file.path(bw_dir, "meta.tsv"))) { message("skip ", g, " (no bigwigs yet)"); next }
  p <- tryCatch(make_browser(g), error=function(e){message(g, " FAILED: ", conditionMessage(e)); NULL})
  if (!is.null(p)) {
    ggsave(file.path(DIR_OUT, sprintf("supp_%s.pdf", g)), p, width=9, height=7, bg="white", limitsize=FALSE)
    message("wrote supp_", g, ".pdf")
  }
}
