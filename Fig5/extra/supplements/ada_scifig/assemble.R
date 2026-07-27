source("build_panels.R")   # builds pA,pB,pC,pD (and writes individual SVGs); defines FONT, DIR_OUT
suppressPackageStartupMessages({library(cowplot)})
row1 <- plot_grid(pA, pB, nrow=1, rel_widths=c(1, 0.97), labels=c("a","b"),
                  label_size=11, label_fontfamily=FONT, label_fontface="bold")
figure <- plot_grid(row1, pC, pD, ncol=1, rel_heights=c(1, 0.86, 0.60),
                    labels=c("","c","d"), label_size=11, label_fontfamily=FONT, label_fontface="bold")
ggsave(file.path(DIR_OUT, "FigureADA.pdf"), figure, width=183, height=170, units="mm", device=cairo_pdf)
ggsave(file.path(DIR_OUT, "FigureADA.svg"), figure, width=183, height=170, units="mm", device=svglite)
message("wrote FigureADA.pdf / .svg (183 x 170 mm)")
