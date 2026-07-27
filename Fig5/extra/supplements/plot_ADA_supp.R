## ============================================================================
## Supplementary Figure — ADA enhancer allelic series (Borzoi deep dive)
##   A  fine-mapping of the independent rare variants in the ADA enhancer
##   B  predicted immune chromatin accessibility (ref vs alt) for the ETS pair
##   C  base-resolution ETS-motif attribution (created vs destroyed)
##   D  predicted immune-TF binding change across the allelic series
## Styling matched to Figure 5.
## ============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(rtracklayer)
  library(ggseqlogo); library(patchwork); library(cowplot)
  library(rcartocolor); library(ggrepel)
})

# Run from the PROJECT ROOT:  Rscript fig5/extra/supplements/plot_ADA_supp.R
source("fig5/extra/supplements/_setup.R")   # fonts, ANC, DIR_*, + folds in SPR browser helpers
theme_set(theme_classic())
theme_update(strip.text = element_text(size = 11), text = element_text(family = "inter"),
             axis.title = element_text(size = 12), axis.text = element_text(size = 10),
             strip.background = element_blank(), plot.tag = element_text(size = 20, face = "bold"),
             plot.title = element_text(size = 12.5, face = "bold"),
             plot.margin = margin(6, 10, 6, 10))

EXP <- file.path(DIR_BORZOI, "exports")
A   <- file.path(DIR_BORZOI, "results/attributions")
GAIN <- "#2A7FB8"; LOSS <- "#C0392B"; REFCOL <- "#2A7FB8"; ALTCOL <- "#C0392B"

FOCAL <- tribble(
  ~pos,      ~ref,~alt, ~pop,   ~aFC,
  44643792,  "C","T", "common", -0.32,
  44646218,  "C","T", "SAS",     0.96,
  44646238,  "C","T", "EAS",    -1.00,
  44646509,  "T","C", "EUR",     0.65,
  44650274,  "C","A", "AFR",    -0.68) %>%
  mutate(tag = as.character(pos), dir = if_else(aFC > 0, "increases", "decreases"))

## ---------------------------------------------------------------------------
## A — ADA fine-mapping (enhancer region), focal variants coloured by direction
## ---------------------------------------------------------------------------
su  <- fread(PATH_SUSIE)
ada <- su %>% filter(gene_name == "ADA", between(position, 44643000, 44651000)) %>%
  mutate(focal = position %in% FOCAL$pos) %>%
  left_join(FOCAL %>% select(pos, dir, pop, aFC), by = c("position" = "pos"))

fm <- ggplot(ada, aes(position / 1e6, pip)) +
  annotate("rect", xmin = 44646.15/1e3, xmax = 44646.6/1e3, ymin = -Inf, ymax = Inf,
           fill = "#F4A261", alpha = 0.14) +
  geom_hline(yintercept = 0.9, linetype = "dashed", linewidth = 0.4) +
  geom_point(data = ~filter(.x, !focal), colour = "grey72", size = 1.8, alpha = 0.7) +
  geom_point(data = ~filter(.x, focal), aes(fill = dir), shape = 23, size = 3.4,
             stroke = 0.5, colour = "black") +
  geom_text_repel(data = ~filter(.x, focal),
                  aes(label = sprintf("%s>%s (%s, %+.2f)", ref, alt, pop, aFC)),
                  size = 3, family = "inter", min.segment.length = 0, box.padding = 0.5,
                  segment.size = 0.3, max.overlaps = Inf, seed = 1) +
  scale_fill_manual(values = c(increases = GAIN, decreases = LOSS),
                    name = "eQTL effect", labels = c(increases = "increases ADA", decreases = "decreases ADA")) +
  scale_y_continuous(limits = c(0, 1.08), breaks = c(0, 0.5, 0.9)) +
  labs(x = "chr20 position (Mb)", y = "PIP", tag = "A",
       title = "ADA fine-mapping — allelic series at the −6 kb lymphoid enhancer") +
  theme(legend.position = c(0.5, 0.42), legend.direction = "horizontal",
        legend.title = element_text(size = 9), legend.text = element_text(size = 9),
        plot.tag.position = c(0.005, 1.02))

## ---------------------------------------------------------------------------
## B — predicted immune accessibility (ref vs alt) for the ETS pair
## ---------------------------------------------------------------------------
load_dnase <- function(tag) {
  dir <- file.path(EXP, paste0("ADA_", tag, "_bigwig")); v <- fread(file.path(dir, "meta.tsv"))$pos[1]
  bind_rows(as.data.frame(import(file.path(dir, "DNase.ref.bw"))) %>% transmute(rel = start - v, score, allele = "ref"),
            as.data.frame(import(file.path(dir, "DNase.alt.bw"))) %>% transmute(rel = start - v, score, allele = "alt")) %>%
    filter(abs(rel) <= 800)
}
acc_panel <- function(tag, ttl, tag_lab = NULL, show_leg = FALSE) {
  ggplot(load_dnase(tag), aes(rel, score)) +
    annotate("rect", xmin = -60, xmax = 60, ymin = -Inf, ymax = Inf, fill = "#F4A261", alpha = 0.22) +
    geom_area(data = ~filter(.x, allele == "ref"), fill = REFCOL, alpha = 0.22) +
    geom_line(aes(colour = allele), linewidth = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey35") +
    scale_colour_manual(values = c(ref = REFCOL, alt = ALTCOL), name = NULL,
                        labels = c(ref = "reference", alt = "alternate")) +
    scale_x_continuous(breaks = seq(-800, 800, 400)) +
    labs(x = "position relative to variant (bp)", y = "predicted immune DNase", title = ttl, tag = tag_lab) +
    theme(plot.title = element_text(size = 10.5, face = "bold"),
          legend.position = if (show_leg) c(0.78, 0.9) else "none",
          legend.text = element_text(size = 8), legend.key.size = unit(3.5, "mm"),
          plot.tag.position = c(0.005, 1.02))
}
b1 <- acc_panel("44646218", "44,646,218 C>T (SAS, +0.96):  accessibility +34%", "B", show_leg = TRUE)
b2 <- acc_panel("44646238", "44,646,238 C>T (EAS, −1.00):  accessibility −50%")

## ---------------------------------------------------------------------------
## C — ETS-motif attribution logos (created vs destroyed)
## ---------------------------------------------------------------------------
ets_logo <- function(tag, allele, ttl, core_rel, tag_lab = NULL) {
  a <- fread(file.path(A, paste0("ada_ets_", tag, ".tsv"))) %>% filter(rel >= -7, rel <= 9)
  bcol <- if (allele == "alt") "alt_base" else "ref_base"; acol <- if (allele == "alt") "attr_alt" else "attr_ref"
  mat <- matrix(0, 4, nrow(a), dimnames = list(c("A", "C", "G", "T"), as.character(a$rel)))
  for (i in seq_len(nrow(a))) mat[a[[bcol]][i], i] <- a[[acol]][i] / 1000
  vcol <- which(a$is_var == 1); cc <- which(a$rel %in% core_rel)
  ggseqlogo(mat, method = "custom") +
    annotate("rect", xmin = min(cc) - 0.5, xmax = max(cc) + 0.5, ymin = -Inf, ymax = Inf, fill = "#F4A261", alpha = 0.28) +
    annotate("rect", xmin = vcol - 0.5, xmax = vcol + 0.5, ymin = -Inf, ymax = Inf, fill = NA, colour = LOSS, linewidth = 0.6) +
    annotate("text", x = mean(cc), y = Inf, label = "ETS core", vjust = 1.3, size = 2.7, family = "inter", colour = "grey25") +
    scale_x_continuous(breaks = seq(1, nrow(a), 4), labels = a$rel[seq(1, nrow(a), 4)]) +
    labs(x = "position relative to variant (bp)", y = "attribution", title = ttl, tag = tag_lab) +
    theme(plot.title = element_text(size = 10.5, face = "bold"),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5), plot.tag.position = c(0.005, 1.02))
}
c1 <- ets_logo("44646218", "alt", "ETS core GAINED (creates TTCC)", -1:2, "C")
c2 <- ets_logo("44646238", "ref", "ETS core LOST (disrupts TTCC)", 1:4)

## ---------------------------------------------------------------------------
## D — predicted immune-TF binding change across the allelic series
## ---------------------------------------------------------------------------
vlab <- c(`44646218` = "44,646,218\nSAS +0.96", `44646238` = "44,646,238\nEAS −1.00",
          `44646509` = "44,646,509\nEUR +0.65", `44650274` = "44,650,274\nAFR −0.68",
          `44643792` = "44,643,792\ncommon −0.32")
alld <- imap_dfr(vlab, ~ fread(file.path(EXP, paste0("ADA_", .y, "_bigwig"), "deltas.tsv")) %>%
                   filter(kind == "TF") %>% mutate(v = .x))
top <- alld %>% group_by(mark) %>% summarise(mx = max(abs(pct))) %>% arrange(desc(mx)) %>% head(14) %>% pull(mark)
hm <- alld %>% filter(mark %in% top) %>%
  mutate(mark = factor(mark, levels = rev(top)), v = factor(v, levels = vlab), pct = pmax(pmin(pct, 60), -60)) %>%
  ggplot(aes(v, mark, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.0f", pct)), size = 2.5, family = "inter") +
  scale_fill_gradient2(low = LOSS, mid = "white", high = GAIN, midpoint = 0, limits = c(-60, 60),
                       name = "TF binding\nΔ (%)") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, tag = "D", title = "Predicted immune-TF binding change (alt vs ref)") +
  theme(axis.text.x = element_text(size = 8.5), axis.text.y = element_text(size = 9),
        axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
        legend.title = element_text(size = 9), plot.tag.position = c(0.005, 1.02))

## ---------------------------------------------------------------------------
figure <- fm / (b1 | b2) / (c1 | c2) / hm + plot_layout(heights = c(0.72, 1, 0.8, 1.75))
ggsave(file.path(DIR_OUT, "FigureS_ADA.pdf"), figure,
       width = 11, height = 14.5, bg = "white", limitsize = FALSE)
message("wrote FigureS_ADA.pdf")
