## ADA allelic series — publication panels (ggplot2 / tidyverse), SVG out for svgutils compose
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(svglite); library(systemfonts)
  library(ggseqlogo); library(patchwork); library(ggrepel)
})

source("../_setup.R")   # ANC palette, DIR_BORZOI, DIR_OUT (run from Fig5/supplements/ada_scifig/)

FONT <- "Helvetica"   # svglite/cairo_pdf render via systemfonts; Helvetica is guaranteed on macOS
EXP <- file.path(DIR_BORZOI, "exports"); ATT <- file.path(DIR_BORZOI, "results/attributions")
## Wong/Okabe-Ito ancestry colors; common = grey (shared ANC palette)

dat <- fread(file.path(ATT,"ada_attr_summary.tsv")) %>%
  mutate(ancestry = c("common","SAS","EAS","EUR","AFR"),
         dnase = c(-26.3, 34.5, -50.4, 42.8, 4.5),
         lab = sprintf("%s>%s", ref, alt),
         poskb = pos/1e3,
         miss = tag=="44650274") %>%
  mutate(ancestry = factor(ancestry, levels=names(ANC)))

base_th <- theme_classic(base_size=8, base_family=FONT) +
  theme(axis.text=element_text(size=7, colour="grey15"), axis.title=element_text(size=8, colour="black"),
        legend.text=element_text(size=7), legend.title=element_text(size=7.5),
        plot.title=element_text(size=8.5, face="bold"), legend.key.size=unit(3,"mm"),
        axis.line=element_line(linewidth=0.3), axis.ticks=element_line(linewidth=0.3))
dir.create("panels", showWarnings=FALSE, recursive=TRUE)   # intermediate per-panel SVGs
sv <- function(f,p,w,h) ggsave(file.path("panels",f), p, width=w, height=h, units="mm", device=svglite)

## ---- A: allelic series lollipop (genomic position x measured aFC) ----
pA <- ggplot(dat, aes(poskb, aFC)) +
  annotate("rect", xmin=44646.1, xmax=44646.6, ymin=-Inf, ymax=Inf, fill="#F4A261", alpha=0.15) +
  annotate("text", x=44646.35, y=-1.2, label="ADA lymphoid enhancer", size=2.2, colour="grey40", family=FONT) +
  geom_hline(yintercept=0, linewidth=0.3, colour="grey50") +
  geom_segment(aes(xend=poskb, yend=0, colour=ancestry), linewidth=0.5) +
  geom_point(aes(colour=ancestry), size=2.6) +
  geom_text_repel(aes(label=lab), size=2.3, family=FONT, box.padding=0.5, point.padding=0.3,
                  min.segment.length=0, seed=1, segment.size=0.2, direction="both") +
  scale_colour_manual(values=ANC, name="ancestry") +
  scale_y_continuous(limits=c(-1.3,1.28), breaks=seq(-1,1,0.5)) +
  labs(x="chr20 position (kb)", y=expression(eQTL~effect~(log[2]~aFC)),
       title="Rare-variant allelic series at the ADA enhancer") +
  base_th + theme(legend.position=c(0.81, 0.58), legend.background=element_blank(),
                  legend.key.height=unit(3.1,"mm"))
sv("A.svg", pA, 90, 64)

## ---- B: concordance — measured aFC vs predicted accessibility (attribution objective) ----
r_all <- cor(dat$aFC, dat$obj_delta); r_c <- cor(dat$aFC[!dat$miss], dat$obj_delta[!dat$miss])
pB <- ggplot(dat, aes(aFC, obj_delta)) +
  geom_hline(yintercept=0, linewidth=0.3, colour="grey70") + geom_vline(xintercept=0, linewidth=0.3, colour="grey70") +
  geom_smooth(data=~filter(.x,!miss), method="lm", se=FALSE, colour="grey55", linewidth=0.5, linetype="dashed") +
  geom_point(aes(colour=ancestry, shape=miss), size=2.8) +
  geom_text_repel(aes(label=paste0(lab," (",ancestry,")")), size=2.2, family=FONT, box.padding=0.35, min.segment.length=0, seed=2, segment.size=0.2) +
  scale_colour_manual(values=ANC, guide="none") +
  scale_shape_manual(values=c(`FALSE`=16, `TRUE`=1), guide="none") +
  annotate("text", x=-0.95, y=170, hjust=0, size=2.4, family=FONT, colour="grey25",
           label=sprintf("concordant (n=4): r = %.2f\nall (n=5): r = %.2f", r_c, r_all)) +
  annotate("text", x=-0.62, y=35, hjust=0, size=2.1, family=FONT, colour="#B8860B", label="Borzoi\nnon-concordant", lineheight=0.9) +
  labs(x=expression(measured~eQTL~effect~(log[2]~aFC)), y="predicted accessibility\nattribution Δ (alt − ref)",
       title="Predicted accessibility tracks the eQTL effect") + base_th
sv("B.svg", pB, 88, 64)

## ---- C: base-resolution accessibility attribution (ChromBPNet-style) — ETS pair ----
logo <- function(tag, allele, ttl, core_rel, showx=TRUE){
  a <- fread(file.path(ATT,paste0("ada_ets_",tag,".tsv"))) %>% filter(rel>=-6, rel<=8)
  bcol <- if(allele=="alt")"alt_base" else "ref_base"; acol <- if(allele=="alt")"attr_alt" else "attr_ref"
  m <- matrix(0,4,nrow(a),dimnames=list(c("A","C","G","T"),as.character(a$rel)))
  for(i in seq_len(nrow(a))) m[a[[bcol]][i],i] <- a[[acol]][i]/1000
  vcol <- which(a$is_var==1); cc <- which(a$rel %in% core_rel)
  ggseqlogo(m, method="custom") +
    annotate("rect", xmin=min(cc)-0.5, xmax=max(cc)+0.5, ymin=-Inf, ymax=Inf, fill="#F4A261", alpha=0.25) +
    annotate("rect", xmin=vcol-0.5, xmax=vcol+0.5, ymin=-Inf, ymax=Inf, fill=NA, colour="#D55E00", linewidth=0.5) +
    scale_x_continuous(breaks=seq(1,nrow(a),3), labels=a$rel[seq(1,nrow(a),3)]) +
    labs(title=ttl, x=if(showx)"position rel. variant (bp)" else NULL, y=allele) +
    theme_classic(base_size=8, base_family=FONT) +
    theme(plot.title=element_text(size=7.6,face="bold"), axis.text.x=element_text(size=6.4),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(),
          axis.title.y=element_text(size=7, angle=0, vjust=0.5), axis.title.x=element_text(size=7),
          axis.line.x=element_line(linewidth=0.3), plot.margin=margin(2,3,1,2))
}
c218 <- (logo("44646218","ref","44,646,218 C>T (SAS, +0.96): ETS core created", -1:2, FALSE) /
         logo("44646218","alt","", -1:2, TRUE))
c238 <- (logo("44646238","ref","44,646,238 C>T (EAS, -1.00): ETS core disrupted", 1:4, FALSE) /
         logo("44646238","alt","", 1:4, TRUE))
pC <- c218 | c238
sv("C.svg", pC, 180, 52)

## ---- D: per-variant predicted accessibility attribution (ranked diverging bars) ----
pD <- dat %>% mutate(v=fct_reorder(paste0(lab,"\n",ancestry), obj_delta)) %>%
  ggplot(aes(obj_delta, v, fill=ancestry)) +
  geom_col(width=0.68) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey40") +
  geom_text(aes(label=sprintf("aFC %+.2f", aFC), hjust=ifelse(obj_delta>0,-0.1,1.1)), size=2.2, family=FONT, colour="grey20") +
  scale_fill_manual(values=ANC, guide="none") +
  scale_x_continuous(limits=c(-360,300), expand=expansion(mult=c(0.02,0.02))) +
  labs(x="predicted accessibility attribution Δ (alt − ref)", y=NULL,
       title="Attribution ranks the allelic series") + base_th +
  theme(axis.text.y=element_text(size=6.6, lineheight=0.85))
sv("D.svg", pD, 180, 42)

cat("panels written:\n"); print(list.files("panels"))
