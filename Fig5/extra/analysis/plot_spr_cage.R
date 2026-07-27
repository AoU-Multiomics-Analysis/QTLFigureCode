source("_setup.R")
VPOS <- 72887998; TSS <- 72887382; SUMMIT <- 72887407
lo <- 72887150; hi <- 72888150     # 1 kb: TSS initiation peak + downstream deletion

d <- read_tsv(file.path(AG,"spr_tracks_window.tsv.gz"), show_col_types=FALSE) %>%
  filter(output_type=="CAGE", strand=="+", pos>=lo, pos<=hi)

tiss <- tribble(~biosample,~grp,~nice,
  "brain","Brain","brain",
  "astrocyte of the cerebellum","Brain","astrocyte (cerebellum)",
  "spinal cord","Brain","spinal cord",
  "putamen","Brain","putamen",
  "GM12878","Blood","GM12878 (lymphoblast)",
  "K562","Blood","K562 (myeloid)")

dd <- d %>% inner_join(tiss, by="biosample") %>% mutate(kb=pos/1e3)
rat <- dd %>% group_by(nice,grp,allele) %>% summarise(mx=max(value),.groups="drop") %>%
  pivot_wider(names_from=allele, values_from=mx) %>%
  mutate(lab=sprintf("%s   (%.0f%% of ref)", nice, 100*alt/ref))
ord <- rat %>% arrange(grp, desc(ref)) %>% pull(nice)
dd <- dd %>% mutate(nice=factor(nice, levels=ord)) %>%
  left_join(rat %>% select(nice,lab), by="nice") %>%
  mutate(lab=factor(lab, levels=rat %>% arrange(match(nice,ord)) %>% pull(lab)))

p <- ggplot(dd, aes(kb, value)) +
  annotate("rect", xmin=(VPOS-2)/1e3, xmax=(VPOS+3)/1e3, ymin=-Inf, ymax=Inf, fill="#F4A261", alpha=0.5) +
  annotate("text", x=(VPOS)/1e3, y=Inf, label="deletion", vjust=1.4, size=2.1, colour="#B8791F", family="inter") +
  geom_vline(xintercept=TSS/1e3, linetype="dashed", linewidth=0.3, colour="grey45") +
  geom_area(data=~filter(.x,allele=="ref"), fill="grey60", colour="grey55", linewidth=0.2, alpha=0.65) +
  geom_area(data=~filter(.x,allele=="alt"), fill="#C0392B", colour="#C0392B", linewidth=0.25, alpha=0.45) +
  facet_wrap(~lab, ncol=2, scales="free_y") +
  scale_x_continuous(breaks=c(72.8873,72.8877,72.8881), labels=c("72.8873","72.8877","72.8881")) +
  labs(x="chr2 position (kb)      dashed = SPR TSS · red = alt · grey = ref",
       y="predicted CAGE (transcription initiation)") +
  theme_classic(base_size=8, base_family="inter") +
  theme(strip.background=element_blank(), strip.text=element_text(size=7.2,face="bold",hjust=0),
        axis.text=element_text(size=6.6,colour="grey20"), axis.title=element_text(size=7.4),
        panel.spacing=unit(3.2,"mm"), plot.margin=margin(4,7,4,6))
ggsave(file.path(FIGDIR,"FigureSPR_CAGE.pdf"), p, width=170, height=125, units="mm", device=cairo_pdf)
ggsave(file.path(FIGDIR,"FigureSPR_CAGE.png"), p, width=170, height=125, units="mm", dpi=300)
cat("wrote FigureSPR_CAGE\n"); print(rat %>% select(grp,nice,ref,alt) %>% mutate(pct=round(100*alt/ref)))
