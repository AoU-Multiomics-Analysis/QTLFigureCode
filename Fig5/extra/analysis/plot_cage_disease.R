source("_setup.R")
VPOS <- 72887998; TSS <- 72887382
lo <- 72887150; hi <- 72888150
d <- arrow::read_parquet(file.path(AG,"spr_tracks_disease.parquet")) %>%
  filter(output_type=="CAGE", strand=="+", pos>=lo, pos<=hi)

reg <- tribble(~biosample,~grp,~nice,
  "substantia nigra","dopamine source","substantia nigra",
  "putamen","striatum","putamen",
  "nucleus accumbens","striatum","nucleus accumbens",
  "cerebellum","cerebellum","cerebellum",
  "astrocyte of the cerebellum","cerebellum","cerebellar astrocyte",
  "brain","reference","brain (bulk)")
dd <- d %>% inner_join(reg, by="biosample") %>% mutate(kb=pos/1e3)
rat <- dd %>% group_by(nice,grp,allele) %>% summarise(mx=max(value),.groups="drop") %>%
  pivot_wider(names_from=allele, values_from=mx) %>%
  mutate(lab=sprintf("%s (%s) — alt %.0f%% of ref", nice, grp, 100*alt/ref))
ord <- c("substantia nigra","putamen","nucleus accumbens","cerebellum","cerebellar astrocyte","brain (bulk)")
rat <- rat %>% mutate(nice=factor(nice,levels=ord)) %>% arrange(nice)
dd <- dd %>% mutate(nice=factor(nice,levels=ord)) %>%
  left_join(rat %>% select(nice,lab),by="nice") %>%
  mutate(lab=factor(lab,levels=rat$lab))

p <- ggplot(dd, aes(kb, value)) +
  geom_vline(xintercept=TSS/1e3, linetype="dashed", linewidth=0.3, colour="grey45") +
  geom_area(data=~filter(.x,allele=="ref"), fill="grey60", colour="grey55", linewidth=0.2, alpha=0.65) +
  geom_area(data=~filter(.x,allele=="alt"), fill="#C0392B", colour="#C0392B", linewidth=0.25, alpha=0.45) +
  facet_wrap(~lab, ncol=2, scales="free_y") +
  scale_x_continuous(breaks=c(72.8873,72.8877,72.8881)) +
  labs(x="chr2 position (kb)      dashed = SPR TSS · grey = ref · red = alt",
       y="predicted CAGE (transcription initiation)") +
  theme_classic(base_size=8, base_family="inter") +
  theme(strip.background=element_blank(), strip.text=element_text(size=7,face="bold",hjust=0),
        axis.text=element_text(size=6.6,colour="grey20"), axis.title=element_text(size=7.4),
        panel.spacing=unit(3.2,"mm"), plot.margin=margin(4,7,4,6))
ggsave(file.path(FIGDIR,"FigureSPR_CAGE_disease.pdf"), p, width=170, height=125, units="mm", device=cairo_pdf)
ggsave(file.path(FIGDIR,"FigureSPR_CAGE_disease.png"), p, width=170, height=125, units="mm", dpi=300)
cat("wrote FigureSPR_CAGE_disease\n"); print(rat %>% select(grp,nice,ref,alt) %>% mutate(pct=round(100*alt/ref)))
