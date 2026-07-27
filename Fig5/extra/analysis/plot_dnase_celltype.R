source("_setup.R")
VPOS<-72887998; TSS<-72887382
lo<-72887650; hi<-72888450   # zoom on the deleted element (summit ~72,888,043)
d <- read_parquet(file.path(AG,"spr_tracks_celltype.parquet")) %>%
  filter(output_type=="DNASE", pos>=lo, pos<=hi)
reg <- tribble(~biosample,~kind,~nice,
  "astrocyte of the cerebellum","cell-type resolved","cerebellar astrocyte",
  "astrocyte of the spinal cord","cell-type resolved","spinal-cord astrocyte",
  "astrocyte of the hippocampus","cell-type resolved","hippocampal astrocyte",
  "glutamatergic neuron","cell-type resolved","glutamatergic neuron",
  "brain","bulk tissue","brain (bulk)",
  "cerebellum","bulk tissue","cerebellum (bulk)")
dd <- d %>% inner_join(reg,by="biosample") %>% mutate(kb=pos/1e3)
pk <- dd %>% filter(allele=="ref") %>% group_by(nice) %>% summarise(rp=max(value),.groups="drop")
rat <- dd %>% group_by(nice,kind) %>% summarise(
  rp=max(value[allele=="ref"]), ap=max(value[allele=="alt"]),.groups="drop") %>%
  mutate(lab=sprintf("%s  (ref peak %.0f → alt %.0f%%)", nice, rp, 100*ap/rp))
ord <- c("cerebellar astrocyte","spinal-cord astrocyte","hippocampal astrocyte","glutamatergic neuron",
         "brain (bulk)","cerebellum (bulk)")
rat <- rat %>% mutate(nice=factor(nice,levels=ord)) %>% arrange(nice)
dd <- dd %>% mutate(nice=factor(nice,levels=ord)) %>% left_join(rat%>%select(nice,lab),by="nice") %>%
  mutate(lab=factor(lab,levels=rat$lab))
p <- ggplot(dd, aes(kb,value)) +
  annotate("rect", xmin=(VPOS-2)/1e3, xmax=(VPOS+3)/1e3, ymin=-Inf, ymax=Inf, fill="#F4A261", alpha=0.5) +
  geom_area(data=~filter(.x,allele=="ref"), fill="grey60", colour="grey55", linewidth=0.2, alpha=0.7) +
  geom_area(data=~filter(.x,allele=="alt"), fill="#C0392B", colour="#C0392B", linewidth=0.3, alpha=0.5) +
  facet_wrap(~lab, ncol=2, scales="fixed") +            # FIXED y -> baseline dilution visible
  scale_x_continuous(breaks=c(72.8879,72.8881,72.8883)) +
  labs(x="chr2 position (kb)   —   orange = deletion   ·   grey = ref DNase, red = alt   (shared y-axis)",
       y="predicted DNase accessibility") +
  theme_classic(base_size=8, base_family="inter") +
  theme(strip.background=element_blank(), strip.text=element_text(size=7,face="bold",hjust=0),
        axis.text=element_text(size=6.6,colour="grey20"), axis.title=element_text(size=7.2),
        panel.spacing=unit(3.2,"mm"), plot.margin=margin(4,7,4,6))
ggsave(file.path(FIGDIR,"FigureSPR_DNase_celltype.pdf"), p, width=172, height=128, units="mm", device=cairo_pdf)
ggsave(file.path(FIGDIR,"FigureSPR_DNase_celltype.png"), p, width=172, height=128, units="mm", dpi=300)
cat("wrote\n"); print(rat %>% select(kind,nice,rp,ap) %>% mutate(pct=round(100*ap/rp)))
