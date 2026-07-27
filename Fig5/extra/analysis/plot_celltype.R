source("_setup.R")
suppressPackageStartupMessages({ library(patchwork) })
VPOS<-72887998; TSS<-72887382
lo<-72887150; hi<-72888150
d <- read_parquet(file.path(AG,"spr_tracks_celltype.parquet")) %>%
  filter(output_type=="CAGE", strand=="+", pos>=lo, pos<=hi)
reg <- tribble(~biosample,~kind,~nice,
  "astrocyte of the cerebellum","cell-type resolved","cerebellar astrocyte",
  "astrocyte of the spinal cord","cell-type resolved","spinal-cord astrocyte",
  "astrocyte of the hippocampus","cell-type resolved","hippocampal astrocyte",
  "glutamatergic neuron","cell-type resolved","glutamatergic neuron",
  "cerebellum","bulk tissue","cerebellum (bulk)",
  "substantia nigra","bulk tissue","substantia nigra (bulk)")
dd <- d %>% inner_join(reg,by="biosample") %>% mutate(kb=pos/1e3)
rat <- dd %>% group_by(nice,kind,allele) %>% summarise(mx=max(value),.groups="drop") %>%
  pivot_wider(names_from=allele,values_from=mx) %>% mutate(pct=100*alt/ref,
    lab=sprintf("%s — alt %.0f%%", nice, pct))
ord <- c("cerebellar astrocyte","spinal-cord astrocyte","hippocampal astrocyte","glutamatergic neuron",
         "cerebellum (bulk)","substantia nigra (bulk)")
rat <- rat %>% mutate(nice=factor(nice,levels=ord)) %>% arrange(nice)
dd <- dd %>% mutate(nice=factor(nice,levels=ord)) %>% left_join(rat%>%select(nice,lab,kind),by=c("nice"))
dd <- dd %>% mutate(lab=factor(lab,levels=rat$lab))
stripcol <- ifelse(grepl("bulk", rat$lab), "#8C6D1F", "#1F5C8C")
p <- ggplot(dd, aes(kb,value)) +
  geom_vline(xintercept=TSS/1e3, linetype="dashed", linewidth=0.3, colour="grey45") +
  geom_area(data=~filter(.x,allele=="ref"), fill="grey60", colour="grey55", linewidth=0.2, alpha=0.65) +
  geom_area(data=~filter(.x,allele=="alt"), fill="#C0392B", colour="#C0392B", linewidth=0.25, alpha=0.45) +
  facet_wrap(~lab, ncol=2, scales="free_y") +
  scale_x_continuous(breaks=c(72.8873,72.8877,72.8881)) +
  labs(x="chr2 position (kb)     dashed = SPR TSS · grey = ref · red = alt",
       y="predicted CAGE (transcription initiation)") +
  theme_classic(base_size=8, base_family="inter") +
  theme(strip.background=element_blank(), strip.text=element_text(size=7,face="bold",hjust=0),
        axis.text=element_text(size=6.6,colour="grey20"), axis.title=element_text(size=7.4),
        panel.spacing=unit(3.2,"mm"), plot.margin=margin(4,7,4,6))
ggsave(file.path(FIGDIR,"FigureSPR_CAGE_celltype.pdf"), p, width=172, height=128, units="mm", device=cairo_pdf)
ggsave(file.path(FIGDIR,"FigureSPR_CAGE_celltype.png"), p, width=172, height=128, units="mm", dpi=300)
cat("wrote\n"); print(rat %>% select(kind,nice,ref,alt,pct) %>% mutate(pct=round(pct)))
