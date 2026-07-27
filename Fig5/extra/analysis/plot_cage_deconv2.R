source("_setup.R")
d1 <- read_parquet(file.path(AG,"spr_cage_celltype.parquet")) %>% filter(strand=="+")
d2 <- read_parquet(file.path(AG,"spr_da_lineage.parquet")) %>% filter(output_type=="CAGE", strand=="+") %>% select(-output_type)
d <- bind_rows(d1,d2) %>% distinct(biosample,allele,pos,.keep_all=TRUE)
cls <- tribble(~biosample,~class,~nice,
 "oligodendrocyte precursor cell","glia / CNS non-neuronal","oligodendrocyte precursor",
 "astrocyte of the cerebellum","glia / CNS non-neuronal","cerebellar astrocyte",
 "astrocyte of the cerebral cortex","glia / CNS non-neuronal","cortical astrocyte",
 "central nervous system pericyte","glia / CNS non-neuronal","CNS pericyte",
 "leptomeningeal cell","glia / CNS non-neuronal","leptomeningeal cell",
 "Schwann cell","glia / CNS non-neuronal","Schwann cell",
 "SK-N-MC","catecholaminergic lineage","SK-N-MC (neuroblastoma)",
 "neuron","neuron","neuron (generic)",
 "neuronal stem cell","neuron","neuronal stem cell",
 "substantia nigra","bulk tissue (diluted)","substantia nigra (bulk)",
 "brain","bulk tissue (diluted)","brain (bulk)",
 "cerebellum","bulk tissue (diluted)","cerebellum (bulk)")
s <- d %>% inner_join(cls,by="biosample") %>% group_by(nice,class) %>%
  summarise(ref=max(value[allele=="ref"]), alt=max(value[allele=="alt"]),.groups="drop") %>%
  mutate(pct=round(100*alt/ref), nice=fct_reorder(nice,ref))
lv <- c("glia / CNS non-neuronal","catecholaminergic lineage","neuron","bulk tissue (diluted)")
s <- s %>% mutate(class=factor(class,levels=lv))
pal <- c("glia / CNS non-neuronal"="#1F6FB2","catecholaminergic lineage"="#0B7A3B",
         "neuron"="#7A4FA3","bulk tissue (diluted)"="#8C6D1F")
p <- ggplot(s, aes(y=nice)) +
  geom_segment(aes(x=alt,xend=ref,yend=nice,colour=class), linewidth=1.1, alpha=0.5) +
  geom_point(aes(x=ref,colour=class), size=2.6) +
  geom_point(aes(x=alt), colour="#C0392B", size=2.6) +
  geom_text(aes(x=ref, label=sprintf("  %d%%",pct)), hjust=0, size=2.3, family="inter", colour="grey30") +
  scale_colour_manual(values=pal, name=NULL) +
  scale_x_continuous(expand=expansion(mult=c(0.02,0.14))) +
  labs(x="predicted SPR TSS initiation (CAGE summit) — coloured = ref, red = alt (deletion)", y=NULL) +
  guides(colour=guide_legend(nrow=2)) +
  theme_classic(base_size=8, base_family="inter") +
  theme(legend.position="top", legend.text=element_text(size=6.9),
        axis.text.y=element_text(size=7.2,colour="grey15"), axis.text.x=element_text(size=6.8),
        axis.title.x=element_text(size=7.2), plot.margin=margin(4,8,4,6),
        legend.key.size=unit(3,"mm"), legend.margin=margin(0,0,0,0))
ggsave(file.path(FIGDIR,"FigureSPR_CAGE_deconv2.pdf"), p, width=162, height=104, units="mm", device=cairo_pdf)
ggsave(file.path(FIGDIR,"FigureSPR_CAGE_deconv2.png"), p, width=162, height=104, units="mm", dpi=300)
cat("wrote\n"); print(s %>% arrange(desc(ref)) %>% select(class,nice,ref,alt,pct))
