source("_setup.R")
d <- read_tsv(file.path(AG,"spr_rna_celltype.tsv"), show_col_types=FALSE) %>%
  mutate(class=recode(class, glia="glia / CNS non-neuronal", catecholaminergic="catecholaminergic lineage",
                      neuron="neuron", bulk="bulk tissue (diluted)"),
         class=factor(class, levels=c("glia / CNS non-neuronal","catecholaminergic lineage","neuron","bulk tissue (diluted)")),
         nice=fct_reorder(biosample, base))
pal <- c("glia / CNS non-neuronal"="#1F6FB2","catecholaminergic lineage"="#0B7A3B",
         "neuron"="#7A4FA3","bulk tissue (diluted)"="#8C6D1F")
p <- ggplot(d, aes(base, nice, colour=class)) +
  geom_segment(aes(x=0,xend=base,yend=nice), linewidth=0.9, alpha=0.5) +
  geom_point(size=2.8) +
  geom_text(aes(label=sprintf("  LFC %+.2f", lfc)), hjust=0, size=2.3, family="inter", colour="grey30") +
  scale_colour_manual(values=pal, name=NULL) +
  scale_x_continuous(expand=expansion(mult=c(0.01,0.28))) +
  labs(x="predicted baseline SPR expression (RNA-seq)   ·   text = variant effect (log2 LFC)", y=NULL) +
  guides(colour=guide_legend(nrow=2)) +
  theme_classic(base_size=8, base_family="inter") +
  theme(legend.position="top", legend.text=element_text(size=6.9),
        axis.text.y=element_text(size=7.2,colour="grey15"), axis.text.x=element_text(size=6.8),
        axis.title.x=element_text(size=7.2), legend.key.size=unit(3,"mm"),
        legend.margin=margin(0,0,0,0), plot.margin=margin(4,8,4,6))
ggsave(file.path(FIGDIR,"FigureSPR_RNA_celltype.pdf"), p, width=160, height=100, units="mm", device=cairo_pdf)
ggsave(file.path(FIGDIR,"FigureSPR_RNA_celltype.png"), p, width=160, height=100, units="mm", dpi=300)
cat("wrote\n")
