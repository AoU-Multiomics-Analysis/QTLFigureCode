# =============================================================================
# Supplemental — proportion of genes with a cross-QTL colocalization, by assay
# pair. NOT part of Figure 3. Produces FigS_CrossQTLPropGeneColoc.pdf.
#
# Run from the project root:  Rscript fig3/supp_prop_gene_coloc.R
# =============================================================================
source("fig3/00_setup.R")

ProportionGeneColocs <- coloc_gene_level() %>%
  group_by(type) %>%
  summarize(PropColoc = mean(Coloc), .groups = "drop") %>%
  ggplot(aes(x = type, y = PropColoc)) +
    geom_col(position = "dodge") +
    labs(x = NULL, y = "Proportion Colocalized") +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())

# Explicit plot = argument: the original relied on ggsave()'s last-plot default,
# which is correct in a notebook but fragile under Rscript.
ggsave(file.path(OUTDIR, "FigS_CrossQTLPropGeneColoc.pdf"), ProportionGeneColocs,
       height = 4, width = 6, dpi = DPI)
message("Wrote ", file.path(OUTDIR, "FigS_CrossQTLPropGeneColoc.pdf"))
