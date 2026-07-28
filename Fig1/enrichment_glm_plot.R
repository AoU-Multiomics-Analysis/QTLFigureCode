# =============================================================================
# Plots the PIP logistic-regression enrichment fitted by fig1/enrichment_glm.R.
# Sourced from that script, but also runnable on its own once the results table
# exists:
#   Rscript fig1/enrichment_glm_plot.R
# =============================================================================
if (!exists("OBJ_DIR")) source("fig1/00_setup.R")
if (!exists("results")) {
  results <- data.table::fread(file.path(OBJ_DIR, "enrichment_glm.tsv"))
  message("Read ", file.path(OBJ_DIR, "enrichment_glm.tsv"))
}

# Every annotation the three tables produce, grouped the way panel H groups its
# twelve. Unlike panel H this is a total mapping: an annotation missing from it
# is an error below rather than a row that silently vanishes.
GLM_GROUPS <- tibble::tribble(
  ~name,                                ~functional_group,
  "coding_sequence_variant",            "Coding",
  "frameshift_variant",                 "Coding",
  "inframe_deletion",                   "Coding",
  "inframe_insertion",                  "Coding",
  "missense_variant",                   "Coding",
  "synonymous_variant",                 "Coding",
  "start_lost",                         "Coding",
  "start_retained_variant",             "Coding",
  "stop_gained",                        "Coding",
  "stop_lost",                          "Coding",
  "stop_retained_variant",              "Coding",
  "incomplete_terminal_codon_variant",  "Coding",
  "splice_acceptor_variant",            "Splice",
  "splice_donor_variant",               "Splice",
  "splice_region_variant",              "Splice",
  "3_prime_UTR_variant",                "Regulatory",
  "5_prime_UTR_variant",                "Regulatory",
  "PLS",                                "Regulatory",
  "pELS",                               "Regulatory",
  "dELS",                               "Regulatory",
  "CTCF.bound",                         "Regulatory",
  "CTCF.only",                          "Regulatory",
  "DNase.H3K4me3",                      "Regulatory",
  "FANTOM5",                            "Regulatory",
  # NMD_transcript_variant is a transcript-biotype flag rather than a
  # consequence on the variant itself; it is grouped here for want of a better
  # home, not because it is a noncoding-RNA annotation.
  "NMD_transcript_variant",             "Noncoding RNA",
  "non_coding_transcript_exon_variant", "Noncoding RNA",
  "non_coding_transcript_variant",      "Noncoding RNA",
  "mature_miRNA_variant",               "Noncoding RNA",
  "intron_variant",                     "Intronic / intergenic",
  "upstream_gene_variant",              "Intronic / intergenic",
  "downstream_gene_variant",            "Intronic / intergenic"
)

GLM_LABELS <- c(
  coding_sequence_variant            = "Coding sequence",
  frameshift_variant                 = "Frameshift",
  inframe_deletion                   = "Inframe deletion",
  inframe_insertion                  = "Inframe insertion",
  missense_variant                   = "Missense",
  synonymous_variant                 = "Synonymous",
  start_lost                         = "Start lost",
  start_retained_variant             = "Start retained",
  stop_gained                        = "Stop gained",
  stop_lost                          = "Stop lost",
  stop_retained_variant              = "Stop retained",
  incomplete_terminal_codon_variant  = "Incomplete terminal codon",
  splice_acceptor_variant            = "Splice acceptor",
  splice_donor_variant               = "Splice donor",
  splice_region_variant              = "Splice region",
  `3_prime_UTR_variant`              = "3' UTR",
  `5_prime_UTR_variant`              = "5' UTR",
  PLS                                = "Promoter (PLS)",
  pELS                               = "Proximal enhancer (pELS)",
  dELS                               = "Distal enhancer (dELS)",
  CTCF.bound                         = "CTCF bound",
  CTCF.only                          = "CTCF only",
  DNase.H3K4me3                      = "DNase / H3K4me3",
  FANTOM5                            = "FANTOM5 enhancer",
  NMD_transcript_variant             = "NMD transcript",
  non_coding_transcript_exon_variant = "Noncoding transcript exon",
  non_coding_transcript_variant      = "Noncoding transcript",
  mature_miRNA_variant               = "Mature miRNA",
  intron_variant                     = "Intron",
  upstream_gene_variant              = "Upstream of gene",
  downstream_gene_variant            = "Downstream of gene"
)

unmapped <- setdiff(unique(results$name), GLM_GROUPS$name)
if (length(unmapped))
  stop("annotations with no functional group: ", paste(unmapped, collapse = ", "))
unlabelled <- setdiff(unique(results$name), names(GLM_LABELS))
if (length(unlabelled))
  stop("annotations with no label: ", paste(unlabelled, collapse = ", "))

group_levels <- c("Coding", "Splice", "Regulatory", "Noncoding RNA",
                  "Intronic / intergenic")

# Slopes are shown exponentiated on a log10 axis, matching panel H: the effect
# reads as an odds ratio between PIP = 1 and PIP = 0, and a depletion sits below
# 1 rather than below 0. exp() is monotone, so the bounds transform directly.
plot_data <- results %>%
  filter(!separated) %>%
  left_join(GLM_GROUPS, by = "name") %>%
  mutate(
    label            = GLM_LABELS[name],
    functional_group = factor(functional_group, levels = group_levels),
    type             = factor(type, levels = c("eQTL", "pQTL", "sQTL")),
    significant      = padj < 0.05,
    or               = exp(estimate),
    or_low           = exp(conf.low),
    or_high          = exp(conf.high)
  )

# Annotations are ordered by their mean slope across assays, and each facet is
# given exactly as much height as it has annotations.
label_order <- plot_data %>%
  group_by(label) %>%
  summarise(m = mean(estimate), .groups = "drop") %>%
  arrange(m) %>%
  pull(label)
plot_data <- plot_data %>% mutate(label = factor(label, levels = label_order))

dodge <- position_dodge(width = 0.7)

pGLM <- ggplot(plot_data, aes(x = or, y = label, fill = type, alpha = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30") +
  geom_linerange(aes(xmin = or_low, xmax = or_high, color = type),
                 position = dodge, show.legend = FALSE, linewidth = 0.8) +
  geom_point(shape = 21, size = 5, color = "black", stroke = 0.4, position = dodge) +
  facet_grid(functional_group ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(labels = scales::label_log()) +
  scale_fill_manual(values = AssayColors, name = "Assay") +
  scale_color_manual(values = AssayColors) +
  scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 1), name = NULL,
                     labels = c(`FALSE` = "FDR ≥ 0.05", `TRUE` = "FDR < 0.05")) +
  labs(x = "Odds ratio of the annotation per unit PIP",
       y = NULL,
       title = "Annotation enrichment as a logistic regression on PIP",
       subtitle = paste("glm(annotated ~ PIP, binomial) over every fine-mapped variant.",
                        "Bars are 95% confidence intervals.")) +
  guides(fill  = guide_legend(override.aes = list(alpha = 1, size = 6)),
         alpha = guide_legend(override.aes = list(fill = "grey30", size = 6))) +
  theme(strip.text.y      = element_text(angle = 0, hjust = 0),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        legend.position   = "bottom",
        plot.title        = element_text(size = SZ_AXIS_TITLE, face = "bold"),
        plot.subtitle     = element_text(size = SZ_LEG_TITLE, color = "grey30"))

ggsave(file.path(OUTDIR, "Figure1_enrichment_glm.pdf"), pGLM,
       width = 16, height = 20, dpi = DPI, limitsize = FALSE)
showtext_opts(dpi = 100)
ggsave(file.path(OUTDIR, "Figure1_enrichment_glm.png"), pGLM,
       width = 16, height = 20, dpi = 100, limitsize = FALSE)
showtext_opts(dpi = DPI)

message("Wrote ", file.path(OUTDIR, "Figure1_enrichment_glm.pdf"), " and .png")
