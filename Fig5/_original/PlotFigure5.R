## ============================================================================
## Figure 5 — AoU rare-variant QTL paper
##
## Panels
##   A  fine-mapped rare ClinGen variants by ancestry (AD/AR) + ClinVar prop.
##   B  |average QTL effect| by variant-effect category (Control/AD/AR)
##   C  log2 aFC of the strongest decreasing variant per gene (AD/AR)
##   D  ACVRL1 fine-mapping locus (PIP) with zoom inset
##   E  ACVRL1 causal promoter variants — PromoterAI base-resolution effect
##   F  SPR (second example) — Borzoi multi-track allelic browser + motif logo
##
## Cleaned to run from local files.  Private AoU gs:// inputs are optional
## (guarded); ClinGen is rebuilt from the public gene-validity download.
## ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(patchwork)
  library(rcartocolor)
  library(colorspace)
  library(showtext); library(sysfonts)
  library(scales)
})

font_add_google("Inter", "inter"); showtext_opts(dpi = 400); showtext_auto()
theme_set(theme_classic())
theme_update(strip.text = element_text(size = 11),
             text = element_text(family = "inter"),
             axis.title = element_text(size = 12),
             axis.text  = element_text(size = 10),
             strip.background = element_blank(),
             plot.tag = element_text(size = 20, face = "bold"),
             plot.margin = margin(9, 12, 9, 12))   # breathing room between panels

## ---------------------------------------------------------------------------
## paths (local)
## ---------------------------------------------------------------------------
DIR_FIG   <- "/Users/evinmpadhi/Fig5_QTL_paper"
PATH_GB   <- "/Users/evinmpadhi/Downloads/41588_2024_1820_MOESM4_ESM.xlsx"      # GeneBayes
PATH_SUSIE<- "/Users/evinmpadhi/Desktop/susie_files/COMB.expression_SusieMerged.annotated.ancestry_skew.tsv.gz"
PATH_PAI  <- "/Users/evinmpadhi/Downloads/promoterAI_tss500.tsv.gz"             # PromoterAI
PATH_CG   <- file.path(DIR_FIG, "clingen_gene_validity.csv")                    # ClinGen (public)
PATH_NOHOM<- file.path(DIR_FIG, "NoHomozygoteVariants.txt")                     # private (optional)

inheritance_colors <- c(Control = "#D43F3A", AD = "#EEA236", AR = "#5CB85C")

## signed-square transform for the PromoterAI y-axis
make_signed_square_trans <- function()
  scales::trans_new("signed_square",
    transform = function(x) sign(x) * x^2,
    inverse   = function(x) sign(x) * sqrt(abs(x)))

## ---------------------------------------------------------------------------
## load reference / annotation tables
## ---------------------------------------------------------------------------
GeneBayes <- read_excel(PATH_GB, sheet = "Supplementary Table 1") %>%
  transmute(ensg, post_mean)

## ClinGen gene-validity -> Definitive/Strong genes with mode of inheritance.
## (rebuilds the private Open Targets clingen parquet from the public download)
clingen_raw <- fread(PATH_CG, skip = 6, header = FALSE, select = 1:7,
  col.names = c("gene","hgnc","disease","mondo","moi","sop","confidence"))

## symbol -> ensembl map (version-stripped) from the susie file, to attach targetId
susie_gene_map <- fread(PATH_SUSIE, select = c("gene_id","gene_name")) %>%
  distinct() %>% mutate(gene_id = sub("\\..*", "", gene_id))

clingen_cleaned_list <- clingen_raw %>%
  filter(confidence %in% c("Definitive","Strong"), moi %in% c("AD","AR","SD")) %>%
  distinct(gene, moi, confidence) %>%
  transmute(targetFromSourceId = gene, allelicRequirements = moi, confidence) %>%
  left_join(susie_gene_map, by = c("targetFromSourceId" = "gene_name")) %>%
  rename(targetId = gene_id)

## private no-homozygote variant list (used to refine Panel C); optional
if (file.exists(PATH_NOHOM)) {
  NoHomozygotes <- fread(PATH_NOHOM) %>% mutate(variant = paste(V1, V2, V3, V4, sep = "_"))
  have_nohom <- TRUE
} else {
  message("NOTE: ", basename(PATH_NOHOM), " not found (private gs:// file) — ",
          "Panel C selects the strongest decreasing variant per gene without it.")
  have_nohom <- FALSE
}

## ---------------------------------------------------------------------------
## eQTL SuSiE fine-mapping (+ GeneBayes shet, ClinGen flag, rarity)
## ---------------------------------------------------------------------------
eQTL_susie <- fread(PATH_SUSIE) %>%
  mutate(gene_id_nov = sub("\\..*", "", gene_id)) %>%
  left_join(GeneBayes, by = c("gene_id_nov" = "ensg")) %>%
  mutate(ClinGen = gene_name %in% clingen_cleaned_list$targetFromSourceId,
         rare    = gvs_max_af < 0.01)

## ---- per-variant rarity / ancestry annotation (used across D-G) -------------
## rare = max sub-population AF < 1%. Labels report the population-max MAF and the
## sub-population it is highest in.
.fmt_af <- function(af) dplyr::case_when(
  is.na(af)   ~ "NA",
  af >= 0.01  ~ sprintf("%.1f%%",  af * 100),
  af >= 0.001 ~ sprintf("%.2f%%",  af * 100),
  TRUE        ~ sprintf("%.3f%%",  af * 100))
make_freq_label <- function(max_af, subpop, rare) {
  pop  <- ifelse(is.na(subpop), "", sprintf(" (%s)", toupper(subpop)))
  base <- sprintf("MAF %s%s", .fmt_af(max_af), pop)
  ifelse(rare, paste0("rare · ", base), base)
}
variant_freq <- eQTL_susie %>%
  transmute(variant_key = gsub("_+", "_", gsub("[:*]", "_", variant)),
            max_af = gvs_max_af, subpop = gvs_max_subpop, rare, log2_aFC) %>%
  distinct(variant_key, .keep_all = TRUE) %>%
  mutate(freq_label = make_freq_label(max_af, subpop, rare),
         afc_label  = sprintf("log₂ aFC %+.1f", log2_aFC))

## ---------------------------------------------------------------------------
## variant effects by inheritance (ClinGen genes vs low-shet control genes)
## ---------------------------------------------------------------------------
gene_assignments <- clingen_cleaned_list %>%
  distinct(targetFromSourceId, allelicRequirements) %>%
  filter(allelicRequirements != "SD")

variant_effect_by_inheritance <- eQTL_susie %>%
  mutate(ControlGene = post_mean < 0.01 & ClinGen == FALSE) %>%
  filter((ClinGen == TRUE | ControlGene == TRUE), pip > 0.9, rare == TRUE) %>%
  select(variant, log2_aFC, gene_name, posterior_mean) %>%
  left_join(gene_assignments, by = c("gene_name" = "targetFromSourceId")) %>%
  distinct() %>%
  mutate(allelicRequirements = factor(replace_na(allelicRequirements, "Control"),
                                      levels = c("Control","AD","AR")))

## ---------------------------------------------------------------------------
## Panel B — |average QTL effect| by direction category
## ---------------------------------------------------------------------------
sign_summary <- variant_effect_by_inheritance %>%
  mutate(dir = if_else(log2_aFC < 0, "Decreasing", "Increasing")) %>%
  group_by(allelicRequirements, dir) %>%
  summarise(MeanEffect = mean(abs(posterior_mean), na.rm = TRUE),
            se = sd(abs(posterior_mean), na.rm = TRUE) / sqrt(n()), .groups = "drop")

aggregate_summary <- variant_effect_by_inheritance %>%
  group_by(allelicRequirements) %>%
  summarise(MeanEffect = mean(abs(posterior_mean), na.rm = TRUE),
            se = sd(abs(posterior_mean), na.rm = TRUE) / sqrt(n()), .groups = "drop") %>%
  mutate(dir = "Aggregate")

average_effect_summary <- bind_rows(sign_summary, aggregate_summary) %>%
  mutate(dir = factor(dir, levels = c("Aggregate","Decreasing","Increasing")),
         conf.low = MeanEffect - 1.96 * se, conf.high = MeanEffect + 1.96 * se)

category_lines <- tibble(x = c(0.60, 1.60, 2.60), xend = c(1.40, 2.40, 3.40))

average_effect_by_gene_category <- ggplot(average_effect_summary,
    aes(dir, MeanEffect, color = allelicRequirements)) +
  geom_segment(data = category_lines, aes(x = x, xend = xend, y = -Inf, yend = -Inf),
               inherit.aes = FALSE, color = "grey72", linewidth = 0.55, lineend = "round") +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high),
                 position = position_dodge(width = 0.5), linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.5), size = 2.3) +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(values = inheritance_colors, drop = FALSE) +
  labs(x = "Variant Effect Category", y = "|Average QTL effect|", color = NULL) +
  coord_cartesian(clip = "off") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_text(margin = margin(t = 1)),
        legend.position = c(0.5, 0.96), legend.direction = "horizontal")

## ---------------------------------------------------------------------------
## Panel C — strongest decreasing variant per gene (AD/AR), lollipop
## ---------------------------------------------------------------------------
n_genes_to_show <- 7
required_genes  <- c("HEXA", "MLC1")

candidate_gene_variants <- variant_effect_by_inheritance %>%
  {if (have_nohom) filter(., variant %in% NoHomozygotes$variant) else .} %>%
  filter(allelicRequirements %in% c("AD","AR"), !is.na(log2_aFC), log2_aFC < 0) %>%
  group_by(gene_name, allelicRequirements) %>%
  slice_min(log2_aFC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(required_gene = gene_name %in% required_genes)

high_effect_genes <- candidate_gene_variants %>%
  group_by(allelicRequirements) %>%
  arrange(desc(required_gene), log2_aFC, .by_group = TRUE) %>%
  slice_head(n = n_genes_to_show) %>% ungroup() %>%
  mutate(allelicRequirements = factor(allelicRequirements, levels = c("AD","AR")),
         display_gene = fct_reorder(gene_name, log2_aFC, .desc = TRUE))

high_effect_plot <- ggplot(high_effect_genes,
    aes(log2_aFC, display_gene, color = allelicRequirements)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = log2_aFC, y = display_gene, yend = display_gene),
               color = "grey78", linewidth = 0.65) +
  geom_point(size = 3.2) +
  facet_grid(allelicRequirements ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(position = "right") +
  scale_color_manual(values = inheritance_colors, breaks = c("AD","AR"), guide = "none") +
  scale_x_continuous(limits = c(-7, 0), breaks = c(-6,-4,-2,0),
                     expand = expansion(mult = c(0.03, 0.02))) +
  labs(x = expression("log"[2]*" aFC (decreasing variants only)"), y = NULL) +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(face = "bold", angle = 0),
        axis.text.y.right = element_text(hjust = 0),
        panel.spacing.y = unit(1, "lines"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(5.5, 20, 5.5, 5.5))

combined_effect_plot <- ((average_effect_by_gene_category + labs(tag = "B")) |
                         (high_effect_plot + labs(tag = "C"))) +
  plot_layout(widths = c(0.85, 1.05))

## ---------------------------------------------------------------------------
## Panel A — fine-mapped rare ClinGen variants by ancestry + ClinVar proportions
## ---------------------------------------------------------------------------
ClinGenOnly <- eQTL_susie %>%
  filter(ClinGen == TRUE) %>%
  left_join(clingen_cleaned_list, by = c("gene_id_nov" = "targetId")) %>%
  mutate(causal = !is.na(log2_aFC))

RareClinGenVariants <- ClinGenOnly %>%
  filter(rare == TRUE) %>%
  mutate(InClinvar = !is.na(clinvar_classification))

base_cols <- setNames(carto_pal(7, "Temps"),
                      c("afr","amr","eas","eur","mid","admix","sas"))
fill_cols <- unlist(lapply(names(base_cols), function(pop)
  c(setNames(lighten(unname(base_cols[pop]), 0.15), paste0(pop, ".FALSE")),
    setNames(darken(unname(base_cols[pop]), 0.10),  paste0(pop, ".TRUE")))))

ancestry_count <- RareClinGenVariants %>%
  count(allelicRequirements, gvs_max_subpop, InClinvar) %>%
  mutate(gvs_max_subpop = str_remove_all(gvs_max_subpop, "gvs_|_af"),
         gvs_max_subpop = if_else(gvs_max_subpop == "oth", "admix", gvs_max_subpop)) %>%
  filter(allelicRequirements != "SD", !is.na(gvs_max_subpop))

clinvar_plot <- ancestry_count %>%
  ggplot(aes(gvs_max_subpop, n, fill = interaction(gvs_max_subpop, InClinvar))) +
  geom_col(color = "transparent") +
  geom_point(aes(color = InClinvar), x = Inf, y = Inf, inherit.aes = FALSE) +
  facet_wrap(~allelicRequirements) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_x_discrete(labels = toupper) +                  # ADMIX / AFR / AMR / ... uppercase
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "grey30"),
                     labels = c(`FALSE` = "Novel Variant", `TRUE` = "ClinVar Annotated"),
                     name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
  geom_segment(aes(x = 0.6, xend = 7.4, y = Inf, yend = Inf), inherit.aes = FALSE,
               color = "grey70", linewidth = 0.6) +
  labs(y = "Fine-mapped Rare Variants") +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), legend.position = c(0.9, 0.85))

ClinVarProportions <- RareClinGenVariants %>%
  mutate(clinvar_class = case_when(
    str_detect(clinvar_classification, "pathogenic") ~ "P/LP",
    str_detect(clinvar_classification, "uncertain significance") ~ "VUS",
    str_detect(clinvar_classification, "benign") ~ "B/LB",
    TRUE ~ "Unannotated")) %>%
  count(allelicRequirements, clinvar_class) %>%
  filter(allelicRequirements != "SD") %>%
  group_by(allelicRequirements) %>% mutate(prop = n / sum(n)) %>% ungroup()

ClinVarPropPlot <- ClinVarProportions %>%
  ggplot(aes(prop, allelicRequirements, fill = clinvar_class)) +
  geom_col(position = "stack", width = 0.85) +
  facet_wrap(~allelicRequirements, nrow = 1, scales = "free_y") +
  paletteer::scale_fill_paletteer_d("futurevisions::mars") +
  ## 0/25/50/75/100% ticks sit directly under the bars, no axis line; extra facet
  ## spacing keeps the 100% (AD) and 0% (AR) labels from colliding
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1),
                     labels = \(x) paste0(x * 100, "%"), expand = expansion(mult = 0)) +
  scale_y_discrete(expand = expansion(add = 0.06)) +   # bar fills the panel -> labels sit on it
  guides(fill = guide_legend(nrow = 1)) +
  labs(fill = NULL, x = NULL) +
  theme(axis.text.y = element_blank(), axis.line.y = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_blank(),                       # <- no x axis line
        axis.text.x = element_text(size = 9.5, margin = margin(t = 0)),
        axis.ticks.x = element_blank(),                      # <- no x ticks (they got in the way)
        strip.text = element_blank(), panel.spacing.x = unit(1.6, "lines"),
        plot.margin = margin(0, 8, 0, 8),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.box.spacing = unit(2, "pt"), legend.margin = margin(t = 0, b = 0),  # tighten legend gap
        legend.text = element_text(size = 9), legend.title = element_text(size = 9),
        legend.key.width = unit(0.45, "cm"), legend.key.height = unit(0.35, "cm"))

panelA <- (clinvar_plot + labs(tag = "A")) / ClinVarPropPlot + plot_layout(heights = c(7, 1))

top_row <- panelA | combined_effect_plot + plot_layout(widths = c(1.5, 1))

## ---------------------------------------------------------------------------
## Panels D/E — ACVRL1 example
## ---------------------------------------------------------------------------
source(file.path(DIR_FIG, "acvrl1_panels.R"))   # defines acvrl1_finemap + acvrl1_promoterAI panels

## D tag/title set inside the fine-mapping panel; shrink D further (white space)
## D and E kept as separate cells so D can share a row with F (guarantees identical size)
acvrl1_mechanism <- acvrl1_promoterAI_plot + labs(tag = "E")

## ---------------------------------------------------------------------------
## Panel F — SPR Borzoi browser (second worked example)
## ---------------------------------------------------------------------------
source(file.path(DIR_FIG, "spr_browser_panel.R"))

SPR_QTL <- eQTL_susie %>% filter(gene_name == "SPR") %>%
  mutate(causal = pip > 0.9, cs_id = factor(cs_id), type = "Expression",
         is_focal = position == SPR_VARIANT_POS)
spr_cs_lv  <- levels(droplevels(SPR_QTL$cs_id))
spr_cs_pal <- setNames(carto_pal(max(3, length(spr_cs_lv)), "Prism")[seq_along(spr_cs_lv)], spr_cs_lv)

## Panel F — SPR fine-map in the same D style: focused regional view (distal causal
## + 72.888 cluster + SPR gene) with a zoom inset on the two clustered causal variants.
spr_ccre <- fread(file.path(DIR_FIG, "spr_ccre_region.tsv"))   # cCREs across the F locus
SPR_region <- eQTL_susie %>%
  filter(gene_name == "SPR", between(position / 1e6, 72.881, 72.8925)) %>%
  mutate(causal = pip > 0.9, cs_id = as.character(cs_id), type = "Expression")

spr_finemap <- plot_finemapping_with_inset(
  SPR_region, pip_thresh = 0.9,
  zoom_center = 72887978, zoom_window_bp = 200,   # tight zoom separates the two causal
  x_limits = c(72.8817, 72.8934),                 # clear margin past the full SPR gene model
  cs_colors = spr_cs_pal, panel_tag = "F", panel_title = "SPR eQTL fine-mapping",
  ccre = spr_ccre)

## SPR browser tracks get their own panel label (G) in a thin header row above the
## tracks (avoids the tag overlapping the first track label); annotate the browser
## variant's identity + rarity/ancestry so G reads on its own
spr_focal <- SPR_QTL %>% filter(is_focal) %>%
  mutate(variant_key = gsub("_+", "_", gsub("[:*]", "_", variant))) %>%
  left_join(variant_freq %>% select(variant_key, freq_label, afc_label), by = "variant_key") %>%
  dplyr::slice(1)
spr_focal_color <- unname(spr_cs_pal[as.character(spr_focal$cs_id)])
spr_focal_vid   <- sprintf("chr2:%s %s>%s",
                           format(spr_focal$position, big.mark = ",", trim = TRUE),
                           spr_focal$ref, spr_focal$alt)
## G header, styled like an E facet strip: variant id / MAF+ancestry / log2 aFC stacked
## and centred, with the browser variant's (ring-free) diamond below. Text matches E
## strips (9 pt bold). Bold "G" panel tag sits at the far left.
## single multi-line label (lineheight/size = E's 9 pt bold strips) so the text
## block height matches E, with the ring-free diamond below
G_header <- ggplot() +
  annotate("text", x = 0.012, y = 0.5, label = "G", hjust = 0, vjust = 0.5,
           fontface = "bold", size = 7, family = "inter") +
  annotate("text", x = 0.5, y = 0.94, hjust = 0.5, vjust = 1, lineheight = 0.9,
           label = paste(spr_focal_vid, spr_focal$freq_label, spr_focal$afc_label, sep = "\n"),
           size = 3.16, fontface = "bold", family = "inter", colour = "grey15") +
  annotate("point", x = 0.5, y = 0.16, shape = 23, size = 3.6,
           fill = spr_focal_color, colour = "black", stroke = 0.8) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void()
spr_browser_g <- cowplot::plot_grid(G_header, make_spr_browser(),
                                    ncol = 1, rel_heights = c(0.14, 1))

spr_browser_wrapped <- wrap_elements(full = spr_browser_g)

## ---------------------------------------------------------------------------
## assemble Figure 5
##   bottom half is a 2x2 grid: fine-map row (D | F) over mechanism row (E | G).
##   Putting D and F in the SAME row guarantees they are exactly the same size.
## ---------------------------------------------------------------------------
finemap_row   <- acvrl1_finemap_plot | spr_finemap
mechanism_row <- acvrl1_mechanism | spr_browser_wrapped
## thin spacer rows give the horizontal bands more breathing room
bottom_half   <- finemap_row / plot_spacer() / mechanism_row +
  plot_layout(heights = c(1, 0.10, 2.33))

figure5 <- top_row / plot_spacer() / bottom_half +
  plot_layout(heights = c(1, 0.08, 2.43))

ggsave(file.path(DIR_FIG, "Figure5.pdf"), figure5,
       width = 18, height = 16, bg = "white", limitsize = FALSE)
message("wrote Figure5.pdf")
