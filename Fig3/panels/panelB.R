# Panel B — signal-level colocalization: percent of genes by the number of
# colocalized signals, with an inset of per-assay credible-set coloc rate.
# Inputs: enloc_cs, susie_eqtl, susie_sqtl, susie_pqtl
source("fig3/00_setup.R")

cs_level        <- coloc_cs_components()
genes_all_three_assays <- susie_finemapping()$genes_all_three_assays

signal_counts <- cs_level %>%
  filter(!is.na(coloc_component)) %>%       # the original used na.omit()
  distinct(Gene, coloc_component) %>%
  dplyr::count(Gene) %>%
  dplyr::count(n, name = "nn")

main <- signal_counts %>%
  ggplot(aes(x = factor(n), y = nn / sum(nn), fill = factor(n))) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = SIGNAL_FILL, guide = "none") +
    # NO hard y limit. scale_y_continuous(limits = ) REMOVES out-of-range
    # values rather than clipping them, and the inherited cap of 0.55 silently
    # deleted the tallest bar once Fig3 moved to Fig2's shared SuSiE tables and
    # the "1 colocalized signal" share rose to 66% (2,298 of 3,483 genes). The
    # y axis text is hidden here, so a truncated bar cannot even be read as
    # truncated -- it just vanishes. See NOTES.md.
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Number of colocalized signals", y = "Percent of genes",
         title = "Signal-level Colocalization") +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = SZ_TITLE, face = "bold"),
          plot.title.position = "plot")

coloc_rate <- cs_level %>%
  filter(Gene %in% genes_all_three_assays) %>%
  group_by(type) %>%
  summarise(n_cs       = n_distinct(cs_id),
            n_coloc    = n_distinct(cs_id[has_coloc]),
            coloc_rate = n_coloc / n_cs,
            .groups = "drop")

inset <- coloc_rate %>%
  ggplot(aes(y = reorder(type, coloc_rate), x = coloc_rate, fill = type)) +
    geom_col(show.legend = FALSE) +
    paletteer::scale_fill_paletteer_d("ltc::minou") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = c(0, 0)) +
    labs(x = "% CS coloc", y = NULL) +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.y = element_line(color = "lightgrey"),
          axis.ticks.y = element_blank(),
          # No plot.tag = element_blank() override here: it would do nothing.
          # wrap_elements(full = ...) below folds `main` and this inset into
          # ONE opaque leaf before tagging happens, so recurse_tags() never
          # walks down into the inset to assign (or theme) a tag on it in the
          # first place -- confirmed empirically by rendering with `&
          # theme(plot.tag = element_text(color = "red"))` and observing no
          # tag, styled or not, on the inset either way.
          plot.margin = margin(0, 0, 0, 0))

pB_raw <- main + inset_element(inset, left = 0.48, bottom = 0.60,
                               right = 0.98, top = 0.98, align_to = "panel",
                               ignore_tag = TRUE)

# ignore_tag = TRUE keeps the inset itself out of the auto-tag sequence, but
# it has a sharp edge: patchwork gives "main + inset_element(...)" the class
# "inset_patch", and its auto-tagging walk (patchwork:::recurse_tags) skips
# ANY inset_patch object whose has_tag() is FALSE -- which, with ignore_tag =
# TRUE, is this whole composite, not just the inset. When pB_raw sits as a
# sibling among other panels (as it does once assemble.R does `P$A + P$B`),
# that means pB_raw's own main plot silently loses ITS tag too, and every
# later panel's letter shifts down by one -- confirmed empirically (see
# tests/test_fig3_assembly.R) before this wrap_elements() call was added.
# wrap_elements(full = ...) folds the whole composite into a single
# "wrapped_patch" leaf, which is not an inset_patch and so is exempt from
# that skip: it gets tagged normally as one panel, and the inset inside it
# still carries no tag of its own.
pB <- wrap_elements(full = pB_raw)

saveRDS(pB, file.path(OBJ_DIR, "B.rds"))
message("Panel B built.")
