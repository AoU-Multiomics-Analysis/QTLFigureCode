# Fig5 pure-seam unit tests. Run from the project root (via tests/run_all.R).
# Sources fig5/00_setup.R for the folded seams; uses in-memory fixture frames so
# no controlled-access data, Excel, or gzip reads are needed.
source("fig5/00_setup.R")
source("tests/fixtures.R")
source("tests/helpers.R")

fx <- make_fig5_fixture()

# --- derive_clingen ---------------------------------------------------------
cg <- derive_clingen(fx$clingen_raw, fx$susie_gene_map)
ok("derive_clingen keeps only Definitive/Strong AD/AR/SD genes",
   setequal(unique(cg$targetFromSourceId), c("HEXA", "MLC1", "GENE_AD", "GENE_BOTH")))
ok("derive_clingen drops Limited-confidence and non-AD/AR/SD rows",
   !any(c("LIMGENE", "XLGENE") %in% cg$targetFromSourceId))
ok("derive_clingen keeps a gene curated under BOTH AD and AR as two rows",
   sum(cg$targetFromSourceId == "GENE_BOTH") == 2L &&
     setequal(cg$allelicRequirements[cg$targetFromSourceId == "GENE_BOTH"], c("AD", "AR")))
ok("derive_clingen maps gene_name -> Ensembl targetId via the SuSiE map",
   cg$targetId[cg$targetFromSourceId == "HEXA"][1] == "ENSGH")

# --- annotate_eqtl ----------------------------------------------------------
clingen_genes <- cg$targetFromSourceId
eqtl <- annotate_eqtl(fx$susie, fx$genebayes, clingen_genes)
ok("annotate_eqtl flags ClinGen genes",
   all(eqtl$ClinGen[eqtl$gene_name %in% c("HEXA", "MLC1", "GENE_AD", "GENE_BOTH")]))
ok("annotate_eqtl leaves non-ClinGen genes unflagged",
   !any(eqtl$ClinGen[eqtl$gene_name %in% c("CTRL1", "NONCLIN")]))
ok("annotate_eqtl sets rare = gvs_max_af < 0.01",
   identical(eqtl$rare, eqtl$gvs_max_af < 0.01) && !eqtl$rare[eqtl$gvs_max_af == 0.02][1])
ok("annotate_eqtl joins GeneBayes post_mean (CTRL1 is low-shet)",
   abs(eqtl$post_mean[eqtl$gene_name == "CTRL1"][1] - 0.005) < 1e-9)

# --- build_variant_freq -----------------------------------------------------
vf <- build_variant_freq(eqtl)
ok("build_variant_freq returns one row per distinct variant", nrow(vf) == 7L)
ok("build_variant_freq normalizes the variant key (colons -> underscores)",
   all(!grepl(":", vf$variant_key)))
ok("build_variant_freq prefixes rare variants with 'rare'",
   all(grepl("^rare", vf$freq_label[vf$rare])))
ok("build_variant_freq formats the aFC label with a sign",
   any(grepl("aFC -6", vf$afc_label)))

# --- variant_effect_by_inheritance (Panels B/C core) ------------------------
vebi <- variant_effect_by_inheritance(eqtl, cg)
ok("variant_effect_by_inheritance assigns Control to low-shet non-ClinGen genes",
   "Control" %in% as.character(vebi$allelicRequirements[vebi$gene_name == "CTRL1"]))
ok("variant_effect_by_inheritance excludes high-shet non-ClinGen genes (NONCLIN)",
   !"NONCLIN" %in% vebi$gene_name)
ok("variant_effect_by_inheritance drops common variants (rare filter keeps aFC != -1)",
   !any(vebi$log2_aFC == -1))
ok("preserved fan-out: a gene curated under BOTH AD and AR counts under both",
   setequal(as.character(vebi$allelicRequirements[vebi$gene_name == "GENE_BOTH"]),
            c("AD", "AR")))

# --- panels build without error (smoke test on fixture frames) --------------
ok("build_panelB returns an untagged ggplot",
   inherits(build_panelB(eqtl, cg), "ggplot") && is.null(build_panelB(eqtl, cg)$labels$tag))
ok("build_panelC returns a ggplot", inherits(build_panelC(eqtl, cg), "ggplot"))
ok("build_panelA returns a patchwork", inherits(build_panelA(eqtl, cg), "patchwork"))

# --- pct_at_element (SPR % change at the regulatory element) ----------------
pc <- pct_at_element(fx$spr_tracks)
ok("pct_at_element returns one row per mark", nrow(pc) == 1L && pc$mark[1] == "CTCF")
ok("pct_at_element computes alt-vs-ref % change (alt = half of ref -> -50%)",
   abs(pc$pct[1] - (-50)) < 1e-6)

# --- id utilities -----------------------------------------------------------
ok("make_variant_id normalizes delimiters and upcases alleles",
   make_variant_id("chr1", 100, "a", "g") == "chr1_100_A_G")
ok("normalize_variant_id collapses repeated separators",
   normalize_variant_id("chr1::100**A") == "chr1_100_A")

cat("test_fig5_helpers.R passed.\n")

# =============================================================================
# Regional window + gene-model resolution (added with the GTF-backed gene track)
# =============================================================================

# regional_window_bp() is a pure seam; redefine nothing, just source the rule.
# Positions: a wide tested cis window, with the causal variants clustered.
pos    <- c(1e6, 1.4e6, 1.499e6, 1.5e6, 1.501e6, 1.6e6, 2e6)
causal <- c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)

w <- regional_window_bp(pos, causal, zoom_center = 1.5e6, zoom_window_bp = 200,
                        flank_bp = 5000)
ok("window is driven by the causal variants, not the tested range",
   w[1] >= 1.49e6 && w[2] <= 1.51e6)
ok("window is far tighter than the full tested span",
   diff(w) < diff(range(pos)) / 10)
ok("flank is applied around the causal span",
   isTRUE(all.equal(w, c(1.499e6 - 5000, 1.501e6 + 5000))))

# The zoom region must stay inside the panel, or the inset callout is orphaned.
wz <- regional_window_bp(pos, causal, zoom_center = 1.5e6, zoom_window_bp = 50000,
                         flank_bp = 100)
ok("window widens to contain the zoom region",
   wz[1] <= 1.5e6 - 50000 && wz[2] >= 1.5e6 + 50000)

# Never draw axis beyond what was actually tested.
wc <- regional_window_bp(c(1e6, 1.1e6), c(TRUE, TRUE), zoom_center = 1.05e6,
                         zoom_window_bp = 10, flank_bp = 5e6)
ok("window is clamped to the tested range",
   wc[1] == 1e6 && wc[2] == 1.1e6)

# No variant reaching the PIP threshold must not produce an empty/NA window.
wn <- regional_window_bp(pos, rep(FALSE, length(pos)), zoom_center = 1.5e6,
                         zoom_window_bp = 200, flank_bp = 5000)
ok("falls back to all tested variants when nothing is causal",
   all(is.finite(wn)) && diff(wn) > 0)

# --- gene-model resolution --------------------------------------------------
# Genes must be selected by chromosome AND coordinate. Before this change the
# track filtered on coordinate alone, which only worked because the two
# published loci do not overlap numerically.
gm <- gene_models_for_window("chr12", 51888000, 51930000)
ok("chr12 window resolves ANKRD33 + ACVRL1",
   setequal(gm$spans$gene, c("ANKRD33", "ACVRL1")))
ok("exons are restricted to the resolved genes",
   all(gm$exons$gene %in% gm$spans$gene))

spr_range <- gene_models_for_window("chr12", 72887408, 72892158)
ok("SPR's coordinates on the WRONG chromosome resolve to nothing",
   nrow(spr_range$spans) == 0)
ok("SPR resolves on its own chromosome",
   identical(gene_models_for_window("chr2", 72887408, 72892158)$spans$gene, "SPR"))

# --- inset clearance --------------------------------------------------------
# The zoom inset sits on an opaque background, so a bottom edge landing mid-marker
# slices the point in half and it reads as clipped. This is what the SPR panel
# used to do to its PIP-0.06 variant.
b_clear <- inset_bottom_clearing_points(x_npc = 0.5, y_npc = 0.9, requested = 0.05,
                                        inset_left = 0.24, inset_right = 0.96)
ok("a bottom with nothing near it is left alone", b_clear == 0.05)

b_one <- inset_bottom_clearing_points(x_npc = 0.5, y_npc = 0.056, requested = 0.05,
                                      inset_left = 0.24, inset_right = 0.96,
                                      marker_npc = 0.022)
ok("a marker straddling the edge pushes the inset above it",
   isTRUE(all.equal(b_one, 0.056 + 0.022)))

# Clearing one point must not park the edge on top of the next one up.
b_lad <- inset_bottom_clearing_points(x_npc = rep(0.5, 3), y_npc = c(0.056, 0.098, 0.110),
                                      requested = 0.05, inset_left = 0.24,
                                      inset_right = 0.96, marker_npc = 0.022)
ok("the edge climbs a ladder of points until every one is clear",
   all(abs(c(0.056, 0.098, 0.110) - b_lad) >= 0.022 - 1e-9))

# Points outside the inset's x-span are not occluded, so they must not move it.
b_out <- inset_bottom_clearing_points(x_npc = c(0.05, 0.99), y_npc = c(0.056, 0.06),
                                      requested = 0.05, inset_left = 0.24,
                                      inset_right = 0.96)
ok("points beside the inset leave its bottom edge alone", b_out == 0.05)

ok("the inset is never pushed past its ceiling",
   inset_bottom_clearing_points(x_npc = rep(0.5, 40), y_npc = seq(0.05, 0.9, length.out = 40),
                                requested = 0.05, inset_left = 0.24,
                                inset_right = 0.96) == 0.4)
