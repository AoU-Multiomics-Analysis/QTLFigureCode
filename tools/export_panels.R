# =============================================================================
# Export each figure's panels as standalone files.
#
# Every figure caches its panels as ggplot objects in <fig>/objects/*.rds, so
# this renders those directly rather than rebuilding anything. Panels must
# already be built (run the figure's run_all.R, or its individual panel
# scripts) — nothing here reads a controlled-access input.
#
# Run from the project root:
#   Rscript tools/export_panels.R              # every figure with built panels
#   Rscript tools/export_panels.R Fig1 Fig3    # just these
#
# Output lands in <fig>/objects/panels/<Fig><Letter>.pdf and .png, which the
# per-figure .gitignore already covers via `objects/`.
#
# Sizes come from PANEL_SIZE below and can be overridden for a run with the
# PANEL_W / PANEL_H environment variables (inches), e.g.
#   PANEL_W=8 PANEL_H=6 Rscript tools/export_panels.R Fig4
#
# Panels are exported WITHOUT their A/B/C tag: tagging happens in each figure's
# assemble.R and belongs to the assembled layout, not to a standalone panel.
# =============================================================================

FIGURES <- c("Fig1", "Fig2", "Fig3", "Fig4", "Fig5", "SuppFig1")

# Roughly the space each panel occupies in its assembled figure, in inches.
# (Assembled sizes: Fig1 25x27, Fig2 40x24, Fig3 18x12, Fig4 15x10, Fig5 18x16,
# SuppFig1 35x23.) These are starting points — override per run if a particular
# panel wants a different aspect.
PANEL_SIZE <- list(
  Fig1     = c(w = 12.5, h = 6.75),
  Fig2     = c(w = 13.0, h = 8.0),
  Fig3     = c(w =  9.0, h = 6.0),
  Fig4     = c(w =  7.5, h = 5.0),
  Fig5     = c(w =  9.0, h = 8.0),
  SuppFig1 = c(w = 17.0, h = 11.0)
)

# Per-panel footprints, for panels whose share of the assembled figure is not
# the figure default. Text is sized for the assembled figure, so a panel
# exported much narrower than its real footprint collides its axis labels —
# which is exactly what a per-figure default does to a wide panel.
#
# Fig1 is 25 x 27 over four equal rows (6.75 in each). Rows 1-3 split evenly at
# 12.5 in; row 4 is `plot_layout(widths = c(1, 2.5))`, so G takes 25/3.5 and H
# takes 25 * 2.5/3.5. The other figures fall back to their default until their
# layouts are encoded here the same way.
PANEL_OVERRIDE <- list(
  Fig1G = c(w =  7.14, h = 6.75),
  Fig1H = c(w = 17.86, h = 6.75)
)

# --- one figure per child process --------------------------------------------
# Each figure's 00_setup.R sets the global ggplot theme, the font and DPI, and
# rebinds OBJ_DIR. Sourcing several of them into one session would leave
# whichever ran last in charge, so a figure that is exported second would be
# drawn with the previous figure's theme. A child process per figure keeps them
# honest.
args <- commandArgs(trailingOnly = TRUE)
targets <- if (length(args)) args else FIGURES
unknown <- setdiff(targets, FIGURES)
if (length(unknown))
  stop("unknown figure(s): ", paste(unknown, collapse = ", "),
       "; choose from ", paste(FIGURES, collapse = ", "))

if (is.na(Sys.getenv("EXPORT_PANELS_CHILD", unset = NA))) {
  failed <- character(0)
  for (fig in targets) {
    message("\n== ", fig)
    status <- system2("Rscript", c("tools/export_panels.R", fig),
                      env = "EXPORT_PANELS_CHILD=1")
    if (!identical(status, 0L)) failed <- c(failed, fig)
  }
  if (length(failed)) {
    message("\nFAILED: ", paste(failed, collapse = ", "))
    quit(status = 1)
  }
  quit(status = 0)
}

# --- child: export one figure -------------------------------------------------
fig   <- targets[1]
lower <- tolower(fig)
suppressMessages(source(file.path(lower, "00_setup.R")))

rds <- sort(list.files(OBJ_DIR, pattern = "^[A-Z][A-Za-z0-9_]*\\.rds$", full.names = TRUE))
if (!length(rds)) {
  message("  no built panels in ", OBJ_DIR, " — run ", lower, "/run_all.R first")
  quit(status = 0)
}

out_dir <- file.path(OBJ_DIR, "panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (f in rds) {
  letter <- tools::file_path_sans_ext(basename(f))
  p      <- readRDS(f)
  name   <- paste0(fig, letter)
  stem   <- file.path(out_dir, name)

  size <- PANEL_OVERRIDE[[name]]
  if (is.null(size)) size <- PANEL_SIZE[[fig]]
  # An explicit PANEL_W / PANEL_H applies to every panel in the run and wins
  # over both the override and the default.
  w <- as.numeric(Sys.getenv("PANEL_W", size[["w"]]))
  h <- as.numeric(Sys.getenv("PANEL_H", size[["h"]]))

  ggsave(paste0(stem, ".pdf"), p, width = w, height = h, dpi = DPI, limitsize = FALSE)
  # Raster preview at full physical size and a lower DPI, with showtext's DPI
  # matched — the same idiom the assemble.R scripts use, and required, or the
  # text is scaled for 300 dpi and comes out enormous.
  showtext_opts(dpi = 100)
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 100, limitsize = FALSE)
  showtext_opts(dpi = DPI)

  message("  ", basename(stem), ".pdf/.png  (", w, " x ", h, " in)")
}

message("  ", length(rds), " panels -> ", out_dir)
