suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
  library(ggh4x)
  library(patchwork)
  library(MatchIt)
  library(readxl)
  library(showtext)
  library(sysfonts)    
  library(paletteer)
})

font_add_google("Inter", "inter")
showtext_opts(dpi = 400)
options(
  repr.plot.width = 15,
  repr.plot.height = 10,
  repr.plot.res = 400
)



data_path <- 'gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/submissions/1a4f6a11-3ec4-48ee-86b7-e9eace3d2e4a/RunFastenloc/86d1e706-b60a-4d80-bcdc-2a66b5ea59d1/call-CollectHighLevelOutputs/high_level_coloc_outputs.tar.gz'
#system(paste0('gsutil cp ',data_path , ' . '))
#system(paste0('tar -xvf ',basename(data_path)))

GROUP_LABELS <- c(
  autoimmune = "AutoImm",
  inflammatory_biomarker = "Inflam",
  allergic_respiratory = "Allerg",
  dermatological = "Derm",
  cardiometabolic = "CardioM",
  lipid_metabolic = "Lipid",
  lipid = "Lipid",
  metabolic = "Metabolic",
  clinical_chemistry = "ClinChem",
  hematological = "Heme",
  anthropometric = "Anthro",
  renal = "Renal",
  renal_urological = "Renal/Urol",
  gastrointestinal = "GI",
  reproductive_endocrine = "ReproEnd",
  musculoskeletal = "MSK",
  cancer = "Cancer",
  neurodegenerative = "Neurodeg",
  neurological = "Neurol",
  neuropsychiatric = "NeuroPsy",
  ophthalmological = "Ophth"
)

GROUP_COLORS <- c(
  autoimmune = "#D62828",
  inflammatory_biomarker = "#F77F00",
  allergic_respiratory = "#2A9D8F",
  dermatological = "#E76F51",
  cardiometabolic = "#277DA1",
  lipid_metabolic = "#6A994E",
  lipid = "#6A994E",
  metabolic = "#80B918",
  clinical_chemistry = "#F4A261",
  hematological = "#A23B72",
  anthropometric = "#E9C46A",
  renal = "#457B9D",
  renal_urological = "#8D99AE",
  gastrointestinal = "#BC6C25",
  reproductive_endocrine = "#E5989B",
  musculoskeletal = "#7F5539",
  cancer = "#9B5DE5",
  neurodegenerative = "#4361EE",
  neurological = "#3A0CA3",
  neuropsychiatric = "#B5179E",
  ophthalmological = "#48CAE4"
)

RATE_COLUMNS <- c(
  stringent = "pct_stringent_union",
  lenient = "pct_lenient_union"
)


clingen_data <- arrow::read_parquet('evidence_clingen/part-00000-006c7897-e219-4036-abbf-af4482a5177e-c000.snappy.parquet') %>%
        select(targetFromSourceId,allelicRequirements,confidence,score,studyId,targetId,allelicRequirements,diseaseFromSource) 
clingen_cleaned_list <- clingen_data %>% 
    distinct(targetId,targetFromSourceId,studyId,confidence,allelicRequirements,diseaseFromSource) %>%
    filter(confidence == 'Definitive' | confidence == 'Strong') %>% 
    unnest(allelicRequirements)
dg <- fread("/home/jupyter/interactions.tsv")
dg[,sym:=toupper(gene_name)]
druggable <- unique(dg[sym!="",sym]); approved <- unique(dg[approved==TRUE&sym!="",sym])
GeneBayes <- read_excel("41588_2024_1820_MOESM4_ESM.xlsx",
                        sheet = "Supplementary Table 1") %>%
  transmute(gene = ensg, shet = as.numeric(post_mean),log_shet = log10(shet)) %>%
  filter(!is.na(shet), shet > 0)

heme_immune_kw <- str_c(c(
  # hematologic
  "an[a]?emia","thrombocyt","neutropen","leuk[a]?emia","lymphoma","h[a]?emoglobin","platelet",
  "h[a]?emolytic","spherocytosis","sideroblastic","myelodysplas","granulocyt","erythrocyt","reticulocyt",
  "coagulation","thromb","bleeding","h[a]?emophilia","neutrophil","eosinophil","basophil","monocyt",
  "myeloid","h[a]?ematopoie","pancytopenia","h[a]?emochromatosis","polycythemia","thalass","macrocyt",
  "blood cell","bone marrow","hemostasis","porphyria","fanconi",
  # immune
  "immunodefic","immune","autoimmun","autoinflammat","immunodysreg","agammaglobulin","lymphoproliferat",
  "complement","immunoglobulin","hypogammaglobulin","granulomatous","histiocyt","hemophagocytic","interferonopathy"),
  collapse = "|")
heme_immune <- clingen_cleaned_list %>%
  filter(str_detect(tolower(diseaseFromSource), heme_immune_kw))   # <- adjust col name if different

heme_immune_genes <- unique(heme_immune$targetId)

read_any <- function(path, ...) {
  args <- list(showProgress = FALSE, ...)
  if (grepl("\\.gz$", path)) args$cmd <- paste("gzip -dc", shQuote(path))
  else                         args$input <- path
  as_tibble(do.call(data.table::fread, args))
}

require_columns <- function(tbl, cols, label) {
  missing <- setdiff(cols, names(tbl))
  if (length(missing) > 0) {
    stop(label, " is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}

as_logical_flag <- function(x) {
  if (is.logical(x)) return(replace_na(x, FALSE))
  if (is.numeric(x)) return(replace_na(x != 0, FALSE))
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

complete_category_colors <- function(categories) {
  categories <- sort(unique(as.character(categories)))
  missing <- setdiff(categories, names(GROUP_COLORS))
  fallback <- if (length(missing) > 0) {
    setNames(scales::hue_pal(l = 55, c = 90)(length(missing)), missing)
  } else {
    character()
  }
  c(GROUP_COLORS, fallback)
}

blend_with_white <- function(hex, amount = 0.65) {
  x <- grDevices::col2rgb(hex) / 255
  y <- x * (1 - amount) + amount
  grDevices::rgb(y[1, ], y[2, ], y[3, ])
}

make_fill <- function(category, intensity, colors) {
  base <- unname(colors[as.character(category)])
  light <- blend_with_white(base, 0.72)
  vapply(seq_along(base), function(i) {
    palette <- grDevices::colorRampPalette(c(light[[i]], base[[i]]))(100)
    palette[pmax(1, pmin(100, round(intensity[[i]] * 99) + 1))]
  }, character(1))
}

reorder_within <- function(x, by, within, fun = mean, sep = "___") {
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  scale_y_discrete(..., labels = function(x) gsub(reg, "", x))
}

build_plot_data <- function(gene, coloc_rate, min_coloc_genes = 1L,
                            rate_type = "stringent") {
  rate_type <- tolower(trimws(rate_type))
  if (length(rate_type) != 1 || is.na(rate_type) ||
      !rate_type %in% names(RATE_COLUMNS)) {
    stop("--rate_type must be one of: ",
         paste(names(RATE_COLUMNS), collapse = ", "), call. = FALSE)
  }
  rate_column <- unname(RATE_COLUMNS[[rate_type]])
  require_columns(gene, c("gene", "trait", "trait_category", "any_coloc"), "--gene")
  require_columns(
    coloc_rate,
    c("trait", "layer", unname(RATE_COLUMNS)),
    "--coloc_rate"
  )
  if (length(min_coloc_genes) != 1 || is.na(min_coloc_genes) ||
      min_coloc_genes < 1 || min_coloc_genes != as.integer(min_coloc_genes)) {
    stop("--min_coloc_genes must be a positive integer.", call. = FALSE)
  }

  gene$any_coloc <- as_logical_flag(gene$any_coloc)
  coloc_counts <- gene %>%
    filter(any_coloc, !is.na(gene), gene != "", !is.na(trait_category), trait_category != "") %>%
    distinct(gene, trait, trait_category) %>%
    count(trait, trait_category, name = "count")

  pct_coloc <- coloc_rate %>%
    filter(layer == "union") %>%
    transmute(
      trait,
      pct_stringent_union,
      pct_lenient_union,
      pct_coloc_union = .data[[rate_column]]
    )
  duplicated_traits <- pct_coloc %>% count(trait) %>% filter(n > 1) %>% pull(trait)
  if (length(duplicated_traits) > 0) {
    stop("--coloc_rate has multiple union rows for trait(s): ",
         paste(duplicated_traits, collapse = ", "), call. = FALSE)
  }

  combined <- coloc_counts %>% left_join(pct_coloc, by = "trait")
  missing_rate <- combined %>% filter(is.na(pct_coloc_union)) %>% pull(trait)
  if (length(missing_rate) > 0) {
    message("Dropped (no union ", rate_type, " coloc rate): ",
            paste(unique(missing_rate), collapse = ", "))
  }
  combined <- combined %>%
    filter(!is.na(pct_coloc_union)) %>%
    mutate(rate_type = .env$rate_type)
  trimmed <- combined %>% filter(count < min_coloc_genes) %>% pull(trait)
  if (length(trimmed) > 0) {
    message("Trimmed ", length(trimmed), " trait(s) with fewer than ",
            min_coloc_genes, " colocalizing genes.")
  }
  combined %>% filter(count >= min_coloc_genes)
}

build_empty_plot <- function() {
  ggplot() +
    annotate(
      "text", x = 0, y = 0,
      label = "No traits have both colocalizing genes and the selected union coloc rate",
      color = "grey35", size = 4
    ) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void()
}

build_plot <- function(combined) {
  colors <- complete_category_colors(combined$trait_category)
  plotdat <- combined %>%
    group_by(trait_category) %>%
    mutate(
      intensity = if (max(count) == min(count)) {
        rep(0.75, n())
      } else {
        (count - min(count)) / (max(count) - min(count))
      },
      intensity = 0.25 + 0.75 * intensity,
      fill_col = make_fill(trait_category, intensity, colors)
    ) %>%
    ungroup() %>%
    mutate(
      trait_category = forcats::fct_drop(
        forcats::fct_reorder(trait_category, count, .fun = sum)
      ),
      trait = reorder_within(trait, count, trait_category),
      base_col = unname(colors[as.character(trait_category)])
    )

  long <- bind_rows(
    plotdat %>% transmute(trait_category, trait, base_col, fill_col,
                          metric = "a_count", value = count),
    plotdat %>% transmute(trait_category, trait, base_col, fill_col,
                          metric = "b_rate", value = pct_coloc_union)
  )
  xline <- -max(plotdat$count) * 0.03
  seg <- plotdat %>% count(trait_category, name = "nbar") %>% mutate(metric = "a_count")
  vref <- tibble(metric = "b_rate", xi = median(plotdat$pct_coloc_union, na.rm = TRUE))
  selected_rate <- unique(plotdat$rate_type)
  metric_labels <- c(
    a_count = "Colocalizing Genes",
    b_rate = paste0("CS coloc rate ")
  )
  category_labeller <- c(GROUP_LABELS, setNames(
    setdiff(as.character(unique(plotdat$trait_category)), names(GROUP_LABELS)),
    setdiff(as.character(unique(plotdat$trait_category)), names(GROUP_LABELS))
  ))

  ggplot(long, aes(y = trait)) +
    geom_col(data = ~ subset(.x, metric == "a_count"),
             aes(x = value, fill = fill_col), width = 0.85, show.legend = FALSE) +
    geom_segment(data = seg, inherit.aes = FALSE,
                 aes(y = 0.5, yend = nbar + 0.5, x = xline, xend = xline),
                 color = "grey60", linewidth = 0.5) +
    geom_vline(data = vref, aes(xintercept = xi),
               linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_segment(data = ~ subset(.x, metric == "b_rate"),
                 aes(x = 0, xend = value, yend = trait, color = base_col), linewidth = 0.4) +
    geom_point(data = ~ subset(.x, metric == "b_rate"),
               aes(x = value, color = base_col), size = 1.6) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_y_reordered() +
    facet_grid(
      trait_category ~ metric,
      scales = "free",
      space = "free_y",
      switch = "y",
      labeller = labeller(
        trait_category = as_labeller(category_labeller, default = label_value),
        metric = as_labeller(metric_labels)
      )
    ) +
    ggh4x::facetted_pos_scales(x = list(
      scale_x_continuous(),
      scale_x_continuous(
        labels = function(x) paste0(x, "%"),
        expand = expansion(mult = c(0, 0.08))
      )
    )) +
    ggh4x::force_panelsizes(cols = c(3, 1.3)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.spacing.x = grid::unit(8, "pt"),
      panel.spacing.y = grid::unit(2, "pt"),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 12),
      strip.text.x = element_text(size = 12)
    )
}


gene_dat <- fread('high_level_coloc_outputs/harmonized/harmonized_coloc.gene.tsv.gz')

combined <- build_plot_data(
    read_any('high_level_coloc_outputs/harmonized/harmonized_coloc.gene.tsv.gz'), 
    read_any('high_level_coloc_outputs/final_summaries/coloc_rate_by_trait.tsv'),
    min_coloc_genes = 1, rate_type ='lenient'
  )
gene_data <- fread('high_level_coloc_outputs/harmonized/harmonized_coloc.gene.tsv.gz')
coloc_rate <- fread('high_level_coloc_outputs/final_summaries/coloc_rate_by_trait.tsv')
coloc_rate %>% arrange(desc(pct_lenient_union)) %>% filter(layer == 'union')
coloc_plot <- combined %>% 
    filter(trait_category != 'neuropsychiatric' & trait_category != 'renal_urological') %>% 
    filter(count > 5) %>% 
    group_by(trait_category) %>% 
    arrange(desc(count)) %>% 
    filter(row_number() < 11) %>% 
    build_plot +
    theme(axis.text.x = element_text(size = 10))
coloc_plot
R2_combined <- fread('R2_TWAS_combined.tsv') %>% 
    mutate(H2_upper = (SE *1.96) + H2) %>% 
    mutate(adj_R2 = r.squared/H2)  %>% 
    filter(Pval < 1e-5) %>% 
    select(phenotype_id,AF_threshold,adj_R2,H2,Pval) %>%
    pivot_wider(names_from = AF_threshold,values_from = adj_R2) %>% 
    mutate(difference = (`0` - `0.01`))
           
R2_dist <- R2_combined %>% 
    #filter(!str_detect(gene_type,'pseudo')) %>% 
    ggplot(aes(x = 100*difference)) + 
        geom_histogram(color = 'white') + 
        ylab('Number of genes') + 
        xlim(-20,40) + 
        theme_classic() +
        theme(axis.text = element_text(size =12),
              axis.title = element_text(size = 14)
             ) + 
        labs(x = expression(R^2 ~ increase ~ "(%)"))
R2_dist
blank_plot <- ggplot() + theme_void()

ConditionalTWAS <- fread('ConditionalTWAS.cleaned.tsv.gz')
prepare_qq_match_source <- function(data) {
  data %>%
    filter(status == "ok") %>%
    distinct(
      gene,
      coloc_trait,
      molecular_trait_id,
      .keep_all = TRUE
    ) %>%
    mutate(
      any_coloc = coalesce(any_coloc, FALSE),
      inflation = abs(z_rare_cond) /
        pmax(abs(rare_zscore), 0.5),
      artifact = coalesce(
        inflation > 3 & abs(Vcor) > 0.3,
        FALSE
      )
    ) %>%
    filter(!artifact) %>%
    group_by(gene, coloc_trait) %>%
    arrange(
      desc(coalesce(GLCP, -Inf)),
      p_rare_cond,
      .by_group = TRUE
    ) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    mutate(
      coloc = as.integer(any_coloc),
      log_tpm = log1p(mean_TPM),
      log_shet = log10(shet),
      log_n_rare = log1p(n_rare),
      abs_vcor = abs(Vcor)
    ) %>%
    filter(
      is.finite(log_tpm),
      is.finite(log_shet),
      is.finite(log_n_rare),
      is.finite(abs_vcor)
    )
}

calculate_qq_data <- function(
    matched_data,
    traits,
    common_cap = 40,
    rare_cap = 12.5
) {
  matched_data %>%
    filter(as.character(coloc_trait) %in% traits) %>%
    select(
      group,
      common = common_pvalue,
      rare = p_rare_cond
    ) %>%
    pivot_longer(
      cols = c(common, rare),
      names_to = "test",
      values_to = "p_value"
    ) %>%
    filter(
      is.finite(p_value),
      between(p_value, 0, 1)
    ) %>%
    group_by(test, group) %>%
    arrange(p_value, .by_group = TRUE) %>%
    mutate(
      rank = row_number(),
      n = n(),
      expected = -log10(rank / (n + 1)),
      observed = -log10(
        pmax(p_value, .Machine$double.xmin)
      ),
      y_cap = case_when(
        test == "common" ~ common_cap,
        test == "rare" ~ rare_cap,
        .default = Inf
      ),
      observed_display = pmin(observed, y_cap),
      capped = observed > y_cap
    ) %>%
    ungroup()
}

match_colocalized_pairs <- function(
    data,
    seed = 20260713,
    ratio = 1,
    caliper = 0.2,
    replace = FALSE
) {
  set.seed(seed)

  match_fit <- MatchIt::matchit(
    coloc ~ log_tpm + log_shet + log_n_rare + abs_vcor,
    data = data,
    method = "nearest",
    distance = "glm",
    exact = ~ coloc_trait,
    ratio = ratio,
    replace = replace,
    caliper = caliper,
    std.caliper = TRUE
  )

  matched_data <- MatchIt::match.data(
    match_fit,
    drop.unmatched = TRUE
  ) %>%
    mutate(
      group = if_else(
        coloc == 1L,
        "Colocalized genes",
        "Matched non-colocalized genes"
      ),
      group = factor(
        group,
        levels = c(
          "Matched non-colocalized genes",
          "Colocalized genes"
        )
      )
    )

  list(
    fit = match_fit,
    data = matched_data
  )
}
plot_qq <- function(qq_data) {
  ggplot(
    qq_data,
    aes(
      x = expected,
      y = observed_display,
      colour = group
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.55,
      colour = "#666666"
    ) +
    geom_point(
      size = 1,
      alpha = 0.8
    ) +
    facet_wrap(
      ~test,
      scales = "free_y",
      labeller = as_labeller(
        c(
          common = "Common",
          rare = "Rare"
        )
      )
    ) +
    scale_colour_manual(
      values = c(
        "Matched non-colocalized genes" = "#2B2B2B",
        "Colocalized genes" = "#275D8E"
      ),
      name = NULL
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(
          shape = 16,
          size = 2.5,
          alpha = 1
        )
      )
    ) +
    labs(
      x = expression(Expected ~ -log[10](P)),
      y = expression(Observed ~ -log[10](P))
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      panel.spacing = grid::unit(2, "lines"),
      legend.position = "top",
      legend.text = element_text(size = 10),
      legend.key.width = grid::unit(1.2, "lines")
    )
}


blood_traits <- c(
  "eosinophil",
  "hematocrit",
  "lymphocyte",
  "MCH",
  "MCHC",
  "monocyte",
  "neutrophil",
  "platelet_count",
  "RDW"
)
qq_match_source <- prepare_qq_match_source(ConditionalTWAS)

match_result <- match_colocalized_pairs(qq_match_source)

match_fit <- match_result$fit
matched_qq <- match_result$data

qq_blood <- calculate_qq_data(
  matched_data = matched_qq,
  traits = blood_traits
)

p_blood_qq <- qq_blood %>% 
    mutate(group = case_when(str_detect(group,'Matched') ~ 'Null',TRUE ~ group)) %>% 
    plot_qq() +
    theme(legend.position = 'none')
p_blood_qq

ColocAdjustedConditionalTWAS <- ConditionalTWAS %>%  
    distinct(gene, coloc_trait, molecular_trait_id, .keep_all = TRUE) %>% 
    mutate(inflation = abs(z_rare_cond) / pmax(abs(rare_zscore), 0.5),
             artifact  = inflation > 3 & abs(Vcor) > 0.3) %>%  
    filter(!artifact)  %>%
    distinct() %>%
    filter(any_coloc == TRUE) %>% 
    group_by(TWAS_path) %>% 
    mutate(fdr_rare_cond = p.adjust(p_rare_cond,method= 'fdr'), 
           fdr_common = p.adjust(common_pvalue,method = 'fdr'),
           rare_zscore_fdr = p.adjust(rare_pvalue,method = 'fdr'),
           fdr_full = p.adjust(full_pvalue,method = 'fdr'),
           artifact = abs(z_rare_cond)/abs(rare_zscore) > 4 & abs(Vcor) > 0.3
          ) %>% 
    mutate(genuine_rare = p_rare_cond < .05 & rare_pvalue < .05,
           t_nominal= p_rare_cond < .05 & rare_pvalue < .05,
           t_fdrcond = fdr_rare_cond < .1 & rare_pvalue < .05,
           t_fdrboth = fdr_rare_cond < .1 & fdr_rare_sig, 
           gene = factor(gene), 
           coloc_trait = factor(coloc_trait),
           druggable = toupper(Description) %in% druggable,
           approved = toupper(Description) %in% approved
            ) 

GeneTWASData <- ColocAdjustedConditionalTWAS %>%  
  group_by(gene, Description) %>%
  summarise(max_rare   = max(abs(z_rare_cond),   na.rm=TRUE),    # strongest independent rare signal
            max_common = max(abs(common_zscore), na.rm=TRUE),
            TPM        = mean(mean_TPM, na.rm=TRUE),
            druggable = first(druggable),
            approved = first(approved),
            ever_common = any(fdr_common < .1,na.rm = T),
            ever_rare = any(t_fdrcond),
            shet = first(shet),
            .groups="drop") %>%
  mutate(log_shet=log10(shet), 
         log_tpm=log10(TPM+1),
         MendelianGene = gene %in% clingen_cleaned_list$targetId,
         HemeImmuneGene = gene %in% heme_immune_genes
        ) %>%
  filter(!is.na(log_shet)) 

# Avoid ambiguity between the external symbol vectors and data columns.
druggable_symbols <- druggable
approved_symbols <- approved

ColocAdjustedConditionalTWAS <- ConditionalTWAS %>%
  distinct(gene, coloc_trait, molecular_trait_id, .keep_all = TRUE) %>%
  mutate(
    inflation = abs(z_rare_cond) / pmax(abs(rare_zscore), 0.5),
    artifact_primary = coalesce(
      inflation > 3 & abs(Vcor) > 0.3,
      FALSE
    )
  ) %>%
  filter(!artifact_primary) %>%
  filter(any_coloc == TRUE) %>%
  group_by(TWAS_path) %>%
  mutate(
    fdr_rare_cond = p.adjust(p_rare_cond, method = "BH"),
    fdr_common = p.adjust(common_pvalue, method = "BH"),
    rare_zscore_fdr = p.adjust(rare_pvalue, method = "BH"),
    fdr_full = p.adjust(full_pvalue, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    genuine_rare = p_rare_cond < 0.05 & rare_pvalue < 0.05,
    t_nominal = p_rare_cond < 0.05 & rare_pvalue < 0.05,
    t_fdrcond = fdr_rare_cond < 0.10 & rare_pvalue < 0.05,
    t_fdrboth = fdr_rare_cond < 0.10 & rare_zscore_fdr < 0.10,

    gene = factor(gene),
    coloc_trait = factor(coloc_trait),
    druggable = toupper(as.character(Description)) %in% druggable_symbols,
    approved = toupper(as.character(Description)) %in% approved_symbols
  )

# Confirm no rows that meet the artifact rule remain.
stopifnot(
  nrow(
    ColocAdjustedConditionalTWAS %>%
      filter(abs(Vcor) > 0.3, inflation > 3)
  ) == 0
)
# ── fit + tidy one enrichment glm ─────────────────────────────────────────────
fit_enrich <- function(outcome, predictors, data = GeneTWASData) {
  glm(reformulate(predictors, response = outcome),
      data = data, family = binomial()) %>%
    broom::tidy(conf.int = TRUE, exponentiate = TRUE)
}

# ── define outcomes and predictor sets once ───────────────────────────────────
outcomes <- c(approved      = "Approved Drug Target",
              druggable     = "Druggable gene",
              MendelianGene = "Mendelian Gene\n(Any)",
              HemeImmuneGene = "Mendelian Gene\n(Immune/Hematological)")

predictor_sets <- list(
  continuous = c("max_rare", "max_common", "log_shet", "log2(TPM)"),
  binary     = c("ever_rare", "ever_common", "log_shet", "log2(TPM)")
)

# ── run every outcome x predictor-set combination ─────────────────────────────
enrichment_results <- expand_grid(outcome = names(outcomes),
                                  model   = names(predictor_sets)) %>%
  mutate(res = map2(outcome, model, ~ fit_enrich(.x, predictor_sets[[.y]]))) %>%
  unnest(res) %>%
  mutate(type = outcomes[outcome]) %>%
  relocate(type, model)

# just the rare-signal effect (the term you care about):
plot_df <- enrichment_results %>%
  filter(term %in% c("max_rare",'max_common')) %>%
  select(type, model, term, OR = estimate, conf.low, conf.high, p.value)  %>% 
  mutate(term = str_remove_all(term, "max_|ever_|TRUE"),
         term = recode(term, rare="Rare", common="Common"),
         term = factor(term, levels=c("Rare","Common")),
         type = factor(type, levels=c("Mendelian Gene\n(Any)",
                                      "Mendelian Gene\n(Immune/Hematological)",
                                      "Druggable gene",
                                      "Approved Drug Target"))
         )


pal <- c(Rare = as.character(paletteer_d("futurevisions::mars"))[1], Common = "grey60")
p <- ggplot(plot_df, aes(OR, reorder(type,OR), colour=term)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey60") +
  geom_pointrange(aes(xmin=conf.low, xmax=conf.high), position=position_dodge(width=.55),
                  linewidth=.8, fatten=2.6) +
  #facet_wrap(~model, scales="free_x") +
  scale_colour_manual(values=pal, name=NULL) +
  theme_classic(base_size=12) +
  theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), strip.background=element_blank(),
        strip.text=element_text(face="bold"), 
        legend.position= c(0.9,0.2),
        panel.spacing=unit(1.4,"lines")) +
    xlab('Odds Ratio (per unit increase in Z-score)') +
    ylab("") 
p

top_row <- (
  (R2_dist + labs(tag = "B")) |
  (p_blood_qq + labs(tag = "C"))
) +
  plot_layout(widths = c(0.8, 1.2))

# Isolate the top row from D's axis-label spacing
top_row <- wrap_elements(full = top_row)

right_side <- top_row /
  (p + labs(tag = "D")) +
  plot_layout(heights = c(1, 0.75))

final_plot <- (
  (
    coloc_plot +
      labs(tag = "A") +
      theme(
        plot.margin = margin(r = 55, unit = "pt"),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 13)
      )
  ) |
  right_side
) +
  plot_layout(widths = c(7, 8)) &
  theme(
    plot.tag = element_text(size = 14),
    plot.tag.position = c(0, 1)
  )

final_plot
ggsave('Fig4.pdf',width = 15,height = 10,dpi = 400 )













