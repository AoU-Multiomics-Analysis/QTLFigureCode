library(tidyverse)
library(data.table)
library(showtext)
library(sysfonts)
library(paletteer)
library(ggdist)
library(ggh4x)
library(rlang)
library(patchwork)
library(readxl)
library(ggrepel)
library(igraph)

font_add_google("Inter", "inter")
showtext_opts(dpi = 400)

showtext_auto()
theme_set(theme_classic())
theme_update( 
            strip.text = element_text(size = 10),
            text = element_text(family = "inter"),
            strip.background = element_blank()
            )

genebayes <- read_xlsx('/home/jupyter/41588_2024_1820_MOESM4_ESM.xlsx',sheet = 2) %>% 
  mutate(shet_q = as.factor(ntile(post_mean, 5))) %>% 
  mutate(
    shet_bin = cut(
      post_mean,
      breaks = c(0, 0.001,0.01, 0.05,1),
      include.lowest = TRUE,
      labels = c('No Selection','Moderate Selection','Strong Selection','Extreme Selection')
    )
  )
system(paste0('gsutil cp gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/LOD_proportions_per_proteins.tsv .'))

LOD_data <- fread('LOD_proportions_per_proteins.tsv')
LOD_data %>% head
median_tpms <- fread('median_tpms.tsv.gz') %>% mutate(Gene = str_remove(Name,'\\..*') ) %>% select(Gene,mean_TPM)
median_tpms %>% head

LoadEnlocGeneLevel <- function(Path) {
x <- readLines(Path)
x <- x[-1]  # drop header

mat <- do.call(
  rbind,
  strsplit(trimws(x), "\\s+")
)

out <- as.data.table(mat)
setnames(out, c(
  "Signal",
  "Num_SNP",
  "CPIP_qtl",
  "CPIP_gwas_marginal",
  "CPIP_gwas_qtl_prior",
  "RCP",
  "LCP"
    )
)
 
}
LoadEnlocVariantLevel <- function(Path) {
x <- readLines(Path)
x <- x[-1]  # drop header

mat <- do.call(
  rbind,
  strsplit(trimws(x), "\\s+")
)

out <- as.data.table(mat)
setnames(out, c(
  "Signal",
  "SNP",
  "PIP_qtl",
  "PIP_gwas_marginal",
  "PIP_gwas_qtl_prior",
  "SCP"
    )
)
}

LoadEnlocCSLevel <- function(Path) {
x <- readLines(Path)
x <- x[-1]  # drop header
mat <- do.call(
  rbind,
  strsplit(trimws(x), "\\s+")
)

out <- as.data.table(mat)
setnames(out, c(
  "Signal",
  "Num_SNP",
  "CPIP_qtl",
  "CPIP_gwas_marginal",
  "CPIP_gwas_qtl_prior",
  "RCP",
  "LCP"
    )
)
}

add_fdr_threshold_label <- function(df, prob_col, control_level = 0.05,
                                    label_col = "Coloc") {
  
  prob_col <- ensym(prob_col)
  
  probs <- df %>%
    pull(!!prob_col) %>%
    as.numeric()
  
  # match manual logic, but remove NAs
  lfdr <- sort(1 - probs, na.last = NA)
  FDR <- cumsum(lfdr) / seq_along(lfdr)
  
  valid <- which(FDR <= control_level)
  
  if (length(valid) == 0) {
    thresh <- Inf
  } else {
    thresh <- 1 - lfdr[max(valid)]
  }
  
  out <- df %>%
    mutate(
      !!label_col := as.numeric(!!prob_col) >= thresh
    )
  
  attr(out, "threshold") <- thresh
  
  out
}



AnnotationPath <- 'gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/checkpoints/AoU.COMB.AC10.AN95.biallelic.v9.annotations.tsv.bgz'
system(paste0('gsutil cp ',AnnotationPath,' .'))
AnnotationData <- fread(basename(AnnotationPath)) %>% 
    mutate(SNP = paste(chrom,pos,ref,alt,sep = '_'))

sQTL_COMB_susie <- 'sQTL_susie.COMB.annotated.cleaned.tsv.gz'
eQTL_COMB_susie <- 'eQTL_susie.COMB.annotated.cleaned.tsv.gz'
pQTL_COMB_susie <- 'pQTL_susie.COMB.annotated.cleaned.tsv.gz'

sQTL <- fread(basename(sQTL_COMB_susie)) %>% 
    mutate(molecular_trait_id = stringr::str_extract(molecular_trait_id, "ENSG[0-9]+")) %>% 
    select(molecular_trait_id,variant,pip,cs_index,posterior_mean,cs_id,posterior_sd) %>% 
    group_by(molecular_trait_id,cs_index) %>% 
    mutate(cpip = sum(pip),number_variants = dplyr::n()) %>%
    ungroup()
eQTL <- fread(basename(eQTL_COMB_susie)) %>% 
    select(molecular_trait_id,variant,pip,cs_index,posterior_mean,cs_id,posterior_sd,gene_name) %>% 
    mutate(molecular_trait_id = str_remove(molecular_trait_id,'\\..*')) %>%
    mutate(variant = str_replace(variant,'chrchr','chr')) %>% 
    group_by(molecular_trait_id,cs_index) %>% 
    mutate(cpip = sum(pip),number_variants = dplyr::n())  %>% 
    ungroup()
pQTL <- fread(basename(pQTL_COMB_susie)) %>% 
    filter(group == 'COMB') %>% 
    mutate(molecular_trait_id = stringr::str_extract(molecular_trait_id, "ENSG[0-9]+")) %>% 
    select(molecular_trait_id,variant,pip,cs_index,posterior_mean,cs_id,posterior_sd) %>% 
    group_by(molecular_trait_id,cs_index) %>% 
    mutate(cpip = sum(pip),number_variants = dplyr::n()) %>%
    ungroup()
eQTL_cs_list <- eQTL %>% distinct(molecular_trait_id,cs_id)
sQTL_cs_list <- sQTL %>% distinct(molecular_trait_id,cs_id)
pQTL_cs_list <- pQTL %>% distinct(molecular_trait_id,cs_id)
cs_list <- bind_rows(eQTL_cs_list %>% mutate(type = 'eQTL'),
                     pQTL_cs_list %>% mutate(type = 'pQTL'),
                     sQTL_cs_list %>% mutate(type = 'sQTL')
                    )

cs_count <- cs_list %>% 
                group_by(molecular_trait_id,type) %>% 
                summarize(num_cs = dplyr::n()) %>% 
                pivot_wider(names_from = type,values_from = num_cs,values_fill  = 0)
n_eQTLs <- eQTL %>% 
    group_by(molecular_trait_id) %>% 
    distinct(cs_id) %>% 
    summarize(n_CS = dplyr::n())
n_sQTLs <- sQTL %>% 
    group_by(molecular_trait_id) %>% 
    distinct(cs_id) %>% 
    summarize(n_CS = dplyr::n())
n_pQTLs <- pQTL %>% 
    group_by(molecular_trait_id) %>% 
    distinct(cs_id) %>% 
    summarize(n_CS = dplyr::n())

QTL_union_gene_list <- eQTL %>% 
    filter(molecular_trait_id %in% pQTL$molecular_trait_id & molecular_trait_id %in% sQTL$molecular_trait_id) %>% 
    distinct(molecular_trait_id) %>% 
    pull(molecular_trait_id)

protein_conversion_list <- pQTL %>% 
    mutate(UniProt = str_remove(cs_id,'_.*')) %>% 
    distinct(UniProt,molecular_trait_id)
gene_conversion_list <- eQTL %>% 
    select(molecular_trait_id,gene_name) %>% 
    distinct() %>% 
    mutate(gene = str_remove(molecular_trait_id,'\\..*')) %>% 
    left_join(protein_conversion_list,by = c('molecular_trait_id'))
#merged_qtl_finemapping <- bind_rows(
#                        pQTL  %>% mutate(type = 'pQTL'),
#                        eQTL %>%  mutate(type = 'eQTL'),
#                        sQTL %>%  mutate(type = 'sQTL')
#                        ) %>% 
#                        select(molecular_trait_id,variant,pip,type)  %>% 
#                        pivot_wider(names_from = type,values_from = pip,values_fill = 0) 
#merged_qtl_finemapping %>% head
#
#
eQTL_pQTL_gene_level <- fread('/home/jupyter/multiomics_qtls/edit/FastEnlocInput/eQTL_pQTL.enloc.gene.out') %>% 
    mutate(type = 'eQTL-pQTL') %>% 
    add_fdr_threshold_label(GLCP)

sQTL_pQTL_gene_level <- fread('/home/jupyter/multiomics_qtls/edit/FastEnlocInput/sQTL_pQTL.enloc.gene.out') %>% 
    mutate(type = 'sQTL-pQTL') %>% 
    add_fdr_threshold_label(GLCP)

eQTL_sQTL_gene_level <-  fread('/home/jupyter/multiomics_qtls/edit/FastEnlocInput/eQTL_sQTL.enloc.gene.out') %>% 
    mutate(type = 'eQTL-sQTL') %>% 
    add_fdr_threshold_label(GLCP)
gene_level_colocs <- bind_rows(
                    eQTL_pQTL_gene_level,
                    sQTL_pQTL_gene_level,
                    eQTL_sQTL_gene_level
                    ) %>% 
                    left_join(genebayes,by = c('Gene'  = 'ensg'))
gene_level_colocs %>% head


eQTL_pQTL_CS_level <- LoadEnlocCSLevel('/home/jupyter/multiomics_qtls/edit/FastEnlocInput/eQTL_pQTL.enloc.sig.out') %>% 
    mutate(type = 'eQTL-pQTL') %>% 
    mutate(RCP = as.numeric(RCP),LCP = as.numeric(LCP))  %>% 
    add_fdr_threshold_label(RCP)


sQTL_pQTL_CS_level <- LoadEnlocCSLevel('/home/jupyter/multiomics_qtls/edit/FastEnlocInput/sQTL_pQTL.enloc.sig.out') %>% 
    mutate(type = 'sQTL-pQTL') %>% 
    mutate(RCP = as.numeric(RCP),LCP = as.numeric(LCP))  %>% 
    add_fdr_threshold_label(RCP)


eQTL_sQTL_CS_level <-  LoadEnlocCSLevel('/home/jupyter/multiomics_qtls/edit/FastEnlocInput/eQTL_sQTL.enloc.sig.out') %>% 
    mutate(type = 'eQTL-sQTL') %>% 
    mutate(RCP = as.numeric(RCP),LCP = as.numeric(LCP)) %>% 
    add_fdr_threshold_label(RCP)
cs_level_colocs <- bind_rows(
                    eQTL_sQTL_CS_level,
                    sQTL_pQTL_CS_level,
                    eQTL_pQTL_CS_level
                    ) %>% 
    separate(Signal,into = c('cs_1','cs_2'),sep = '\\(@\\)',remove = FALSE)  %>% 
    mutate(cs_1 = str_remove(cs_1,'^ENSG[0-9]+:'),cs_2 = str_remove(cs_2,'^ENSG[0-9]+:'))


eQTL_sQTL_coloc <- LoadEnlocVariantLevel('/home/jupyter/multiomics_qtls/edit/FastEnlocInput/eQTL_sQTL.enloc.snp.out') %>%
    mutate(type = 'eQTL-sQTL') %>% 
    mutate(across(contains('PIP'), ~as.numeric(.)),
           Gene = str_remove(Signal,':.*')) %>% 
    mutate(SCP = as.numeric(SCP)) %>% 
    add_fdr_threshold_label(SCP)

eQTL_pQTL_coloc <- LoadEnlocVariantLevel('/home/jupyter/multiomics_qtls/edit/FastEnlocInput/eQTL_pQTL.enloc.snp.out') %>%
    mutate(type = 'eQTL-pQTL') %>% 
    mutate(across(contains('PIP'), ~as.numeric(.)),
           Gene = str_remove(Signal,':.*')) %>% 
    mutate(SCP = as.numeric(SCP)) %>% 
    add_fdr_threshold_label(SCP)

       
sQTL_pQTL_coloc <- LoadEnlocVariantLevel('/home/jupyter/multiomics_qtls/edit/FastEnlocInput/sQTL_pQTL.enloc.snp.out') %>%
    mutate(type = 'sQTL-pQTL') %>% 
    mutate(across(contains('PIP'), ~as.numeric(.)),
           Gene = str_remove(Signal,':.*'))  %>% 
    mutate(SCP = as.numeric(SCP)) %>% 
    add_fdr_threshold_label(SCP)


variant_level_coloc <- bind_rows(
                                eQTL_sQTL_coloc,
                                eQTL_pQTL_coloc,
                                sQTL_pQTL_coloc
                                )
annotated_variant_level_coloc <- variant_level_coloc%>% 
    left_join(AnnotationData,by = 'SNP') %>% 
    mutate(rare = gvs_max_af < .01)  
annotated_variant_level_coloc %>% write_tsv('CrossQTLVariantLevelColoc.Annotated.tsv.gz')

#annotated_variant_level_coloc <- fread('CrossQTLVariantLevelColoc.Annotated.tsv.gz')
fixed_afs <- annotated_variant_level_coloc %>% 
    select(contains('gvs') & contains('af')) %>% 
    select(-gvs_all_af,-gvs_max_af,-contains('ac'),-contains('an')) %>% 
    mutate(across(everything(),~case_when(. > .5 ~ 1 - .,TRUE ~ .))) %>% 
    mutate(max_af = pmax(gvs_afr_af,gvs_amr_af,gvs_eas_af,gvs_eur_af,gvs_mid_af,gvs_sas_af,gvs_oth_af))  %>% 
    select(max_af)
temp_pqtl <- pQTL %>% select(variant,molecular_trait_id,posterior_mean,posterior_sd)
temp_eqtl <- eQTL %>% select(variant,molecular_trait_id,posterior_mean,posterior_sd)
temp_sqtl <- sQTL %>% select(variant,molecular_trait_id,posterior_mean,posterior_sd)

variant_level_cleaned <- annotated_variant_level_coloc %>% 
    bind_cols(max_af= fixed_afs$max_af) %>% 
    mutate(rare = max_af < .01) %>% 
    mutate(log_maf = -log10(max_af)) %>%
    mutate(log_maf = scale(log_maf)) %>% 
    left_join(genebayes, by = c('Gene' = 'ensg'))  %>% 
    filter(Coloc == TRUE) 

options(repr.plot.width = 10, repr.plot.height = 6, repr.plot.res = 400)
library(patchwork)
# Sup Fig X - Plot posterior distributions across each colocalizaiton type and feature

GLCP_dist <- gene_level_colocs %>% 
    ggplot(aes(x = GLCP,fill = type)) +
        geom_density(alpha = 0.7,show.legend = FALSE) +
        theme(axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank()
             ) + 
            labs(tag = "A")
RCP_dist <- cs_level_colocs %>% 
    ggplot(aes(x = RCP,fill = type)) +
        geom_density(show.legend = FALSE) + 
        theme(axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank()
             )
SCP_dist <- annotated_variant_level_coloc %>% 
    ggplot(aes(x = SCP,fill = type)) + 
        geom_density() + 
        theme(axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank()
             )
fastenloc_posterior_dists <- GLCP_dist + RCP_dist + SCP_dist
ggsave('CrossQTLPosteriorDistributions.pdf',height = 7,width= 7,plot = fastenloc_posterior_dists)
system(paste0('gsutil cp CrossQTLPosteriorDistributions.pdf gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/PaperFigures/SupplementalFigures/CrossQTL'))
fastenloc_posterior_dists


gene_discovery_across_thresholds <- gene_level_colocs %>% 
    group_by(type) %>% 
    summarise(
        n_tested = n(),
        n_0.01 = sum(GLCP > 0.01, na.rm = TRUE),
        n_0.1  = sum(GLCP > 0.1, na.rm = TRUE),
        n_0.5  = sum(GLCP > 0.5, na.rm = TRUE),
        n_0.9  = sum(GLCP > 0.9, na.rm = TRUE)
      ) %>% 
    mutate(posterior = 'GLCP')

cs_discovery_across_thresholds <- cs_level_colocs %>% 
    group_by(type) %>% 
    summarise(
        n_tested = n(),
        n_0.01 = sum(RCP > 0.01, na.rm = TRUE),
        n_0.1  = sum(RCP > 0.1, na.rm = TRUE),
        n_0.5  = sum(RCP > 0.5, na.rm = TRUE),
        n_0.9  = sum(RCP > 0.9, na.rm = TRUE)
      ) %>% 
    mutate(posterior = 'RCP')

variant_discovery_across_thresholds <- annotated_variant_level_coloc %>% 
    group_by(type) %>% 
    summarise(
        n_tested = n(),
        n_0.01 = sum(SCP > 0.01, na.rm = TRUE),
        n_0.1  = sum(SCP > 0.1, na.rm = TRUE),
        n_0.5  = sum(SCP > 0.5, na.rm = TRUE),
        n_0.9  = sum(SCP > 0.9, na.rm = TRUE)
      ) %>% 
    mutate(posterior = 'SCP')


discovery_by_type <- bind_rows(variant_discovery_across_thresholds,
                               cs_discovery_across_thresholds,
                               gene_discovery_across_thresholds)
CrossQTLDiscoveyrAcrossthresholds <- discovery_by_type %>% 
    pivot_longer(!c(type,n_tested,posterior))  %>% 
    mutate(threshold = as.numeric(str_remove(name,'n_')))  %>% 
    ggplot(aes(x = threshold,y = log10(value),group = type,color = type)) + 
        geom_point() +  
        geom_line() + 
        facet_wrap(~posterior,nrow = 1) +
        ylab('log10(Number Colocs)') +
        labs(tag = "B")

SuppFigPlotCrossQTLQC <-  fastenloc_posterior_dists /CrossQTLDiscoveyrAcrossthresholds 
ggsave('CrossQTLPosteriorDistributions.pdf',height = 7,width= 10,plot = SuppFigPlotCrossQTLQC)
system(paste0('gsutil cp CrossQTLPosteriorDistributions.pdf gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/PaperFigures/SupplementalFigures/CrossQTL'))
SuppFigPlotCrossQTLQC  


gene_level_colocs %>% 
    group_by(Gene) %>% 
    summarize(any_coloc = any(Coloc)) %>% 
    summarize(number_coloc = sum(any_coloc),prop_coloc = mean(any_coloc))

gene_level_colocs %>% 
    group_by(type) %>% 
    summarize(number_coloc = sum(Coloc),prop_coloc = mean(Coloc),total = dplyr::n())

all_proportions <- gene_level_colocs %>% 
    filter(Gene %in% QTL_union_gene_list) %>% 
    group_by(Gene) %>% 
    summarize(NumColocs = sum(Coloc)) %>% 
    dplyr::count(NumColocs) %>% 
    mutate(prop = n/sum(n)) 
all_proportions

gene_level_colocs %>% 
    filter(Coloc == TRUE) %>% 
    filter(str_detect(type,'pQTL')) %>%
    distinct(Gene) %>% 
    nrow


all_proportions %>% 
    filter(NumColocs >0) %>% 
    summarize(sum(prop))

# Sup Fig X - Evaluate proportion of genes with a cross QTL colocalization
options(repr.plot.width = 6, repr.plot.height = 4, repr.plot.res = 400)

ProportionGeneColocs <- gene_level_colocs %>%  
    group_by(type) %>% 
    summarize(PropColoc = mean(Coloc)) %>% 
    ggplot(aes(x = type,y = PropColoc)) + 
        geom_col(position = 'dodge') +
        theme(axis.line.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text= element_text(size = 9),
              axis.title= element_text(size = 11)
             )  +
        ylab('Proportion Colocalized')
ggsave('CrossQTLPropGenesColoc.pdf')
ProportionGeneColocs
system(paste0('gsutil cp CrossQTLPropGenesColoc.pdf gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/PaperFigures/SupplementalFigures/CrossQTL/'))


# Collapse credible sets by colocalization

cleaned_eqtl_sqtl <- cs_level_colocs %>% 
    filter(Coloc == TRUE & type == 'eQTL-sQTL') %>% 
    separate(Signal, into = c('sig_1','sig_2'),sep = '\\(@\\)') %>% 
    mutate(sig_2_list = str_split(sig_2, "_(?=ENSG[0-9]+:chr)")) %>% 
    unnest(sig_2_list) %>% 
    mutate(gene_1 = str_extract(sig_1,"^ENSG[0-9]+"), gene_2 =str_extract(sig_2_list,"^ENSG[0-9]+") ) %>% 
    filter(gene_1 == gene_2) %>% 
    select(cs_1,sig_2_list,RCP,LCP,Coloc,type,gene_1)


cleaned_eqtl_pqtl <- cs_level_colocs %>% 
    filter(Coloc == TRUE & type == 'eQTL-pQTL')  %>% 
    separate(Signal, into = c('sig_1','sig_2'),sep = '\\(@\\)') %>% 
    mutate(sig_2_list = str_split(sig_2, "_(?=ENSG[0-9]+:)")) %>% 
    unnest(sig_2_list) %>% 
    mutate(gene_1 = str_extract(sig_1,"^EN[<43;65;42MSG[0-9]+"), gene_2 =str_extract(sig_2_list,"^ENSG[0-9]+") ) %>% 
    filter(gene_1 == gene_2) %>% 
    select(cs_1,sig_2_list,RCP,LCP,Coloc,type,gene_1)


cleaned_sqtl_pqtl <-  cs_level_colocs %>% 
    filter(Coloc == TRUE & type == 'sQTL-pQTL') %>% 
    separate(Signal, into = c('sig_1','sig_2'),sep = '\\(@\\)') %>% 
    mutate(sig_2_list = str_split(sig_2, "_(?=ENSG[0-9]+:)")) %>% 
    unnest(sig_2_list) %>% 
    mutate(gene_1 = str_extract(sig_1,"^ENSG[0-9]+"), gene_2 =str_extract(sig_2_list,"^ENSG[0-9]+") ) %>% 
    filter(gene_1 == gene_2) %>% 
    select(cs_1,sig_2_list,RCP,LCP,Coloc,type,gene_1)


cleaned_cs_level_colocs <- bind_rows(
                                cleaned_eqtl_pqtl,
                                cleaned_eqtl_pqtl,
                                cleaned_sqtl_pqtl
                                ) %>% 
                                mutate(cs_2 = str_remove(sig_2_list,'^ENSG.*:')) %>% 
                                select(gene_1,RCP,LCP,type,cs_1,cs_2)  %>%
                                dplyr::rename('gene' = 'gene_1')




edges <- cleaned_cs_level_colocs %>%
  transmute(
    from = cs_1,
    to = cs_2,
    type = type
  ) %>%
  filter(!is.na(from), !is.na(to)) %>%
  distinct(from, to)

g <- graph_from_data_frame(edges, directed = FALSE)

components_df <- tibble(
  cs_id = names(components(g)$membership),
  coloc_component = as.integer(components(g)$membership)
)

cs_level <- cs_list %>%
  rename(Gene = molecular_trait_id, cs_id = cs_id) %>%
  left_join(components_df, by = "cs_id") %>%
  mutate(
    has_coloc = !is.na(coloc_component)
  )

coloc_summary <- cs_level %>% 
    group_by(Gene) %>% 
    mutate(component = case_when(is.na(coloc_component) ~ row_number(),TRUE ~ coloc_component)) %>% 
    select(Gene,has_coloc,component) %>% 
    distinct() %>% 
    summarize(PropColoc = mean(has_coloc),NumColoc = sum(has_coloc),Total = dplyr::n()) %>% 
    left_join(median_tpms,by = c('Gene')) %>% 
    left_join(genebayes,by = c('Gene' = 'ensg')) %>% 
    filter(is.finite(mean_TPM) & is.finite(post_mean))  %>% 
    left_join(gene_conversion_list,by = c('Gene' = 'gene')) %>% 
    left_join(LOD_data,by = c('UniProt' ))

coloc_summary

proportion_colocs <- all_proportions %>% 
    mutate(NumColocs = factor(NumColocs, levels = c(0, 1, 2, 3))) %>%
    ggplot(aes(x = NumColocs, y = prop, fill = NumColocs)) + 
        geom_col(width = 0.75) +
        scale_fill_manual(
            values = c(
                "0" = "#D9D9D9",
                "1" = "#A9C1D9",
                "2" = "#7B9AB8",
                "3" = "#577590"
            ),
            guide = "none"
        ) +
        annotate(
            "segment",
            x = 1.5,
            xend = 4.5,
            y = 0.455,
            yend = 0.455,
            color = "grey50",
            linewidth = 0.25
        ) +
        annotate(
            "text",
            x = 3,
            y = 0.465,
            label = paste0("86.2%", " shared across ≥1 molecular layer"),
            size = 3,
            color = "grey20"
        ) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),breaks = seq(0, 0.5, by = 0.1),limits = c(0, 0.55)) +
        ylab("Percent of Genes") +
        xlab("Number of Colocalized Molecular Layers") +
        theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text = element_text(size = 9),
            axis.title = element_text(size = 11)
        )
cs_reds <- c(
  "0" = "#E9DCDC",
  "1" = "#D8BBBB",
  "2" = "#C89A9A",
  "3" = "#B37979",
  "4" = "#9E5E5E",
  "5" = "#8A4D4D",
  "6" = "#763F3F",
  "7" = "#643232"
)

cs_coloc_plot <- cs_level %>% 
  na.omit() %>% 
  distinct(Gene, coloc_component) %>% 
  count(Gene) %>%
  count(n, name = "nn") %>% 
  ggplot(aes(x = factor(n), y = nn / sum(nn), fill = factor(n))) + 
      geom_col(width = 0.75) + 
      scale_fill_manual(values = cs_reds, guide = "none") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1),    limits = c(0, 0.55)) +
      xlab("Number of colocalized signals") + 
      ylab("Percent of genes") +
      theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 11)
      )

cs_coloc_rate <- cs_level %>% 
    filter(Gene %in% QTL_union_gene_list) %>% 
  group_by(type) %>%  
  summarise(
    n_cs = n_distinct(cs_id),
    n_coloc = n_distinct(cs_id[has_coloc]),
    coloc_rate = n_coloc / n_cs
  )

cs_coloc_rate_plot <- cs_coloc_rate %>% 
    ggplot(aes(y = reorder(type,coloc_rate),x = coloc_rate,fill = type)) + 
        geom_col(show.legend = FALSE) + 
        theme(axis.line.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.y = element_line(color = 'lightgrey'),
              axis.ticks.y = element_blank()
             ) +
        xlab('% CS coloc') +
        ylab("") +
        paletteer::scale_fill_paletteer_d("ltc::minou")   +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1),expand = c(0,0 )) + 
        labs(tag = NULL)
p1 <- proportion_colocs +
  ggtitle("Gene-level Colocalization") +
  labs(tag = "A") +
  theme(
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    plot.title.position = "plot",
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0, 1)
  )

p2 <- cs_coloc_plot +
  ggtitle("Signal-level Colocalization") +
  labs(tag = "B") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0, 1)
  )

inset_plot <- cs_coloc_rate_plot +
  labs(tag = NULL) +
  theme(
    plot.tag = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

p2_inset <- p2 +
  inset_element(
    inset_plot,
    left   = 0.48,
    bottom = 0.60,
    right  = 0.98,
    top    = 0.98,
    align_to = "panel"
  )

coloc_count_summary <-
  p1 + p2_inset +
  plot_layout(widths = c(1, 1))

coloc_count_summary

eQTL_sQTL_model_input <- variant_level_cleaned %>% 
    filter(type == 'eQTL-sQTL') %>% 
    left_join(temp_eqtl,by = c('SNP' = 'variant','Gene' = 'molecular_trait_id')) %>% 
    left_join(temp_sqtl,by = c('SNP' = 'variant','Gene' = 'molecular_trait_id')) 

eQTL_pQTL_model_input <- variant_level_cleaned %>% 
    filter(type == 'eQTL-pQTL') %>% 
    left_join(temp_eqtl,by = c('SNP' = 'variant','Gene' = 'molecular_trait_id')) %>% 
    left_join(temp_pqtl,by = c('SNP' = 'variant','Gene' = 'molecular_trait_id')) 

sQTL_pQTL_model_input <- variant_level_cleaned %>% 
    filter(type == 'sQTL-pQTL') %>% 
    left_join(temp_sqtl,by = c('SNP' = 'variant','Gene' = 'molecular_trait_id')) %>% 
    left_join(temp_pqtl,by = c('SNP' = 'variant','Gene' = 'molecular_trait_id')) 


eQTL_sQTL_model_res <- eQTL_sQTL_model_input %>% 
    do(broom::tidy(lm(abs(posterior_mean.y) ~ abs(posterior_mean.x)*log_maf,data = .),conf.int = TRUE)) %>% 
    filter(!str_detect(term,'Int')) %>% 
    mutate(type = 'eQTL-sQTL')


eQTL_pQTL_model_res <- eQTL_pQTL_model_input %>% 
    do(broom::tidy(lm(log1p(posterior_mean.y) ~ (log1p(posterior_mean.x))*log_maf,data = .),conf.int = TRUE)) %>% 
    filter(!str_detect(term,'Int')) %>% 
    mutate(type = 'eQTL-pQTL')


sQTL_pQTL_model_res <- sQTL_pQTL_model_input %>% 
    do(broom::tidy(lm(abs(posterior_mean.y) ~ abs(posterior_mean.x)*log_maf,data = .),conf.int = TRUE)) %>% 
    filter(!str_detect(term,'Int')) %>% 
    mutate(type = 'sQTL-pQTL')



allele_frequncy_interaction_model_res <- bind_rows(
                                            eQTL_sQTL_model_res,
                                            eQTL_pQTL_model_res,
                                            sQTL_pQTL_model_res
                                            )  
                

allele_frequncy_interaction_model_res 
#model <- lm(z_x ~ z_y + log_maf  ,data = model_input)

#binned_correlation <- model_input %>%
#    filter(!is.na(max_af)) %>%   
#    mutate(
#        maf_decile = ntile(log_maf, 10),
#        maf_bin = case_when(max_af >= .05 ~ 'Common',max_af < .05 & max_af > .01 ~ 'low-freq',TRUE ~'rare' )
#      ) 
rare_gene_list <- annotated_variant_level_coloc %>% filter(Coloc == TRUE & type == 'sQTL-pQTL') %>% distinct(Gene)


 annotated_variant_level_coloc  %>% 
    filter( type == 'eQTL-pQTL' ) %>% 
    #filter(Gene %in% rare_gene_list$Gene)   %>% 
    mutate(gvs_max_af = as.numeric(-log10(gvs_max_af))) %>% 
    #select(Gene,gvs_max_af,Coloc) %>% 
    do(broom::tidy(glm(Coloc ~ gvs_max_af + PIP_qtl + PIP_gwas_marginal,family ='binomial',data = .)))

eQTL_sQTL_model_res <- eQTL_sQTL_model_input %>% 
    filter(!is.na(posterior_mean.y) & !is.na(posterior_mean.x)) %>% 
    do(broom::tidy(lm(log1p(abs(posterior_mean.y)) ~ log1p(abs(posterior_mean.x))*log_maf,data = .),conf.int = TRUE)) %>% 
    filter(!str_detect(term,'Int')) %>% 
    mutate(type = 'eQTL-sQTL')


eQTL_pQTL_model_res <- eQTL_pQTL_model_input %>% 
    filter(!is.na(posterior_mean.y) & !is.na(posterior_mean.x)) %>% 
    do(broom::tidy(lm(log1p(posterior_mean.y) ~ (log1p(posterior_mean.x))*log_maf,data = .),conf.int = TRUE)) %>% 
    filter(!str_detect(term,'Int')) %>% 
    mutate(type = 'eQTL-pQTL')


sQTL_pQTL_model_res <- sQTL_pQTL_model_input %>% 
    filter(!is.na(posterior_mean.y) & !is.na(posterior_mean.x)) %>% 
    do(broom::tidy(lm(log1p(abs(posterior_mean.y)) ~ log1p(abs(posterior_mean.x))*log_maf,data = .),conf.int = TRUE)) %>% 
    filter(!str_detect(term,'Int')) %>% 
    mutate(type = 'sQTL-pQTL')



allele_frequncy_interaction_model_res <- bind_rows(
                                            eQTL_sQTL_model_res,
                                            eQTL_pQTL_model_res,
                                            sQTL_pQTL_model_res
                                            )  
allele_frequncy_interaction_model_res

MAF_cross_ome_correlation<- allele_frequncy_interaction_model_res %>% 
    filter(str_detect(term,':')) %>% 
    ggplot(aes(x = estimate,y = type)) + 
        geom_linerange(aes(xmin = conf.low,xmax = conf.high),color ='grey40') + 
        geom_point(size = 4,color = '#364B63') + 
        geom_vline(xintercept = 0,linetype = 'dashed',color = 'grey40')  + 
        theme(axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              plot.margin = margin(l = 60, unit = "pt")
             ) + 
        xlab('Effect of -log(MAF) on cross-ome correlation') + 
        labs(tag = "C") 
MAF_cross_ome_correlation

options(repr.plot.width = 18, repr.plot.height = 6, repr.plot.res = 400)

top_row <- (coloc_count_summary | MAF_cross_ome_correlation) + 
        plot_layout(widths = c(2.5,1)) 
ggsave('Fig3TopHalf.pdf',height = 6,width = 12)
system(paste0('gsutil cp Fig3TopHalf.pdf  gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/PaperFigures/'))

options(repr.plot.width = 6, repr.plot.height = 3, repr.plot.res = 400)

binned_data <- eQTL_pQTL_model_input %>%
  mutate(maf_bin = case_when(max_af >= .05 ~ 'Common',max_af < .05 & max_af >= .01 ~ 'Low-Freq',TRUE ~'Rare')) %>%
  filter(!is.na(maf_bin))  %>% 
  mutate(posterior_mean.x = log1p(posterior_mean.x),posterior_mean.y = log1p(posterior_mean.y))

lims <- binned_data %>%
  summarise(
    lim = quantile(abs(c(posterior_mean.x, posterior_mean.y)), 0.97, na.rm = TRUE)
  ) %>%
  pull(lim)

slopes <- binned_data %>%
  group_by(maf_bin) %>%
  summarise(
    n = n(),
    slope = coef(lm(posterior_mean.y ~ posterior_mean.x ))[["posterior_mean.x"]],
    .groups = "drop"
  ) %>%
  mutate(label = paste0("β = ", round(slope, 2), "\nN = ", scales::comma(n)))

effetct_size_comparison <- ggplot(binned_data , aes(x = posterior_mean.x,y =  posterior_mean.y)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.35, color = "grey65") +

  geom_point(
    data = . %>% filter(Coloc == TRUE),
    alpha = 0.6,
    size = 4,
    color = '#6EC5C3'
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.7,
    color = "red",
    fill = "red",
    alpha = 0.12
  ) +
  #scale_color_manual(values = c('grey40','#42B7B9FF')) +
  geom_label(
    data = slopes,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.15,
    size = 2.7,
    label.size = 0,
    fill = "white",
    alpha = 0.85
  ) +
  facet_wrap(~maf_bin, nrow = 1) +
  coord_cartesian(xlim = c(-lims, lims), ylim = c(-lims, lims)) +
  labs(
    x = "Expression Posterior",
    y = "Protein Posterior"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    panel.spacing.x = unit(0.8, "lines"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )  
  
effetct_size_comparison

options(repr.plot.width = 18, repr.plot.height = 12, repr.plot.res = 400)
blank_plot <- ggplot() + theme_void()
top_row/(effetct_size_comparison + blank_plot + plot_layout(widths = c(1.3,1))) + 
        plot_annotation(tag_levels = "A")






