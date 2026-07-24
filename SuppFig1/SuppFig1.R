library(tidyverse)
library(data.table)
library(rcartocolor)
library(ggrepel)
library(plyranges)
library(patchwork)
library(enrichR)


theme_set(theme_classic())
theme_update(panel.border = element_rect(fill = NA,linewidth = .9), 
            axis.line = element_blank(),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 20))


#gencode_V48_GTF <- '/home/jupyter/gencode.v48.annotation.gtf.gz'
#gencode_v48_tpms <- 'gs://prod-drc-broad/aou_rnaseq/final_10k_rna_20250721/rnaseqc2/aou_rnaseq_20250717.gene_tpm.gct.gz'
gencode_V45_tpms <- 'gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/QTL_5500_subset/expression/aou_rnaseq_121624.gene_tpm.fixed.gct.gz'
data_dir <- 'gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/RNA_seq/aou_rnaseq/1000291_gencode_v48/'
#gene_results <- 'gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/RNA_seq/aou_rnaseq/1000291_gencode_v48/1000291.gene_tpm.gct.gz'
#isoform_results <- 'gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/RNA_seq/aou_rnaseq/1000291_gencode_v48/1000291.rsem.isoforms.results.gz'


# pull TF list data 
TF_list_link <- 'https://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt'
system(paste0('wget ',TF_list_link))
TF_list <- fread(basename(TF_list_link),header = FALSE)


compute_mean_TPM <- function(TPM_df,gene_annotations){
annotation_dat <- gene_annotations %>% distinct()
TPM_means <- TPM_df %>% 
    pivot_longer(!c(Name,Description)) %>% 
    group_by(Description,Name) %>% 
    summarize(mean_TPM = mean(value),SE = sd(value)/dplyr::n())  %>% 
    left_join(annotation_dat,by = c('Name' = 'gene_id'))   
TPM_means
}
plot_gene_comparison <- function(TPM_df,point_size = 3,detection_threshold = 1, gene_type = '',hex = FALSE){

axis_label <- paste0('Unique ',gene_type, ' genes')
cor_values <- cor(log2(TPM_df$mean_TPM_AoU+1),log2(TPM_df$mean_TPM_GTEx+1),use = 'complete.obs')
plot_text <- paste0('Cor = ',round(cor_values,digits = 2))

cor_plot <- TPM_df %>% 
    ggplot(aes(x = log2(mean_TPM_AoU + 1),y = log2(mean_TPM_GTEx +1))) + 
        xlab('log2(AoU mean TPM)') + 
        ylab('log2(GTEx mean TPM)') +
        annotate("text",x = 1.8, y =19,label = plot_text,size = 7) 

if (hex == FALSE){
cor_plot <- cor_plot + 
    geom_point(shape = 21,size =3,fill ='black',alpha = 0.7) +  
    geom_abline(slope = 1 ,color = 'red') 
}
else{
cor_plot <- cor_plot + 
        geom_hex(aes(fill = stat(log10(count))))    +
        scale_fill_viridis_c() + 
        geom_abline(slope = 1 ,color = 'red')  

}  
 
detection_plot <- TPM_df %>% 
    filter((mean_TPM_AoU > detection_threshold & mean_TPM_GTEx < detection_threshold)|(mean_TPM_AoU < detection_threshold & mean_TPM_GTEx > detection_threshold)  ) %>% 
    mutate(type = case_when(mean_TPM_AoU > mean_TPM_GTEx ~ 'AoU' ,TRUE ~ 'GTEx')) %>% 
    dplyr::count(type)  %>% 
    ggplot(aes(x = type,y = n,fill = type)) + 
        geom_col(show.legend = FALSE) +
        ylab(axis_label) +
        xlab('Cohort') + 
        scale_fill_carto_d(palette = 'Safe')
    
outplot <- cor_plot + detection_plot  
outplot
}
gencode_V45_GTF <- 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz'
system(paste0('gsutil cp ',gencode_V45_GTF,' .'))

gtf_df <- rtracklayer::import(gencode_V45_GTF)
gene_annotations <- gtf_df %>% data.frame() %>% select(gene_id,gene_type) %>% distinct()


gencode_V45_tpms <- 'gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/QTL_5500_subset/expression/aou_rnaseq_121624.gene_tpm.fixed.gct.gz'
GTEx_whole_blood_path <- 'https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_whole_blood.gct.gz'

#system(paste0('gsutil cp ',gencode_V45_tpms, ' .'))
#GTEx_whole_blood <- fread(basename(GTEx_whole_blood_path),skip =2 )
#AoU_gencode_V45_TPMs_df <- fread(basename(gencode_V45_tpms),skip = 2 )


# load AoU GENCODE 45 TPMs and compute average TPM per gene
AoU_mean_TPMs_v45  <- AoU_gencode_V45_TPMs_df %>% 
            compute_mean_TPM(gene_annotations) %>% 
            distinct() %>% 
            ungroup() %>% 
            mutate(gene_id = str_remove(Name,'\\..*'))

# load GTEx GENCODE 45 TPMs and compute average TPM per gene
GTEx_mean_TPMs_v45 <- GTEx_whole_blood %>% 
            compute_mean_TPM(gene_annotations) %>% 
            distinct() %>% 
            ungroup() %>% 
            mutate(gene_id = str_remove(Name,'\\..*'))

merged_TPMs_averages <- AoU_mean_TPMs_v45 %>% 
    select(-Name) %>% 
    inner_join(GTEx_mean_TPMs_v45 %>% select(-Description,-gene_type,-Name), by = c('gene_id'),suffix = c('_AoU','_GTEx') ) %>% 
    mutate(mean_TPM_GTEx = case_when(is.na(mean_TPM_GTEx) ~ 0,TRUE ~ mean_TPM_GTEx),
           mean_TPM_AoU = case_when(is.na(mean_TPM_AoU) ~ 0,TRUE ~ mean_TPM_AoU) )


dbs <- c("GO_Molecular_Function_2025", "GO_Cellular_Component_2025", 
     "GO_Biological_Process_2025")

AoU_specific_genes  <- merged_TPMs_averages %>% 
    filter(gene_type == 'protein_coding') %>% 
    filter(mean_TPM_AoU >1  &  mean_TPM_GTEx < 1) %>% 
    pull(Description)

GTEx_specific_genes  <- merged_TPMs_averages %>% 
    filter(gene_type == 'protein_coding') %>% 
    filter(mean_TPM_AoU < 1  &  mean_TPM_GTEx >1) %>% 
    pull(Description)

message('Number of AoU Specific genes:')
AoU_specific_genes %>% length()

message('Number of GTEx specific genes:')
GTEx_specific_genes %>% length()

AoU_GO_enr <- enrichr(AoU_specific_genes,dbs) %>% bind_rows()
GTEx_GO_enr <- enrichr(GTEx_specific_genes,dbs) %>% bind_rows()

options(repr.plot.width = 25, repr.plot.height = 15, repr.plot.res = 200)

AoU_GTEx_enrichment_plot <- bind_rows(
    AoU_GO_enr %>% arrange(`Adjusted.P.value`) %>% filter(row_number() < 10) %>% mutate(type = 'AoU') ,
    GTEx_GO_enr %>% arrange(`Adjusted.P.value`) %>% filter(row_number() < 10)%>% mutate(type = 'GTEx')
         ) %>% 
    mutate(Term = str_remove(Term,'\\(GO.*')) %>% 
    mutate(neglog10_p = -log10(`Adjusted.P.value`)) %>% 
    ggplot(aes(x = -log10(`Adjusted.P.value`),y = reorder(Term,neglog10_p))) + 
            geom_col() + 
            facet_wrap(~type,scales = 'free_y') +
            scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
            xlab('-log10(FDR)') +
            ylab('Term') + 
            geom_vline(xintercept = 1.3) +
            theme(panel.spacing = unit(1, "cm"),strip.text = element_text(size = 20)) +
            ylab("")

AoU_GTEx_enrichment_plot

options(repr.plot.width = 12, repr.plot.height = 7, repr.plot.res = 200)

protein_coding_genes <- merged_TPMs_averages %>% 
    filter(gene_type == 'protein_coding') %>% 
    filter(mean_TPM_AoU > 1  |  mean_TPM_GTEx > 1)



protein_coding_comparison <- protein_coding_genes %>% plot_gene_comparison(gene_type = 'coding',hex = TRUE)
protein_coding_comparison

options(repr.plot.width = 7, repr.plot.height = 5, repr.plot.res = 200)

TF_genes <- merged_TPMs_averages %>% 
    filter(gene_type == 'protein_coding') %>% 
    filter(Description %in% TF_list$V1)  %>% 
    filter(mean_TPM_GTEx > 1 | mean_TPM_AoU > 1) %>% 
    mutate(across(contains('TPM'), ~log2(.+1 ))) 



TF_genes %>% 
    mutate(across(contains('TPM'), ~log2(.+1 ))) %>% 
    do(broom::tidy(lm(mean_TPM_GTEx~mean_TPM_AoU,data =. )))
TF_gene_comparison <- TF_genes %>% 
    select(Description,contains('TPM')) %>% 
    pivot_longer(!Description) %>% 
    mutate(name = str_remove(name,'mean_TPM_')) %>% 
    ggplot(aes(x = name , y= log2(value + 1 ),fill = name)) + 
        geom_boxplot(show.legend = FALSE) +
        ylab('TF Expression \n log2(TPM + 1)') + 
        scale_fill_carto_d(palette = 'Safe')
TF_gene_comparison
dat <- TF_genes %>% 
    mutate(across(contains('TPM'), ~log2(. +1)))

t.test(dat$mean_TPM_AoU,dat$mean_TPM_GTEx,paired = TRUE)

options(repr.plot.width = 12, repr.plot.height = 7, repr.plot.res = 200)

lnc_RNA_genes <- merged_TPMs_averages %>% 
    filter(gene_type == 'lncRNA') 

lncRNA_comparison <- lnc_RNA_genes %>% 
    filter(mean_TPM_AoU > .1 |  mean_TPM_GTEx > .1) %>% 
    plot_gene_comparison(gene_type = 'lncRNA',hex = TRUE)
dat <- lnc_RNA_genes %>% 
    mutate(across(contains('TPM'), ~log2(. +1)))

t.test(dat$mean_TPM_AoU,dat$mean_TPM_GTEx,paired = TRUE)
lncRNA_comparison

options(repr.plot.width = 35, repr.plot.height = 25, repr.plot.res = 300)

top_row <- protein_coding_comparison | (AoU_GTEx_enrichment_plot + theme(plot.margin = margin(l = '2',unit = 'cm')))
bottom_row <- TF_gene_comparison | lncRNA_comparison    & 
                                    theme(plot.margin = margin(l = '2',t = '5',unit = 'cm'))


SuppFig1 <- top_row/bottom_row +   plot_annotation(tag_levels = "A") & theme(
    plot.tag = element_text(size = 25, face = "bold")
  )
SuppFig1
ggsave('SupplementalFig1.pdf',height = 23,width = 35)

snRNA_genes <- merged_TPMs_averages %>% 
    filter(gene_type == 'snRNA') %>% 
    filter(mean_TPM_AoU > .1 |  mean_TPM_GTEx > .1) 


snRNA_genes %>% plot_gene_comparison(gene_type = 'snRNA',hex = FALSE)

miRNA_genes <- merged_TPMs_averages %>% 
    filter(gene_type == 'miRNA') %>% 
    filter(mean_TPM_AoU > .1 |  mean_TPM_GTEx > .1)


miRNA_genes %>% plot_gene_comparison(gene_type = 'miRNA',hex = FALSE)

options(repr.plot.width = 7, repr.plot.height = 7, repr.plot.res = 200)

rRNA_genes <- merged_TPMs_averages %>% 
    filter(gene_type == 'rRNA')
rRNA_genes %>%  select(Description,mean_TPM_AoU,mean_TPM_GTEx) %>% 
    pivot_longer(!Description) %>% 
        mutate(name = str_remove(name,'mean_TPM_')) %>% 
        ggplot(aes(x = value,fill = name)) + 
            geom_histogram(position = 'identity',alpha= 0.7) +
            xlab('Mean rRNA TPM') + 
            ylab('Number genes') + 
            scale_fill_carto_d(palette = 'Safe')

merged_TPMs_averages %>% 
    filter(str_detect(Description,'^HB') & gene_type == 'protein_coding') %>% 
    ggplot(aes(x = log2(mean_TPM_AoU),y = log2(mean_TPM_GTEx))) + 
        geom_point(size = 6 ,shape = 21,fill = 'black',alpha = 0.6) + 
        geom_abline(slope = 1 ,color = 'red')  +
        ylim(-5,20) + 
        xlim(-5,20) + 
        xlab('log2(Mean AoU TPM)') + 
        ylab('log2(Mean GTEx TPM)') + 
        geom_text_repel(aes(label = Description))



options(repr.plot.width = 15, repr.plot.height = 7, repr.plot.res = 200)

gene_detection_summary <- AoU_mean_TPMs_v45 %>% 
    filter(mean_TPM > .1) %>% 
    dplyr::count(gene_type) %>% 
    mutate(gene_type = case_when(gene_type == 'protein_coding' ~ 'Coding',TRUE~gene_type)) %>% 
    filter(!str_detect(gene_type,'IG|TR|pseudo|artifact|Mt|TEC|rRNA'))

small_rnas_list <- c('misc_RNA','sRNA','vault_RNA','scaRNA','rRNA','ribozyme')

other_small_rnas <- gene_detection_summary %>% 
    filter(gene_type %in% small_rnas_list) %>% 
    summarize(gene_type = 'other small RNAs',n = sum(n))


gene_detection_summary %>%
    bind_rows(other_small_rnas) %>% 
    filter(!gene_type %in% small_rnas_list) %>% 
    ggplot(aes(x = gene_type,y = n,fill = gene_type,color = gene_type)) +
        geom_col(show.legend = FALSE,alpha = 0.7) +
        ylab('Number detected \n (TPM > 0.1)') +
        scale_y_log10() +
        scale_fill_carto_d(palette = 'Safe') +
        scale_color_carto_d(palette = 'Safe') +
        geom_text(aes(label = n),vjust = -.25,color = 'black') +
        xlab('Type of gene')


options(repr.plot.width = 7, repr.plot.height = 5, repr.plot.res = 200)

isoform_df <- fread(basename(isoform_results))
isoform_df %>% 
    filter(TPM > 0.1) %>% 
    group_by(gene_id) %>% 
    summarize(n_isoforms = dplyr::n()) %>% 
    ggplot(aes(x = n_isoforms)) + 
        geom_histogram() +
        xlab('Isoforms detected per gene \n (TPM > 0.1)') +
        ylab('Number genes') 
system(paste0('gsutil cp gs://fc-secure-b8771cfd-5455-4292-a720-8533eb501a93/annotations/GO_annotations.csv .'))

GO_data <- fread('GO_annotations.csv')

GO_data %>% head












