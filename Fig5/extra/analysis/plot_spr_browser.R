source("_setup.R")
VPOS <- 72887998; TSS <- 72887382
d <- read_tsv(file.path(AG,"spr_tracks_window.tsv.gz"), show_col_types=FALSE) %>%
  mutate(mark = case_when(output_type=="CHIP_HISTONE"~histone_mark,
                          output_type=="CHIP_TF"~tf, TRUE~output_type)) %>%
  select(allele,pos,value,biosample,strand,mark)

pick <- tribble(
  ~layer,        ~context, ~biosample,                         ~mark,      ~strand,
  "DNase",       "Brain",  "brain",                            "DNASE",    ".",
  "DNase",       "Blood",  "GM12878",                          "DNASE",    ".",
  "CAGE (init.)","Brain",  "brain",                            "CAGE",     "+",
  "CAGE (init.)","Blood",  "GM12878",                          "CAGE",     "+",
  "H3K27ac",     "Brain",  "dorsolateral prefrontal cortex",   "H3K27ac",  ".",
  "H3K27ac",     "Blood",  "GM12878",                          "H3K27ac",  ".",
  "H3K27me3",    "Brain",  "brain",                            "H3K27me3", ".",
  "H3K27me3",    "Blood",  "peripheral blood mononuclear cell","H3K27me3", ".",
  "CTCF",        "Brain",  "brain",                            "CTCF",     ".",
  "CTCF",        "Blood",  "GM12878",                          "CTCF",     ".") %>%
  mutate(rid=row_number())

dat <- pick %>% left_join(d, by=c("biosample","mark","strand"),
                          relationship="many-to-many")

eff <- dat %>% filter(pos>VPOS-250, pos<VPOS+250) %>%
  group_by(rid,allele) %>% summarise(s=sum(value),.groups="drop") %>%
  pivot_wider(names_from=allele, values_from=s) %>% mutate(l2=log2((alt+1)/(ref+1)))

pick2 <- pick %>% left_join(eff,by="rid") %>%
  mutate(bs=recode(biosample, "dorsolateral prefrontal cortex"="cortex BA9",
                   "peripheral blood mononuclear cell"="PBMC"),
         lab=sprintf("%s — %s  (log2 %+.1f)", layer, bs, l2),
         lab=fct_reorder(lab, rid))

plotdat <- dat %>% left_join(pick2 %>% select(rid,lab), by="rid") %>%
  mutate(kb=pos/1e3, lab=factor(lab, levels=pick2 %>% arrange(rid) %>% pull(lab)))

p <- ggplot(plotdat, aes(kb, value)) +
  annotate("rect", xmin=(VPOS-2)/1e3, xmax=(VPOS+3)/1e3, ymin=-Inf, ymax=Inf, fill="#F4A261", alpha=0.35) +
  geom_vline(xintercept=TSS/1e3, linetype="dashed", linewidth=0.3, colour="grey45") +
  geom_area(data=~filter(.x, allele=="ref"), fill="grey55", colour=NA, alpha=0.55) +
  geom_line(data=~filter(.x, allele=="alt"), colour="#C0392B", linewidth=0.45) +
  facet_wrap(~lab, ncol=2, scales="free_y", dir="h") +
  scale_x_continuous(breaks=scales::pretty_breaks(4)) +
  labs(x="chr2 position (kb)   —   dashed = SPR TSS, orange band = deletion (GCGCC>G)",
       y="predicted signal   (ref = grey fill, alt = red line)") +
  theme_classic(base_size=8, base_family="inter") +
  theme(strip.background=element_blank(), strip.text=element_text(size=7.2,face="bold",hjust=0),
        axis.text=element_text(size=6.6,colour="grey20"), axis.title=element_text(size=7.4),
        panel.spacing=unit(3,"mm"), plot.margin=margin(4,6,4,6))
ggsave(file.path(FIGDIR,"FigureSPR_browser.pdf"), p, width=180, height=190, units="mm", device=cairo_pdf)
ggsave(file.path(FIGDIR,"FigureSPR_browser.png"), p, width=180, height=190, units="mm", dpi=300)
cat("wrote FigureSPR_browser.pdf/png\n"); print(pick2 %>% select(layer,context,bs,ref,alt,l2))
