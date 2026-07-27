source("_setup.R")
suppressPackageStartupMessages({ library(patchwork) })

## ---- VANGL2: CTCF boundary collapse across neural cells ----
V<-160428770
d<-read_parquet(file.path(AG,"tracks_VANGL2.parquet")) %>% filter(output_type=="CHIP_TF", tf=="CTCF")
labs<-tibble(biosample=c("choroid plexus epithelial cell","brain","glutamatergic neuron",
                         "brain microvascular endothelial cell","astrocyte of the spinal cord","neural cell"),
             nice=c("choroid plexus epith.","brain","glutamatergic neuron",
                    "brain microvasc. endoth.","spinal-cord astrocyte","neural cell"))
d<-d %>% inner_join(labs,by="biosample") %>% mutate(kb=pos/1e3, nice=factor(nice,levels=labs$nice))
rat<-d %>% group_by(nice) %>% summarise(pct=round(100*max(value[allele=="alt"])/max(value[allele=="ref"])),.groups="drop")
d<-d %>% left_join(rat,by="nice") %>% mutate(lab=sprintf("%s  (CTCF alt %d%%)",nice,pct),
        lab=factor(lab,levels=unique(sprintf("%s  (CTCF alt %d%%)",labs$nice,rat$pct[match(labs$nice,rat$nice)]))))
pV<-ggplot(d,aes(kb,value))+
  annotate("rect",xmin=(V-2)/1e3,xmax=(V+2)/1e3,ymin=-Inf,ymax=Inf,fill="#F4A261",alpha=0.55)+
  geom_area(data=~filter(.x,allele=="ref"),fill="grey60",colour="grey55",linewidth=0.2,alpha=0.7)+
  geom_area(data=~filter(.x,allele=="alt"),fill="#C0392B",colour="#C0392B",linewidth=0.25,alpha=0.5)+
  facet_wrap(~lab,ncol=3,scales="free_y")+
  scale_x_continuous(breaks=scales::pretty_breaks(3))+
  labs(x="chr1 position (kb) — orange = variant · grey = ref CTCF, red = alt   (VANGL2 TSS 28 kb away)",
       y="predicted CTCF ChIP",title="VANGL2 chr1:160,428,770 G>C — CTCF/cohesin boundary collapse (aFC +6.05)")+
  theme_classic(base_size=8,base_family="inter")+
  theme(strip.background=element_blank(),strip.text=element_text(size=6.9,face="bold",hjust=0),
        plot.title=element_text(size=8.5,face="bold"),axis.text=element_text(size=6.3,colour="grey20"),
        axis.title=element_text(size=7.2),panel.spacing=unit(2.6,"mm"),plot.margin=margin(4,6,4,6))
ggsave(file.path(FIGDIR,"FigureVANGL2_CTCF.pdf"),pV,width=180,height=105,units="mm",device=cairo_pdf)
ggsave(file.path(FIGDIR,"FigureVANGL2_CTCF.png"),pV,width=180,height=105,units="mm",dpi=300)

## ---- CFH: kidney/liver promoter loss (CAGE + DNase) ----
C<-196651994
cf<-read_parquet(file.path(AG,"tracks_CFH.parquet"))
sel<-tribble(~biosample,~assay,~out,~nice,
 "renal cortical epithelial cell","CAGE","CAGE","CAGE — renal cortical epith.",
 "epithelial cell of proximal tubule","CAGE","CAGE","CAGE — proximal tubule",
 "HepG2","CAGE","CAGE","CAGE — HepG2 (liver)",
 "kidney epithelial cell","DNASE","DNASE","DNase — kidney epithelial",
 "epithelial cell of proximal tubule","DNASE","DNASE","DNase — proximal tubule",
 "HepG2","DNASE","DNASE","DNase — HepG2 (liver)")
cd<-cf %>% inner_join(sel,by=c("biosample","output_type"="out")) %>%
   filter(output_type!="CAGE" | strand=="+") %>% mutate(kb=pos/1e3, nice=factor(nice,levels=sel$nice))
rat2<-cd %>% group_by(nice) %>% summarise(pct=round(100*max(value[allele=="alt"])/max(value[allele=="ref"])),.groups="drop")
cd<-cd %>% left_join(rat2,by="nice") %>% mutate(lab=sprintf("%s (alt %d%%)",nice,pct),
        lab=factor(lab,levels=sprintf("%s (alt %d%%)",sel$nice,rat2$pct[match(sel$nice,rat2$nice)])))
pC<-ggplot(cd,aes(kb,value))+
  annotate("rect",xmin=(C-2)/1e3,xmax=(C+2)/1e3,ymin=-Inf,ymax=Inf,fill="#F4A261",alpha=0.55)+
  geom_area(data=~filter(.x,allele=="ref"),fill="grey60",colour="grey55",linewidth=0.2,alpha=0.7)+
  geom_area(data=~filter(.x,allele=="alt"),fill="#2A6FB2",colour="#2A6FB2",linewidth=0.25,alpha=0.5)+
  facet_wrap(~lab,ncol=3,scales="free_y")+
  scale_x_continuous(breaks=scales::pretty_breaks(3))+
  labs(x="chr1 position (kb) — orange = variant · grey = ref, blue = alt   (240 bp from CFH TSS)",
       y="predicted signal",title="CFH chr1:196,651,994 G>T — kidney/liver promoter loss (aFC -2.05)")+
  theme_classic(base_size=8,base_family="inter")+
  theme(strip.background=element_blank(),strip.text=element_text(size=6.9,face="bold",hjust=0),
        plot.title=element_text(size=8.5,face="bold"),axis.text=element_text(size=6.3,colour="grey20"),
        axis.title=element_text(size=7.2),panel.spacing=unit(2.6,"mm"),plot.margin=margin(4,6,4,6))
ggsave(file.path(FIGDIR,"FigureCFH_promoter.pdf"),pC,width=180,height=78,units="mm",device=cairo_pdf)
ggsave(file.path(FIGDIR,"FigureCFH_promoter.png"),pC,width=180,height=78,units="mm",dpi=300)
cat("wrote VANGL2 + CFH figures\n"); print(rat); print(rat2)
