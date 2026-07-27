import os, sys
import pandas as pd, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
pd.set_option("display.width",200); pd.set_option("display.max_columns",30)
OUT=config.AG_DIR
df=pd.read_parquet(f"{OUT}/ada_ag_scores_full.parquet")
order=[44643792,44646218,44646238,44646509,44650274]
lab={44643792:"common -0.32",44646218:"SAS +0.96",44646238:"EAS -1.00",44646509:"EUR +0.65",44650274:"AFR -0.68"}

IMM=r"T-cell|T cell|B cell|natural killer|monocyte|macrophage|lymph|GM128|K562|Jurkat|CD4|CD8|dendritic|hematopoietic|thymus|leukocyte|myeloid|NAMALWA|DND-41|Karpas|OCI-LY|MM.1S|RPMI8226|NCI-H929"
df["imm"]=df["biosample_name"].fillna("").str.contains(IMM,case=False,regex=True)
df["ot"]=df["output_type"].astype(str)

def blk(title): print("\n"+"="*70+"\n"+title+"\n"+"="*70)

# 1) accessibility: strongest signed immune effect per variant (DNASE+ATAC)
acc=df[df.imm & df.ot.str.contains("DNASE|ATAC")]
blk("Accessibility (immune) — most extreme signed quantile per variant")
for p in order:
    d=acc[acc.pos==p]
    if len(d)==0: print(f"{lab[p]:14s}: none"); continue
    i=d["quantile_score"].abs().idxmax()
    r=d.loc[i]
    print(f"{lab[p]:14s}: {r['ot']:6s} q={r['quantile_score']:+.3f} raw={r['raw_score']:+.3f}  {r['biosample_name']}")

# 2) ETS family TF ChIP (immune), signed effect per variant
ets=df[df.imm & (df.ot.str.contains("CHIP_TF")) & df["transcription_factor"].fillna("").str.contains(r"^(GABPA|GABPB1|ELF1|ELF4|ETS1|ETV6|SPI1|FLI1|ELK1|ETV1)$",regex=True)]
blk("ETS-family TF ChIP (immune) — mean signed quantile per variant")
piv=ets.groupby(["pos"]).apply(lambda d: pd.Series({
    "mean_q":d["quantile_score"].mean(),
    "min_q":d["quantile_score"].min(),
    "max_q":d["quantile_score"].max(),
    "n":len(d)})).reset_index()
piv["lab"]=piv.pos.map(lab); print(piv[["lab","mean_q","min_q","max_q","n"]].to_string(index=False))

# 3) per-variant top immune TF hits (any TF) by |quantile|
blk("Top immune TF-ChIP hits per variant (|quantile|)")
tf=df[df.imm & df.ot.str.contains("CHIP_TF")]
for p in order:
    d=tf[tf.pos==p].copy()
    if len(d)==0: continue
    d["aq"]=d["quantile_score"].abs()
    top=d.sort_values("aq",ascending=False).head(6)
    hits="; ".join(f"{t}({q:+.2f})" for t,q in zip(top["transcription_factor"],top["quantile_score"]))
    print(f"{lab[p]:14s}: {hits}")

# 4) GM12878 specific: accessibility + H3K27ac + key ETS
blk("GM12878 — accessibility, H3K27ac, ETS TFs (signed quantile)")
gm=df[df["biosample_name"].fillna("").str.contains("GM12878")]
for p in order:
    d=gm[gm.pos==p]
    def q(sel):
        s=d[sel]; return s["quantile_score"].mean() if len(s) else np.nan
    dn=q(d.ot.str.contains("DNASE")); at=q(d.ot.str.contains("ATAC"))
    k27=q(d["histone_mark"].fillna("").eq("H3K27ac"))
    etsq=q(d.ot.str.contains("CHIP_TF") & d["transcription_factor"].fillna("").isin(["GABPA","ELF1","ETS1","ETV6","SPI1"]))
    print(f"{lab[p]:14s}: DNase q={dn:+.3f}  ATAC q={at:+.3f}  H3K27ac q={k27:+.3f}  ETS q={etsq:+.3f}")

# 5) ADA RNA gene-level: max |LFC| immune track + gtex whole blood
blk("ADA RNA — immune-track mean LFC and GTEx Whole Blood")
ada=df[(df.gene_name=="ADA") & df.ot.str.contains("RNA")]
for p in order:
    d=ada[ada.pos==p]
    imm=d[d.imm]["raw_score"].mean()
    wb=d[d["gtex_tissue"].fillna("").str.contains("Blood")]["raw_score"]
    wbv=wb.mean() if len(wb) else np.nan
    print(f"{lab[p]:14s}: immune meanLFC={imm:+.4f}  GTEx-WholeBlood LFC={wbv:+.4f}  (n_imm={d['imm'].sum()})")
