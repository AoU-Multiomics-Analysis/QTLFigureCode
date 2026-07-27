import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome
pd.set_option("display.width",220); pd.set_option("display.max_columns",40)
KEY=config.require_key(); OUT=config.AG_DIR
model=dna_client.create(KEY)

# SPR promoter-proximal rare deletion (GRCh38), measured whole-blood eQTL log2 aFC = -2.76
var=genome.Variant("chr2", 72887998, "GCGCC", "G")   # 4-bp deletion of a GC-box
iv=var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
out=model.score_variant(interval=iv, variant=var, variant_scorers=scorers)
df=variant_scorers.tidy_scores(out)
for c in df.columns:
    if df[c].dtype==object and df[c].map(lambda x:not isinstance(x,(str,int,float,bool,type(None)))).any():
        df[c]=df[c].astype(str)
df["ot"]=df["output_type"].astype(str)
df.to_parquet(f"{OUT}/spr_ag_scores_full.parquet")
print("total track-rows:", len(df), "| output types:", sorted(df.ot.unique()))

# ---- SPR gene-level RNA (LFC) ----
rna=df[(df.gene_name=="SPR") & df.ot.str.contains("RNA")]
print("\n=== SPR gene-level RNA LFC (alt vs ref) : distribution across tissues ===")
print(rna["raw_score"].describe()[["mean","min","max"]].to_string())
top=rna.reindex(rna["raw_score"].abs().sort_values(ascending=False).index).head(10)
print(top[["biosample_name","gtex_tissue","Assay title","raw_score","quantile_score"]].to_string(index=False))

# ---- 'actual effects' across ALL tracks: |quantile|>=0.9 ----
TH=0.90
eff=df[df["quantile_score"].abs()>=TH].copy()
print(f"\n=== tracks with strong effect (|quantile|>= {TH}): {len(eff)} of {len(df)} ===")
print(eff.groupby("ot").size().sort_values(ascending=False).to_string())

def layer(name, mask, extra):
    d=eff[mask].copy()
    if not len(d): print(f"\n--- {name}: none ---"); return
    d["aq"]=d["quantile_score"].abs()
    d=d.sort_values("aq",ascending=False).head(20)
    print(f"\n--- {name}  (top {len(d)}; signed quantile) ---")
    print(d[extra+["biosample_name","raw_score","quantile_score"]].to_string(index=False))

layer("Chromatin accessibility (DNase/ATAC)", eff.ot.str.contains("DNASE|ATAC"), ["ot","Assay title"])
layer("TSS activity (CAGE / PRO-cap)", eff.ot.str.contains("CAGE|PROCAP|PRO_CAP|PROCap|PRO-cap|PROcap"), ["ot","Assay title"])
layer("Histone ChIP", eff.ot.str.contains("HISTONE"), ["histone_mark"])
layer("TF ChIP", eff.ot.str.contains("CHIP_TF"), ["transcription_factor"])

# TF summary: which factors most often show strong effect and net direction
tf=eff[eff.ot.str.contains("CHIP_TF")]
if len(tf):
    g=tf.groupby("transcription_factor")["quantile_score"].agg(["mean","min","max","count"]).sort_values("count",ascending=False)
    print("\n=== TF ChIP factors with strong effects (by frequency; net signed) ===")
    print(g.head(25).to_string())
