import os, sys
import pandas as pd, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
pd.set_option("display.width",200); pd.set_option("display.max_columns",30)
OUT=config.AG_DIR
df=pd.read_parquet(f"{OUT}/spr_ag_scores_full.parquet"); df["ot"]=df["output_type"].astype(str)

rna=df[(df.gene_name=="SPR") & df.ot.str.contains("RNA")]
print("=== SPR gene-level RNA LFC (alt vs ref) across tissues ===")
print("n tissues:", len(rna), "| mean LFC %.3f  min %.3f  max %.3f  median %.3f"%(
    rna.raw_score.mean(), rna.raw_score.min(), rna.raw_score.max(), rna.raw_score.median()))
wb=rna[rna['gtex_tissue'].fillna('').str.contains('Blood')]
print("GTEx Whole Blood LFC:", wb['raw_score'].tolist())
print("\nMost negative SPR RNA tracks:")
print(rna.sort_values('raw_score').head(8)[['biosample_name','gtex_tissue','Assay title','raw_score','quantile_score']].to_string(index=False))

acc=df[df.ot.str.contains("DNASE|ATAC") & (df.quantile_score.abs()>=0.9)]
print(f"\n=== Accessibility tracks with strong effect: {len(acc)} ===")
print("sign split: neg=%d pos=%d"%((acc.quantile_score<0).sum(),(acc.quantile_score>0).sum()))
acc2=acc.reindex(acc.quantile_score.abs().sort_values(ascending=False).index).head(15)
print(acc2[['ot','biosample_name','raw_score','quantile_score']].to_string(index=False))

# overall directionality summary by layer
print("\n=== net direction by layer (strong-effect tracks, |q|>=0.9) ===")
s=df[df.quantile_score.abs()>=0.9].copy()
def cat(o):
    if 'DNASE' in o or 'ATAC' in o: return 'accessibility'
    if 'CAGE' in o or 'PROCAP' in o: return 'TSS (CAGE/PRO-cap)'
    if 'HISTONE' in o: return 'histone'
    if 'CHIP_TF' in o: return 'TF ChIP'
    if 'RNA' in o: return 'RNA'
    if 'SPLICE' in o or 'JUNC' in o: return 'splice'
    return o
s['layer']=s.ot.map(cat)
print(s.groupby('layer')['quantile_score'].agg(n='size', mean='mean',
      pct_down=lambda x:(x<0).mean()*100).round(3).to_string())
