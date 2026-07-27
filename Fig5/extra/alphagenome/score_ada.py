import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome

KEY = config.require_key()
OUT = config.AG_DIR
model = dna_client.create(KEY)

VARIANTS = [
    dict(pos=44643792, ref="C", alt="T", pop="common", aFC=-0.32),
    dict(pos=44646218, ref="C", alt="T", pop="SAS",    aFC= 0.96),
    dict(pos=44646238, ref="C", alt="T", pop="EAS",    aFC=-1.00),
    dict(pos=44646509, ref="T", alt="C", pop="EUR",    aFC= 0.65),
    dict(pos=44650274, ref="C", alt="A", pop="AFR",    aFC=-0.68),
]
scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())

frames=[]
for v in VARIANTS:
    vid=f"{v['pos']}_{v['ref']}>{v['alt']}_{v['pop']}"
    var=genome.Variant("chr20", v["pos"], v["ref"], v["alt"])
    iv=var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    out=model.score_variant(interval=iv, variant=var, variant_scorers=scorers)
    td=variant_scorers.tidy_scores(out)
    td["variant_id"]=vid; td["pop"]=v["pop"]; td["aFC"]=v["aFC"]; td["pos"]=v["pos"]
    frames.append(td); print("scored", vid, len(td), flush=True)

df=pd.concat(frames, ignore_index=True)
# stringify any non-primitive object columns (Interval etc.)
for c in df.columns:
    if df[c].dtype==object and df[c].map(lambda x: not isinstance(x,(str,int,float,bool,type(None)))).any():
        df[c]=df[c].astype(str)
df.to_parquet(f"{OUT}/ada_ag_scores_full.parquet")
print("saved full parquet:", df.shape)

# ---- immune / blood filter ----
IMM=r"T-cell|T cell|B cell|natural killer|monocyte|macrophage|lymph|GM128|GM129|K562|Jurkat|CD4|CD8|dendritic|hematopoietic|thymus|leukocyte|myeloid|NAMALWA|DND-41|Karpas|OCI-LY|MM.1S|RPMI8226|NCI-H929"
bios=df["biosample_name"].fillna("")
gtex=df["gtex_tissue"].fillna("") if "gtex_tissue" in df else ""
imm_mask=bios.str.contains(IMM, case=False, regex=True) | df["gtex_tissue"].fillna("").str.contains("Blood", case=False)

# ADA gene-level RNA (LFC) across ALL tracks -> then immune subset
ada=df[(df["gene_name"]=="ADA") & (df["output_type"].astype(str).str.contains("RNA", case=False))].copy()
ada.to_csv(f"{OUT}/ada_gene_rna_all.tsv", sep="\t", index=False)
ada_imm=ada[ada["biosample_name"].fillna("").str.contains(IMM,case=False,regex=True) | ada["gtex_tissue"].fillna("").str.contains("Blood",case=False)]
print("\n=== ADA gene-level RNA LFC (immune/blood tracks) : median raw_score per variant ===")
print(ada_imm.groupby(["pos","pop","aFC"])["raw_score"].median().reset_index().to_string(index=False))

# immune epigenetic tracks (non-gene layers)
epi=df[imm_mask & df["output_type"].astype(str).str.contains("DNASE|ATAC|CHIP|CAGE", case=False, regex=True)].copy()
epi.to_csv(f"{OUT}/ada_immune_epigenetic.tsv", sep="\t", index=False)
print("\nimmune epigenetic rows:", len(epi), "| ADA RNA immune rows:", len(ada_imm))

# GM12878 focused snapshot
gm=df[df["biosample_name"].fillna("").str.contains("GM12878")].copy()
gm.to_csv(f"{OUT}/ada_gm12878.tsv", sep="\t", index=False)
print("GM12878 rows:", len(gm))
