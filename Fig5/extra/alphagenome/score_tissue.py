import os, sys
import re, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome
KEY=config.require_key(); model=dna_client.create(KEY); OUT=config.AG_DIR
cand=pd.read_csv(f"{OUT}/screen_candidates.tsv",sep="\t")
scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())

CATS={
 'blood/immune': r"blood|lymph|T-cell|T cell|B cell|monocyte|macrophage|natural killer|CD4|CD8|GM12|K562|Jurkat|myeloid|dendritic|thymus|spleen|PBMC|leukocyte|erythro|megakaryo|NAMALWA|DND-41|granulocyte|mast cell",
 'brain/neural': r"brain|neuron|neural|astrocy|glia|oligodendr|cortex|cerebell|hippocamp|spinal|nigra|SK-N|BE2C|Daoy|caudate|putamen|accumbens|pericyte|leptomening|Schwann|microglia|neuroblast|medulloblast",
 'heart/cardiac': r"heart|cardiac|ventricle|atri|myocard|cardiomyocyte",
 'vascular/SMC': r"aorta|artery|arterial|smooth muscle|endotheli|vascular|coronary|HUVEC",
 'kidney': r"kidney|renal|nephron|glomerul|podocyte|HK-2|tubule",
 'liver': r"liver|hepato|HepG2|Huh",
 'lung': r"lung|bronch|alveol|pulmonary|airway|Calu|A549|NCI-H",
 'skin/epith': r"skin|keratinoc|epiderm|melanocyte|fibroblast|foreskin",
 'GI': r"intestin|colon|stomach|gastric|sigmoid|duoden|rectum|esophag|Caco|HCT|colorect",
 'muscle': r"skeletal muscle|myocyte|myotube|myoblast|psoas",
 'reprod/other': r"testis|ovary|uter|prostate|placenta|mammary|breast|MCF|pancrea|adrenal|thyroid|pituitary|eye|retina|iris",
}
def cat_of(b):
    b=str(b)
    for c,pat in CATS.items():
        if re.search(pat,b,re.I): return c
    return 'other'

rows=[]
for i,r in cand.iterrows():
    try:
        v=genome.Variant(str(r.chrom),int(r.position),str(r.ref),str(r.alt))
        iv=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
        td=variant_scorers.tidy_scores(model.score_variant(interval=iv,variant=v,variant_scorers=scorers))
        sc=td['variant_scorer'].astype(str)
        acc=td[sc.str.contains("DIFF_LOG2_SUM") & td.output_type.astype(str).str.contains("DNASE|ATAC|CAGE")].copy()
        acc['cat']=acc.biosample_name.map(cat_of)
        acc['aq']=acc.quantile_score.abs()
        # per-category robust effect = 90th pct |quantile|
        g=acc.groupby('cat')['aq'].quantile(0.90)
        d=dict(variant=r.variant,gene=r.gene_name,aFC=round(float(r.aFC),2),maxAF=r.maxAF,subpop=r.subpop,
               pLI=round(float(r.pLI),3) if pd.notna(r.pLI) else np.nan, csq=r.consequence)
        for c in CATS: d[c]=round(float(g.get(c,np.nan)),3) if c in g and pd.notna(g.get(c)) else np.nan
        bl=d['blood/immune'] if pd.notna(d['blood/immune']) else 0
        nonblood={c:d[c] for c in CATS if c!='blood/immune' and pd.notna(d[c])}
        if nonblood:
            topc=max(nonblood,key=nonblood.get); d['top_nonblood']=topc; d['top_nonblood_q']=round(nonblood[topc],3)
            d['nonblood_minus_blood']=round(nonblood[topc]-bl,3)
        else: d['top_nonblood']=""; d['top_nonblood_q']=np.nan; d['nonblood_minus_blood']=np.nan
        rows.append(d)
        print(f"[{i+1}/{len(cand)}] {r.gene_name:10s} blood={bl:.2f} top_nonblood={d['top_nonblood']}({d['top_nonblood_q']})",flush=True)
    except Exception as e:
        print(f"[{i+1}/{len(cand)}] {r.gene_name} FAIL {e}",flush=True)
res=pd.DataFrame(rows); res.to_csv(f"{OUT}/screen_tissue.tsv",sep="\t",index=False)
print("saved screen_tissue.tsv",len(res))
