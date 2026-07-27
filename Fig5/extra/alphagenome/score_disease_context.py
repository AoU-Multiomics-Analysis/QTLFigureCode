import os, sys
import re, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome
KEY=config.require_key(); model=dna_client.create(KEY); OUT=config.AG_DIR
scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())

CATS={
 'brain/neural': r"brain|neuron|neural|astrocy|glia|oligodendr|cortex|cerebell|hippocamp|spinal|nigra|SK-N|BE2C|Daoy|caudate|putamen|accumbens|pericyte|leptomening|Schwann|microglia|neuroblast|medulloblast",
 'vascular/SMC': r"aorta|artery|arterial|smooth muscle|endotheli|vascular|coronary|HUVEC|pericyte",
 'heart/cardiac': r"heart|cardiac|ventricle|atri|myocard|cardiomyocyte",
 'skin/epith': r"skin|keratinoc|epiderm|melanocyte|foreskin|esophag|epithelial",
 'liver': r"liver|hepato|HepG2|Huh",
 'kidney': r"kidney|renal|nephron|glomerul|podocyte|HK-2|tubule",
 'lung': r"lung|bronch|alveol|pulmonary|airway|Calu|A549|NCI-H",
}
# variant -> (chrom,pos,ref,alt,gene,aFC, [disease categories])
CAND=[
 ("chr17",39539868,"C","A","PPP1R1B", 2.36,["brain/neural"]),
 ("chr10",88990820,"AC","A","ACTA2", -1.86,["vascular/SMC","heart/cardiac"]),
 ("chr6", 7609474,"T","C","DSP",     3.09,["heart/cardiac","skin/epith"]),
 ("chr1", 196651994,"G","T","CFH",  -2.05,["liver","kidney"]),
 ("chr1", 160428770,"G","C","VANGL2",6.05,["brain/neural","skin/epith"]),
 ("chr12",51912316,"C","T","ACVRL1",-6.64,["vascular/SMC","lung"]),
 ("chr8", 142611326,"G","T","ADGRB1",1.61,["brain/neural"]),
 ("chr5", 180675167,"C","G","FLT4",  1.52,["vascular/SMC","lung"]),
 ("chr1", 24304391,"C","T","GRHL3",  2.91,["skin/epith"]),
 ("chr12",12267039,"G","T","LRP6",   1.50,["vascular/SMC","liver"]),
]
def catmatch(b,cats):
    b=str(b)
    for c in cats:
        if re.search(CATS[c],b,re.I): return True
    return False

rows=[]
for chrom,pos,ref,alt,gene,aFC,cats in CAND:
    try:
        v=genome.Variant(chrom,pos,ref,alt); iv=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
        td=variant_scorers.tidy_scores(model.score_variant(interval=iv,variant=v,variant_scorers=scorers))
        sc=td['variant_scorer'].astype(str); ot=td['output_type'].astype(str)
        dm=td['biosample_name'].apply(lambda b: catmatch(b,cats))
        diff=td[sc.str.contains("DIFF_LOG2_SUM")]
        def eff(mask):
            s=diff[mask & dm]
            if not len(s): return (np.nan,np.nan,"")
            i=s['raw_score'].abs().idxmax(); r=s.loc[i]
            return (round(s['raw_score'].mean(),2), round(float(r['raw_score']),2), str(r['biosample_name'])[:20])
        acc=eff(ot.str.contains("DNASE|ATAC")); cage=eff(ot.str.contains("CAGE"))
        # histone: top signed
        h=diff[ot.str.contains("HISTONE") & dm]
        hmark,hval="",np.nan
        if len(h): i=h['raw_score'].abs().idxmax(); hmark=h.loc[i,'histone_mark']; hval=round(float(h.loc[i,'raw_score']),2)
        # TF: top signed
        tf=diff[ot.str.contains("CHIP_TF") & dm]
        tfn,tfv="",np.nan
        if len(tf): i=tf['quantile_score'].abs().idxmax(); tfn=tf.loc[i,'transcription_factor']; tfv=round(float(tf.loc[i,'quantile_score']),2)
        # transcriptomics: RNA baseline (Active) + LFC for the gene in disease tissue
        act=td[sc.str.contains("GeneMaskActiveScorer") & (td.gene_name==gene) & dm]['raw_score']
        lfc=td[sc.str.contains("GeneMaskLFCScorer") & (td.gene_name==gene) & dm]['raw_score']
        base=round(float(act.max()),2) if len(act) else np.nan
        rlfc=round(float(lfc.median()),3) if len(lfc) else np.nan
        rows.append(dict(gene=gene,aFC=aFC,disease_tissue="+".join(cats),
            acc_mean=acc[0],acc_top=acc[1],acc_where=acc[2],
            cage_mean=cage[0], rna_base=base, rna_lfc=rlfc,
            top_histone=f"{hmark}{hval:+.1f}" if hmark else "",
            top_TF=f"{tfn}({tfv:+.2f})" if tfn else ""))
        print(f"{gene:8s} {cats[0]:14s} acc{acc[0]}/{acc[1]}@{acc[2]:18s} cage{cage[0]} RNAbase{base} lfc{rlfc} {hmark}{hval:+.1f} TF={tfn}{tfv}",flush=True)
    except Exception as e:
        print(gene,"FAIL",e,flush=True)
res=pd.DataFrame(rows); res.to_csv(f"{OUT}/disease_context.tsv",sep="\t",index=False)
print("\nsaved disease_context.tsv")
