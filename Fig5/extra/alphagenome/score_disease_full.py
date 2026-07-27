import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome
KEY=config.require_key(); model=dna_client.create(KEY); OUT=config.AG_DIR
scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
CAND=[
 ("chr17",39539868,"C","A","PPP1R1B", 2.36),("chr10",88990820,"AC","A","ACTA2",-1.86),
 ("chr6", 7609474,"T","C","DSP",3.09),("chr1",196651994,"G","T","CFH",-2.05),
 ("chr1",160428770,"G","C","VANGL2",6.05),("chr12",51912316,"C","T","ACVRL1",-6.64),
 ("chr8",142611326,"G","T","ADGRB1",1.61),("chr5",180675167,"C","G","FLT4",1.52),
 ("chr1",24304391,"C","T","GRHL3",2.91),("chr12",12267039,"G","T","LRP6",1.50),
]
keep=['gene_name','output_type','variant_scorer','track_name','biosample_name',
      'transcription_factor','histone_mark','Assay title','raw_score','quantile_score']
frames=[]
for chrom,pos,ref,alt,gene,aFC in CAND:
    v=genome.Variant(chrom,pos,ref,alt); iv=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    td=variant_scorers.tidy_scores(model.score_variant(interval=iv,variant=v,variant_scorers=scorers))
    sc=td['variant_scorer'].astype(str)
    # keep DIFF (epigenetic/CAGE) rows + gene-level rows for the target gene
    sub=td[sc.str.contains("DIFF_LOG2_SUM") | ((sc.str.contains("GeneMask")) & (td.gene_name==gene))].copy()
    sub=sub[[c for c in keep if c in sub.columns]]
    sub['V_gene']=gene; sub['V_aFC']=aFC
    frames.append(sub); print("scored",gene,len(sub),flush=True)
alld=pd.concat(frames,ignore_index=True)
alld.to_parquet(f"{OUT}/disease_context_full.parquet")
print("saved disease_context_full.parquet",alld.shape)
