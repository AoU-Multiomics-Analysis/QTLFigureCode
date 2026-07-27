import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from alphagenome.models import dna_client
from alphagenome.data import genome
KEY=config.require_key(); model=dna_client.create(KEY); OT=dna_client.OutputType
OUT=config.AG_DIR
VPOS=72887998; TSS=72887382
var=genome.Variant("chr2", VPOS, "GCGCC", "G")
iv=var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_500KB)   # 500kb context for accuracy
CURIES=["UBERON:0000955","UBERON:0002626","UBERON:0001874","CL:0002603","CL:0000679",
        "UBERON:0009834","UBERON:0002240","UBERON:0013756","CL:2000001","EFO:0002784","EFO:0002067"]
outs=[OT.DNASE, OT.ATAC, OT.CAGE, OT.CHIP_HISTONE, OT.CHIP_TF, OT.RNA_SEQ]
vo=model.predict_variant(interval=iv, variant=var, requested_outputs=outs, ontology_terms=CURIES)

WIN=1600; lo, hi = VPOS-WIN, VPOS+WIN
rows=[]
for attr,otype in [('dnase','DNASE'),('atac','ATAC'),('cage','CAGE'),
                   ('chip_histone','CHIP_HISTONE'),('chip_tf','CHIP_TF'),('rna_seq','RNA_SEQ')]:
    tref=getattr(vo.reference, attr, None); talt=getattr(vo.alternate, attr, None)
    if tref is None: continue
    for td,allele in [(tref,'ref'),(talt,'alt')]:
        vals=np.asarray(td.values); res=td.resolution; st=td.interval.start
        pos=st+np.arange(vals.shape[0])*res
        keep=np.where((pos>=lo)&(pos<=hi))[0]
        md=td.metadata.reset_index(drop=True)
        for j in range(vals.shape[1]):
            m=md.iloc[j]
            for k in keep:
                rows.append((otype, allele, int(pos[k]), float(vals[k,j]),
                             m.get('biosample_name',''), str(m.get('strand','.')),
                             str(m.get('histone_mark','')), str(m.get('transcription_factor','')),
                             str(m.get('Assay title','')), str(m.get('ontology_curie',''))))
df=pd.DataFrame(rows, columns=['output_type','allele','pos','value','biosample','strand',
                               'histone_mark','tf','assay','ontology'])
df.to_parquet(f"{OUT}/spr_tracks_window.parquet")
df.to_csv(f"{OUT}/spr_tracks_window.tsv.gz", sep="\t", index=False)
print("saved", df.shape, "window", lo, hi, "res per otype:")
print(df.groupby('output_type')['pos'].agg(['min','max','nunique']).to_string())
print("biosamples:", sorted(df.biosample.unique()))
