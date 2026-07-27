import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from alphagenome.models import dna_client
from alphagenome.data import genome
KEY=config.require_key(); model=dna_client.create(KEY); OT=dna_client.OutputType
OUT=config.AG_DIR; VPOS=72887998
var=genome.Variant("chr2", VPOS, "GCGCC", "G")
iv=var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_500KB)
CURIES=["CL:0002603","CL:0002606","CL:0002604","CL:0000127","CL:0000679","CL:0011020","CL:0000047",  # cell types
        "UBERON:0002037","UBERON:0002038","UBERON:0000955"]                                          # bulk contrast
outs=[OT.DNASE, OT.CAGE]
vo=model.predict_variant(interval=iv, variant=var, requested_outputs=outs, ontology_terms=CURIES)
WIN=1600; lo,hi=VPOS-WIN,VPOS+WIN; rows=[]
for attr,otype in [('dnase','DNASE'),('cage','CAGE')]:
    tr=getattr(vo.reference,attr,None); ta=getattr(vo.alternate,attr,None)
    if tr is None: continue
    for td,al in [(tr,'ref'),(ta,'alt')]:
        vals=np.asarray(td.values); pos=td.interval.start+np.arange(vals.shape[0])*td.resolution
        keep=np.where((pos>=lo)&(pos<=hi))[0]; md=td.metadata.reset_index(drop=True)
        for j in range(vals.shape[1]):
            m=md.iloc[j]
            for k in keep:
                rows.append((otype,al,int(pos[k]),float(vals[k,j]),m.get('biosample_name',''),
                             str(m.get('strand','.')),str(m.get('ontology_curie',''))))
df=pd.DataFrame(rows,columns=['output_type','allele','pos','value','biosample','strand','ontology'])
df.to_parquet(f"{OUT}/spr_tracks_celltype.parquet")
print("saved",df.shape); print(sorted(df.biosample.unique()))
