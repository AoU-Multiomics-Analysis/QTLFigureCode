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
CURIES=["CL:0002603","CL:0002605","CL:0002453","CL:0002573","CL:0000540","CL:0000047",
        "CL:0000708","CL:0002575","CL:0002590","CL:0000706","UBERON:0002037","UBERON:0000955"]
vo=model.predict_variant(interval=iv, variant=var, requested_outputs=[OT.CAGE], ontology_terms=CURIES)
WIN=1600; lo,hi=VPOS-WIN,VPOS+WIN; rows=[]
for td,al in [(vo.reference.cage,'ref'),(vo.alternate.cage,'alt')]:
    vals=np.asarray(td.values); pos=td.interval.start+np.arange(vals.shape[0])*td.resolution
    keep=np.where((pos>=lo)&(pos<=hi))[0]; md=td.metadata.reset_index(drop=True)
    for j in range(vals.shape[1]):
        m=md.iloc[j]
        for k in keep:
            rows.append((al,int(pos[k]),float(vals[k,j]),m.get('biosample_name',''),str(m.get('strand','.'))))
df=pd.DataFrame(rows,columns=['allele','pos','value','biosample','strand'])
df.to_parquet(f"{OUT}/spr_cage_celltype.parquet")
# summit reduction (+ strand) per cell type
print("cell-type CAGE at TSS summit (+ strand): ref, alt, %ref")
for bs in sorted(df.biosample.unique()):
    r=df[(df.biosample==bs)&(df.allele=='ref')&(df.strand=='+')]; a=df[(df.biosample==bs)&(df.allele=='alt')&(df.strand=='+')]
    if len(r)==0: continue
    rr=r.value.max(); aa=a.value.max()
    print(f"  {bs:44s} ref={rr:6.1f} alt={aa:6.1f}  ({100*aa/max(rr,.1):3.0f}%)")
