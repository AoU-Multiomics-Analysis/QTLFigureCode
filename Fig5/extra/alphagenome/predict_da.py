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
# catecholaminergic/dopaminergic-lineage neuroblastoma + glia & neuron reference
CURIES=["EFO:0002860","EFO:0003072","EFO:0005725","EFO:0005721",   # SK-N-MC, SK-N-SH, BE2C, SK-N-DZ
        "CL:0002603","CL:0002605","CL:0000540","CL:0000047","UBERON:0002038"]  # cerebellar+cortical astro, neuron, NSC, substantia nigra
vo=model.predict_variant(interval=iv, variant=var, requested_outputs=[OT.CAGE,OT.DNASE], ontology_terms=CURIES)
WIN=1600; lo,hi=VPOS-WIN,VPOS+WIN; rows=[]
for attr,ot in [('cage','CAGE'),('dnase','DNASE')]:
    for td,al in [(getattr(vo.reference,attr),'ref'),(getattr(vo.alternate,attr),'alt')]:
        vals=np.asarray(td.values); pos=td.interval.start+np.arange(vals.shape[0])*td.resolution
        keep=np.where((pos>=lo)&(pos<=hi))[0]; md=td.metadata.reset_index(drop=True)
        for j in range(vals.shape[1]):
            m=md.iloc[j]
            for k in keep: rows.append((ot,al,int(pos[k]),float(vals[k,j]),m.get('biosample_name',''),str(m.get('strand','.'))))
df=pd.DataFrame(rows,columns=['output_type','allele','pos','value','biosample','strand'])
df.to_parquet(f"{OUT}/spr_da_lineage.parquet")
print("=== SPR element in dopaminergic-lineage vs glia vs neuron ===")
for ot,strand in [('CAGE','+'),('DNASE','.')]:
    print(f"\n-- {ot} summit (ref -> alt %ref) --")
    for bs in sorted(df.biosample.unique()):
        r=df[(df.output_type==ot)&(df.biosample==bs)&(df.strand==strand)&(df.allele=='ref')]
        a=df[(df.output_type==ot)&(df.biosample==bs)&(df.strand==strand)&(df.allele=='alt')]
        if len(r)==0: continue
        rr,aa=r.value.max(),a.value.max()
        print(f"  {bs:32s} ref={rr:6.1f} alt={aa:6.1f} ({100*aa/max(rr,.1):3.0f}%)")
