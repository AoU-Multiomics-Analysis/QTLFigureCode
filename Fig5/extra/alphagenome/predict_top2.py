import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from alphagenome.models import dna_client
from alphagenome.data import genome
KEY=config.require_key(); model=dna_client.create(KEY); OT=dna_client.OutputType
OUT=config.AG_DIR
def run(tag, chrom,pos,ref,alt, outs, curies, keep_tf=None, win=2000):
    v=genome.Variant(chrom,pos,ref,alt); iv=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_100KB)
    vo=model.predict_variant(interval=iv, variant=v, requested_outputs=outs, ontology_terms=curies)
    lo,hi=pos-win,pos+win; rows=[]
    amap={'dnase':'DNASE','atac':'ATAC','cage':'CAGE','chip_tf':'CHIP_TF'}
    for attr,ot in amap.items():
        tr=getattr(vo.reference,attr,None); ta=getattr(vo.alternate,attr,None)
        if tr is None: continue
        for td,al in [(tr,'ref'),(ta,'alt')]:
            vals=np.asarray(td.values); p=td.interval.start+np.arange(vals.shape[0])*td.resolution
            keep=np.where((p>=lo)&(p<=hi))[0]; md=td.metadata.reset_index(drop=True)
            for j in range(vals.shape[1]):
                m=md.iloc[j]; tf=str(m.get('transcription_factor',''))
                if ot=='CHIP_TF' and keep_tf and tf not in keep_tf: continue
                for k in keep:
                    rows.append((ot,al,int(p[k]),float(vals[k,j]),m.get('biosample_name',''),str(m.get('strand','.')),tf))
    df=pd.DataFrame(rows,columns=['output_type','allele','pos','value','biosample','strand','tf'])
    df.to_parquet(f"{OUT}/tracks_{tag}.parquet"); print(tag,"saved",df.shape,"bios",sorted(df.biosample.unique())[:6])

# VANGL2 insulator: CTCF+RAD21 + accessibility in neural cells
run("VANGL2", "chr1",160428770,"G","C",
    [OT.CHIP_TF, OT.DNASE, OT.ATAC],
    ["CL:0000706","CL:0000679","UBERON:0000955","CL:0002606","CL:2000044","CL:0002319","CL:0000047"],
    keep_tf={"CTCF","RAD21","SMC3","STAG1"})
# CFH promoter: CAGE + DNase + Pol II in kidney/liver
run("CFH", "chr1",196651994,"G","T",
    [OT.CAGE, OT.DNASE, OT.CHIP_TF],
    ["CL:0002584","CL:0002306","CL:1000494","CL:0002188","CL:0002518","EFO:0001187","CL:0000182","UBERON:0002107"],
    keep_tf={"POLR2A","POLR2G","CTCF","GABPB1","IRF2"})
