import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from alphagenome.models import dna_client
from alphagenome.data import genome
KEY=config.require_key(); model=dna_client.create(KEY)
OT=dna_client.OutputType
VPOS=72887998; TSS=72887382
var=genome.Variant("chr2", VPOS, "GCGCC", "G")
iv=var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)   # 16kb context, fast, local promoter well inside
print("interval:", iv, "width", iv.width)

CURIES=["UBERON:0000955","UBERON:0002626","UBERON:0001874","CL:0002603","CL:0000679",
        "UBERON:0009834","UBERON:0002240","UBERON:0013756","CL:2000001","EFO:0002784","EFO:0002067"]
outs=[OT.DNASE, OT.ATAC, OT.CAGE, OT.CHIP_HISTONE, OT.CHIP_TF, OT.RNA_SEQ]
vo=model.predict_variant(interval=iv, variant=var, requested_outputs=outs, ontology_terms=CURIES)
print("VariantOutput attrs:", [a for a in dir(vo) if not a.startswith('_')])
ref, alt = vo.reference, vo.alternate
print("reference attrs:", [a for a in dir(ref) if not a.startswith('_')])
for name in ['dnase','atac','cage','chip_histone','chip_tf','rna_seq']:
    td=getattr(ref, name, None)
    if td is None: print(name, "-> None"); continue
    print(f"\n=== {name} === values shape {np.asarray(td.values).shape} resolution {td.resolution}")
    print("  interval:", td.interval)
    md=td.metadata
    print("  metadata cols:", list(md.columns))
    print(md[[c for c in ['name','biosample_name','strand','histone_mark','transcription_factor','ontology_curie'] if c in md.columns]].to_string(index=False))
