import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome

KEY = config.require_key()
model = dna_client.create(KEY)

v = genome.Variant(chromosome="chr20", position=44646238,
                   reference_bases="C", alternate_bases="T")
interval = v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
print("interval:", interval)

scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
print("n recommended scorers:", len(scorers))
for name, s in variant_scorers.RECOMMENDED_VARIANT_SCORERS.items():
    print("  scorer:", name, "->", type(s).__name__)

out = model.score_variant(interval=interval, variant=v, variant_scorers=scorers)
print("\nscore_variant returned:", type(out), "len", len(out))
for i, ad in enumerate(out):
    print(f"\n--- output {i}: {type(ad).__name__} shape={getattr(ad,'shape',None)} ---")
    try:
        print("  var columns:", list(ad.var.columns))
        print("  uns keys:", list(ad.uns.keys()) if hasattr(ad,'uns') else None)
        print("  var head:\n", ad.var.head(3).to_string())
        import numpy as np
        print("  X dtype/shape:", ad.X.shape, "sample:", np.asarray(ad.X).ravel()[:3])
    except Exception as e:
        print("  introspect err:", repr(e))
