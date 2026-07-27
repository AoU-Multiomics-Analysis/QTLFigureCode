import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome
KEY=config.require_key(); model=dna_client.create(KEY); OUT=config.AG_DIR
cand=pd.read_csv(f"{OUT}/screen_candidates.tsv",sep="\t")
scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())
rows=[]
for i,r in cand.iterrows():
    try:
        v=genome.Variant(str(r.chrom), int(r.position), str(r.ref), str(r.alt))
        iv=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
        td=variant_scorers.tidy_scores(model.score_variant(interval=iv, variant=v, variant_scorers=scorers))
        sc=td['variant_scorer'].astype(str); ot=td['output_type'].astype(str)
        diff=td[sc.str.contains("DIFF_LOG2_SUM")]; lfc=td[sc.str.contains("GeneMaskLFCScorer")]
        # gene-level RNA LFC for target gene
        gl=lfc[lfc.gene_name==r.gene_name]['raw_score']
        rna=gl.median() if len(gl) else np.nan
        # strongest signed effect across DIFF layers (quantile)
        d=diff.copy(); d['aq']=d['quantile_score'].abs()
        best=d.loc[d['aq'].idxmax()] if len(d) else None
        def layer(o):
            o=str(o)
            return ('access' if ('DNASE'in o or 'ATAC'in o) else 'CAGE' if('CAGE'in o or 'PROCAP'in o)
                    else 'TF' if 'CHIP_TF'in o else 'histone' if 'HISTONE'in o else o)
        acc=d[d.output_type.astype(str).str.contains("DNASE|ATAC")]
        accq=acc['quantile_score'].abs().max() if len(acc) else np.nan
        cage=d[d.output_type.astype(str).str.contains("CAGE")]
        cageq=cage['quantile_score'].abs().max() if len(cage) else np.nan
        tf=d[d.output_type.astype(str).str.contains("CHIP_TF")].copy()
        toptf=""; toptfq=np.nan
        if len(tf):
            t=tf.loc[tf['quantile_score'].abs().idxmax()]; toptf=t['transcription_factor']; toptfq=t['quantile_score']
        conc = (not np.isnan(rna)) and (np.sign(rna)==np.sign(r.aFC))
        rows.append(dict(variant=r.variant, gene=r.gene_name, aFC=round(float(r.aFC),2), maxAF=r.maxAF, subpop=r.subpop,
            pLI=round(float(r.pLI),3) if pd.notna(r.pLI) else np.nan, csq=r.consequence,
            rna_lfc=round(float(rna),3) if pd.notna(rna) else np.nan,
            best_layer=layer(best['output_type']) if best is not None else "",
            best_q=round(float(best['quantile_score']),3) if best is not None else np.nan,
            best_track=str(best['biosample_name']) if best is not None else "",
            accq=round(float(accq),3) if pd.notna(accq) else np.nan,
            cageq=round(float(cageq),3) if pd.notna(cageq) else np.nan,
            top_tf=toptf, top_tf_q=round(float(toptfq),3) if pd.notna(toptfq) else np.nan,
            concordant=conc))
        print(f"[{i+1}/{len(cand)}] {r.gene_name:10s} aFC{float(r.aFC):+.1f} rna{rna:+.2f} acc_q{accq:.2f} tf={toptf}", flush=True)
    except Exception as e:
        print(f"[{i+1}/{len(cand)}] {r.gene_name} FAILED: {e}", flush=True)
res=pd.DataFrame(rows); res.to_csv(f"{OUT}/screen_results.tsv",sep="\t",index=False)
print("\nsaved screen_results.tsv:", len(res))
