#!/usr/bin/env python3
"""Broaden chromatin with histones and add a deduplicated membrane/signaling union."""
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--paper-root', type=Path, required=True)
parser.add_argument('--prepared-paper-root', type=Path, required=True)
args = parser.parse_args()
ROOT = args.paper_root.resolve()
PREPARED = args.prepared_paper_root.resolve()
if ROOT == PREPARED:
    raise ValueError('Prepared outputs must use a separate directory.')
PARENT = ROOT / 'output/figure1_novel_gene_categories'
BASE = PREPARED / 'output/figure1_novel_gene_categories/reconciled_20261027'
OUT = PREPARED / 'output/figure1_novel_gene_categories/seven_categories_20261027'
INVESTIGATION = ROOT / 'output/uncategorized_novel_egenes_20261027'
OUT.mkdir(parents=True, exist_ok=True)
CATEGORIES = ['cancer_driver', 'mendelian_disease', 'lof_constrained', 'transcription_chromatin',
              'rna_processing_splicing', 'plasma_membrane_signaling', 'drug_target']
ASSAYS = ['eQTL', 'sQTL', 'pQTL']


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(name, rows, fields=None):
    with (OUT / name).open('w') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=fields or list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def flag(value):
    return 'TRUE' if value else 'FALSE'


membership = read(BASE / 'novel_coding_gene_annotations.tsv')
assert Counter(r['assay'] for r in membership) == {'eQTL':1207, 'sQTL':2035, 'pQTL':1143}
genes = {r['gene_id']:r['gene_name'] for r in membership}
assert len(genes) == 4146
hgnc_path = PARENT / 'annotations/hgnc_complete_set_2026-10-15.txt'
hgnc = [r for r in read(hgnc_path) if r['status'] == 'Approved']
by_id, by_symbol, protein_hgnc = defaultdict(list), {}, defaultdict(set)
gencode_path = PREPARED / 'output/qtl_circos_simplified_options/gene_reference.tsv'
gencode_symbols = defaultdict(set)
for row in read(gencode_path):
    gencode_symbols[row['gene_name']].add(row['gene_id'])
gencode_symbols = {symbol:next(iter(ids)) for symbol,ids in gencode_symbols.items() if len(ids)==1}
for row in hgnc:
    if row['ensembl_gene_id'] not in ('', 'NA'):
        by_id[row['ensembl_gene_id']].append(row)
    assert row['symbol'] not in by_symbol
    by_symbol[row['symbol']] = row
    for protein in row['uniprot_ids'].split('|'):
        if protein not in ('', 'NA'): protein_hgnc[protein].add(row['hgnc_id'])
identity, hgnc_genes, histones = [], defaultdict(set), set()
for gene, symbol in sorted(genes.items()):
    candidates = by_id.get(gene, [])
    if len(candidates) == 1:
        chosen, method = candidates[0], 'Exact Ensembl ID'
    elif not candidates and symbol in by_symbol:
        chosen, method = by_symbol[symbol], 'Unique HGNC approved symbol'
    else:
        chosen, method = {}, 'No unambiguous HGNC mapping'
    groups = [x.strip() for x in chosen.get('gene_group', '').split('|')]
    qualifying = [x for x in groups if re.fullmatch(r'H[1-4](?:A|B)? histones', x)]
    if qualifying: histones.add(gene)
    if chosen: hgnc_genes[chosen['hgnc_id']].add(gene)
    identity.append(dict(gene_id=gene, gene_name=symbol, hgnc_id=chosen.get('hgnc_id', ''),
                         HGNC_symbol=chosen.get('symbol', ''), mapping=method,
                         HGNC_ensembl_id=chosen.get('ensembl_gene_id', ''),
                         HGNC_locus_type=chosen.get('locus_type', ''),
                         histone=flag(bool(qualifying)), histone_families=';'.join(qualifying)))
write('gene_identity_mapping.tsv', identity)
write('histone_gene_mapping.tsv', [r for r in identity if r['histone']=='TRUE'])

print('Read the matched GO release and construct the two descendant sets.')
ontology_path = INVESTIGATION / 'annotations/go-basic_GO_release_2026-08-05.json.gz'
gaf_path = INVESTIGATION / 'annotations/HUMAN-uniprot_GO_release_2026-08-05.gaf.gz'
with gzip.open(ontology_path, 'rt') as handle:
    graph = json.load(handle)['graphs'][0]
assert '/2026-07-26/' in graph['meta']['version']
gid = lambda value: 'GO:'+m.group(1) if (m := re.search(r'GO[_:](\d{7})$', value or '')) else None
terms, alt = {}, {}
for node in graph['nodes']:
    term = gid(node['id'])
    if not term or node.get('type') != 'CLASS' or node.get('meta', {}).get('deprecated'): continue
    terms[term] = node.get('lbl', '')
    for item in node.get('meta', {}).get('basicPropertyValues', []):
        if item['pred'].endswith('#hasAlternativeId') and gid(item['val']):
            alt[gid(item['val'])] = term
children = defaultdict(set)
for edge in graph['edges']:
    child, parent = gid(edge['sub']), gid(edge['obj'])
    if child in terms and parent in terms and edge['pred'] in ('is_a', 'http://purl.obolibrary.org/obo/BFO_0000050'):
        children[parent].add(child)
roots = {'GO:0005886':'plasma_membrane', 'GO:0007165':'signal_transduction'}
descendants = {}
for root in roots:
    seen, stack = set(), [root]
    while stack:
        term = stack.pop()
        if term in seen: continue
        seen.add(term); stack.extend(children[term])
    descendants[root] = seen
assert terms['GO:0005886']=='plasma membrane' and terms['GO:0007165']=='signal transduction'
root_genes = {root:set() for root in roots}
evidence, header, mapping_methods = [], [], Counter()
with gzip.open(gaf_path, 'rt') as handle:
    for line in handle:
        if line.startswith('!'):
            header.append(line.rstrip()); continue
        row = line.rstrip('\n').split('\t')
        if len(row)<15 or row[0]!='UniProtKB' or row[12].split('|')[0]!='taxon:9606': continue
        if 'NOT' in row[3].split('|') or row[6]=='ND': continue
        term = alt.get(row[4], row[4])
        matching_roots = [root for root in roots if term in descendants[root]]
        if not matching_roots: continue
        product = re.sub(r'-\d+$', '', row[1])
        ids, method = protein_hgnc.get(product, set()), 'HGNC UniProt accession'
        if not ids and row[2] in by_symbol:
            ids, method = {by_symbol[row[2]]['hgnc_id']}, 'HGNC approved symbol'
        local_genes = set.union(*(hgnc_genes.get(hgnc_id,set()) for hgnc_id in ids)) if ids else set()
        if not ids and row[2] in gencode_symbols and gencode_symbols[row[2]] in genes:
            local_genes = {gencode_symbols[row[2]]}
            method = 'Unique full-reference GENCODE symbol'
        for gene in local_genes:
            for root in matching_roots:
                root_genes[root].add(gene)
                evidence.append(dict(gene_id=gene, gene_name=genes[gene], root=root,
                    annotation_go_id=term, product_id=row[1], symbol=row[2], evidence=row[6],
                    qualifier=row[3], reference=row[5], mapping=method))
                mapping_methods[method] += 1
assert any('2026-07-26' in line for line in header)
with gzip.open(OUT / 'membrane_signaling_annotation_evidence.tsv.gz', 'wt') as handle:
    writer = csv.DictWriter(handle, delimiter='\t', fieldnames=list(evidence[0]))
    writer.writeheader(); writer.writerows(evidence)
membrane = set.union(*root_genes.values())
write('membrane_signaling_gene_union.tsv', [dict(gene_id=gene,gene_name=genes[gene],
      plasma_membrane=flag(gene in root_genes['GO:0005886']),
      signal_transduction=flag(gene in root_genes['GO:0007165'])) for gene in sorted(membrane)])

print('Expand the chromatin union and add the membrane/signaling union.')
for row in membership:
    gene = row['gene_id']
    row['transcription_chromatin_previous'] = row['transcription_chromatin']
    row['histone'] = flag(gene in histones)
    row['transcription_chromatin'] = flag(row['transcription_chromatin_previous']=='TRUE' or gene in histones)
    row['plasma_membrane'] = flag(gene in root_genes['GO:0005886'])
    row['signal_transduction'] = flag(gene in root_genes['GO:0007165'])
    row['plasma_membrane_signaling'] = flag(gene in membrane)
write('novel_coding_gene_annotations.tsv', membership)
definitions = {r['category']:r for r in read(BASE / 'category_definitions.tsv')}
definitions['transcription_chromatin'] = dict(category='transcription_chromatin',
    label='Transcription and chromatin', definition='Union of the previous Human TFs v1.01 / EpiFactors 2.1 regulatory-protein set with HGNC H1, H2A, H2B, H3 and H4 histone families; HGNC snapshot 2026-10-15; deduplicate stable Ensembl gene IDs')
definitions['plasma_membrane_signaling'] = dict(category='plasma_membrane_signaling',
    label='Plasma membrane / signal transduction', definition='Union of positive human protein annotations to GO:0005886 or GO:0007165 and is_a/part_of descendants; matched official GO release 2026-08-05, ontology 2026-07-26; exclude NOT and ND, include IEA; deduplicate stable Ensembl gene IDs across roots')
write('category_definitions.tsv', [definitions[c] for c in CATEGORIES])
lists, counts, distinct, overlaps, coverage = [], [], [], [], []
for category in CATEGORIES:
    across = set()
    for assay in ASSAYS:
        subset = [r for r in membership if r['assay']==assay]
        selected = [r for r in subset if r[category]=='TRUE']
        ids = {r['gene_id'] for r in selected}; across.update(ids)
        assert len(ids)==len(selected)
        counts.append(dict(assay=assay,category=category,n_genes=len(ids),
                           total_novel_coding=len(subset),label=definitions[category]['label']))
        lists.extend(dict(assay=assay,category=category,gene_id=r['gene_id'],gene_name=r['gene_name']) for r in selected)
    distinct.append(dict(category=category,n_genes=len(across),label=definitions[category]['label']))
write('category_counts_by_assay.tsv', counts)
write('distinct_category_counts.tsv', distinct)
write('category_gene_lists.tsv', lists)
write('distinct_category_gene_lists.tsv', [dict(category=c,gene_id=g,gene_name=genes[g])
    for c in CATEGORIES for g in sorted({r['gene_id'] for r in lists if r['category']==c})])
for assay in ASSAYS+['All assays']:
    subset = [r for r in membership if r['assay']==assay] if assay in ASSAYS else membership
    total = {r['gene_id'] for r in subset}
    previous = {r['gene_id'] for r in subset if any(r[c]=='TRUE' for c in CATEGORIES if c not in ('transcription_chromatin','plasma_membrane_signaling')) or r['transcription_chromatin_previous']=='TRUE'}
    covered = {r['gene_id'] for r in subset if any(r[c]=='TRUE' for c in CATEGORIES)}
    coverage.append(dict(assay=assay,total_novel_coding=len(total),previous_six_category_union=len(previous),
                         new_seven_category_union=len(covered),newly_covered=len(covered-previous),
                         uncovered=len(total-covered),percent_covered=100*len(covered)/len(total)))
    p = {r['gene_id'] for r in subset if r['plasma_membrane']=='TRUE'}
    s = {r['gene_id'] for r in subset if r['signal_transduction']=='TRUE'}
    old_reg = {r['gene_id'] for r in subset if r['transcription_chromatin_previous']=='TRUE'}
    h = {r['gene_id'] for r in subset if r['histone']=='TRUE'}
    overlaps.append(dict(assay=assay,plasma_membrane=len(p),signal_transduction=len(s),
                         membrane_signaling_overlap=len(p&s),membrane_signaling_union=len(p|s),
                         previous_transcription_chromatin=len(old_reg),histones=len(h),
                         histones_already_in_regulation=len(h&old_reg),expanded_chromatin=len(h|old_reg)))
write('coverage_summary.tsv', coverage)
write('union_deduplication_audit.tsv', overlaps)
write('assay_totals.tsv',[dict(assay=a,total_novel_coding=sum(r['assay']==a for r in membership)) for a in ASSAYS])
new_counts = {r['category']:r['n_genes'] for r in distinct}
sentence = (f"Across QTL types, novel protein-coding targets included {new_counts['cancer_driver']:,} cancer drivers, "
    f"{new_counts['mendelian_disease']:,} Mendelian disease genes, {new_counts['lof_constrained']:,} genes under strong loss-of-function constraint, "
    f"{new_counts['transcription_chromatin']:,} transcription and chromatin genes (including histones), "
    f"{new_counts['rna_processing_splicing']:,} genes involved in RNA processing and splicing, "
    f"{new_counts['plasma_membrane_signaling']:,} genes annotated to plasma membrane or signal transduction, "
    f"and {new_counts['drug_target']:,} drug targets (Fig. 1C).")
(OUT / 'manuscript_sentence.txt').write_text(sentence+'\n\nEach gene is counted once per category across eQTLs, sQTLs and pQTLs. Categories can overlap. The drug-target category includes approved or clinical-stage drugs.\n')
audit = dict(previous_data=str(BASE),total_distinct_coding_genes=4146,categories=CATEGORIES,
    distinct_counts=new_counts,coverage=coverage,union_audit=overlaps,GO_roots=roots,
    GAF_header=header,GO_release='2026-08-05',ontology_version=graph['meta']['version'],
    annotation_mapping='Full HGNC protein accession mapping; unique approved HGNC symbol then unique full-reference GENCODE symbol fallback; exact Ensembl identity preferred for local genes',
    evidence_rows=len(evidence),mapping_methods=dict(mapping_methods),
    histone_families=['H1 histones','H2A histones','H2B histones','H3 histones','H4 histones'],
    stable_gene_ID_deduplication=True,previous_draft_preserved=True,
    input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                  [BASE/'novel_coding_gene_annotations.tsv',hgnc_path,gencode_path,ontology_path,gaf_path]})
(OUT / 'data_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(dict(distinct_counts=new_counts,per_assay_union_audit=overlaps,coverage=coverage),indent=2))
