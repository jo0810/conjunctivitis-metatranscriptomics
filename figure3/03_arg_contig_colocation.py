"""Figure 3 inputs, step 3: contig-level ARG-host co-location (outlined cells in Figure 3C).

Uses high-confidence ARG-host assignments on assembled contigs (ARG open reading frame called by RGI
and host sequence on the same contig; host confirmed against NCBI nr and nt). Assignments are matched
to the high-frequency ARGs shown in Figure 3B and summarized per species-ARG pair.
Run after 01_arg_quantification.py.
"""
import os
import re

import pandas as pd

DATA = 'data'                              # inputs (run from the repository root)
OUT = os.path.join('output', 'figure3')    # outputs
os.makedirs(OUT, exist_ok=True)

h = pd.read_csv(os.path.join(DATA, 'arg_contig_host_assignments.tsv'), sep='\t', dtype=str).fillna('')
h = h[~h['host_taxon'].str.contains('synthetic', case=False, na=False)]
top = pd.read_csv(os.path.join(OUT, 'arg_top35.csv'))


def norm(s):
    """CARD best-hit name -> bare gene name."""
    s = str(s).strip()
    s = re.sub(r'\s+with\s+mutation[s]?\s+conferring.*$', '', s, flags=re.I)
    s = re.sub(r'\s+conferring\s+resistance.*$', '', s, flags=re.I)
    s = re.sub(r'^[A-Z][a-z]+\s+[a-z]+\s+', '', s)
    return s.strip().lower()


genes = {g.lower(): g for g in top['gene_name']}


def match_top(g):
    if g in genes:
        return genes[g]
    for k, v in genes.items():                 # species-prefixed names such as Ecol_acrA
        if k.split('_')[-1] == g:
            return v
    return None


h['g'] = h['best_hit_aro'].apply(norm)
h['gene_name'] = h['g'].apply(match_top)
hit = h[h['gene_name'].notna()]
pairs = (hit.groupby(['host_taxon', 'gene_name'])
            .agg(n_contigs=('contig', 'nunique'), n_libraries=('sample', 'nunique'))
            .reset_index().rename(columns={'host_taxon': 'species'}))
pairs.to_csv(os.path.join(OUT, 'arg_contig_colocated_pairs.csv'), index=False)

print(f'species-ARG pairs with contig support: {len(pairs)}')
