"""Figure 3 inputs, step 2: species-ARG Spearman correlations (Figure 3C).

ARGs detected in >= 2 libraries are correlated with the relative abundance of each prioritized species
across all analysed libraries (two-sided Spearman test); Benjamini-Hochberg q values across all tests are
added. Run after 01_arg_quantification.py.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.dirname(HERE)]          # this folder and the repository root
from config import OUTPUT_DIR, in_file, out_file  # noqa: E402

MOD = 'figure3'

dm_tab = pd.read_csv(os.path.join(OUTPUT_DIR, MOD, 'arg_dpm_matrix.tsv'), sep='\t')
libs = [c for c in dm_tab.columns if c != 'feature_id']
dm = dm_tab.set_index('feature_id')
meta = pd.read_csv(os.path.join(OUTPUT_DIR, MOD, 'arg_feature_summary.tsv'), sep='\t').set_index('feature_id')

S = pd.read_csv(in_file('species_relative_abundance.tsv'), sep='\t', index_col=0)
S.index = [str(i).split('s__')[-1].replace('_', ' ') for i in S.index]
S = S.reindex(columns=libs)
sp_names = sorted(S.index.unique())

nlib = (dm > 0).sum(axis=1)
ge2 = [f for f in dm.index if nlib[f] >= 2]
rows = []
for s in sp_names:
    w = S.loc[s].astype(float).values
    smax = np.nanmax(w) if len(w) else np.nan
    for f in ge2:
        v = dm.loc[f].values.astype(float)
        if np.std(v) == 0 or np.std(w) == 0:
            rho, pv = np.nan, np.nan
        else:
            rho, pv = spearmanr(v, w)
        rows.append((s, f, meta.loc[f, 'gene_name'], meta.loc[f, 'Primary_Drug_Class'],
                     meta.loc[f, 'CARD_Drug_Class'], rho, pv, len(libs), smax,
                     v.max(), int((v > 0).sum())))
sp_out = pd.DataFrame(rows, columns=['species', 'ARG_feature_id', 'ARG_gene', 'Primary_Drug_Class',
                                     'CARD_Drug_Class', 'spearman_rho', 'p_value', 'n_libraries',
                                     'species_max_RA_pct', 'ARG_max_DPM', 'ARG_n_passed_libraries'])
m = sp_out['p_value'].notna()
q = np.full(len(sp_out), np.nan)
pv = sp_out.loc[m, 'p_value'].values
o = np.argsort(pv)
n = len(pv)
adj = np.minimum.accumulate((pv[o] * n / (np.arange(n) + 1))[::-1])[::-1]
qq = np.empty(n)
qq[o] = np.minimum(adj, 1)
q[np.where(m)[0]] = qq
sp_out['q_value_BH_all_tests'] = q
sp_out.to_csv(out_file(MOD, 'arg_species_spearman.tsv'), sep='\t', index=False)
print(f'species={len(sp_names)}  ARGs in >=2 libraries={len(ge2)}  tests={len(sp_out)}')
