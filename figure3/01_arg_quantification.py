"""Figure 3 inputs, step 1: ARG quantification.

Retains ARGs with coverage breadth >= 50% and either DPM > 1 or read count > 10 in a library, and
writes the per-library DPM matrix, the per-ARG annotation table read by figure3.R and the 35
high-frequency ARGs (prevalence >= 50% and mean DPM >= 3) in the display order of Figures 3B-3C.
"""
import os

import pandas as pd

DATA = 'data'                              # inputs (run from the repository root)
OUT = os.path.join('output', 'figure3')    # outputs
os.makedirs(OUT, exist_ok=True)

BREADTH_MIN = 50    # coverage breadth: % of reference protein residues covered by >= 1 read
DPM_MIN = 1         # an ARG is retained in a library if DPM > 1 ...
READS_MIN = 10      # ... or read count > 10

# Display order of primary drug classes in Figure 3B
CLASS_ORDER = ["fluoroquinolone", "aminoglycoside", "multiple_drug_classes",
               "multidrug_efflux_or_multiclass", "tetracycline", "peptide/polymyxin",
               "aminocoumarin", "nitroimidazole"]

libs = [str(s) for s in pd.read_csv(os.path.join(DATA, 'library_metadata.csv'))['sample']]
N = len(libs)
a = pd.read_csv(os.path.join(DATA, 'arg_read_level_metrics.tsv'), sep='\t')
a = a[a['sample'].astype(str).isin(libs)].copy()
p = a[(a['coverage_breadth_pct'] >= BREADTH_MIN) &
      ((a['depth_per_million'] > DPM_MIN) | (a['read_count'] > READS_MIN))].copy()
ann = pd.read_csv(os.path.join(DATA, 'arg_feature_annotation.tsv'), sep='\t').set_index('feature_id')

# ---------- per-library DPM matrix (ARGs x libraries; rows and columns sorted by identifier) ----------
keep = sorted(p['feature_id'].unique())
dm = (p.pivot_table(index='feature_id', columns='sample', values='depth_per_million',
                    aggfunc='sum', fill_value=0)
        .reindex(index=keep, columns=sorted(libs), fill_value=0))
dm.reset_index().to_csv(os.path.join(OUT, 'arg_dpm_matrix.tsv'), sep='\t', index=False)

# ---------- annotation of the retained ARGs ----------
fs = ann.reset_index()
fs = fs[fs['feature_id'].isin(keep)].copy()
fs['n_libraries_detected'] = [(dm.loc[f] > 0).sum() for f in fs['feature_id']]
fs['prevalence_pct'] = fs['n_libraries_detected'] / N * 100
fs['mean_DPM'] = [dm.loc[f].mean() for f in fs['feature_id']]
fs = fs.sort_values(['n_libraries_detected', 'mean_DPM'], ascending=False)
fs[['feature_id', 'gene_name', 'ARO_ID', 'Primary_Drug_Class', 'CARD_Drug_Class', 'CARD_Resistance_Mechanism',
    'CARD_AMR_Gene_Family', 'n_libraries_detected', 'prevalence_pct', 'mean_DPM']].to_csv(
    os.path.join(OUT, 'arg_feature_summary.tsv'), sep='\t', index=False)

# ---------- 35 high-frequency ARGs in display order ----------
cand = fs[(fs['prevalence_pct'] >= 50) & (fs['mean_DPM'] >= 3)].copy()
cand = cand.sort_values(['n_libraries_detected', 'mean_DPM'], ascending=False).head(35).copy()
cand['class_rank'] = cand['Primary_Drug_Class'].apply(
    lambda c: CLASS_ORDER.index(c) if c in CLASS_ORDER else len(CLASS_ORDER))
cand = cand.sort_values(['class_rank', 'mean_DPM', 'n_libraries_detected'],
                        ascending=[True, False, False])
cand['plot_order'] = range(1, len(cand) + 1)
top = pd.DataFrame({
    'plot_order': cand['plot_order'].values, 'feature_id': cand['feature_id'].values,
    'gene_name': cand['gene_name'].values, 'ARO_ID': cand['ARO_ID'].values,
    'Primary_Drug_Class': cand['Primary_Drug_Class'].values,
    'CARD_Drug_Class': cand['CARD_Drug_Class'].values,
    'CARD_Resistance_Mechanism': cand['CARD_Resistance_Mechanism'].values,
    'CARD_AMR_Gene_Family': cand['CARD_AMR_Gene_Family'].values,
    'n_passed_samples': cand['n_libraries_detected'].values,
    'prevalence_pct': cand['prevalence_pct'].values,
    'mean_DPM_all': cand['mean_DPM'].values})
top.to_csv(os.path.join(OUT, 'arg_top35.csv'), index=False)
print(f'libraries={N}  retained ARGs={len(keep)}  high-frequency ARGs shown={len(top)}')
