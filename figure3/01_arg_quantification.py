"""Figure 3 inputs, step 1: ARG quantification.

Retains ARGs with coverage breadth >= 50% and either DPM > 1 or read count > 10 in a library, and
writes the per-library DPM matrix, the per-ARG annotation table read by figure3.R and the 35
high-frequency ARGs (prevalence >= 50% and mean DPM >= 3) in the display order of Figures 3B-3C.
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.dirname(HERE)]          # this folder and the repository root
from config import out_file  # noqa: E402
from arg_common import annotation, filtered_alignments, libraries  # noqa: E402

MOD = 'figure3'

# Display order of primary drug classes in Figure 3B
CLASS_ORDER = ["fluoroquinolone", "aminoglycoside", "multiple_drug_classes",
               "multidrug_efflux_or_multiclass", "tetracycline", "peptide/polymyxin",
               "aminocoumarin", "nitroimidazole"]

libs, _ = libraries()
N = len(libs)
p = filtered_alignments(libs)
ann = annotation()

# ---------- per-library DPM matrix (ARGs x libraries; rows and columns sorted by identifier) ----------
keep = sorted(p['feature_id'].unique())
dm = (p.pivot_table(index='feature_id', columns='sample', values='depth_per_million',
                    aggfunc='sum', fill_value=0)
        .reindex(index=keep, columns=sorted(libs), fill_value=0))
dm.reset_index().to_csv(out_file(MOD, 'arg_dpm_matrix.tsv'), sep='\t', index=False)

# ---------- annotation of the retained ARGs ----------
fs = ann.reset_index()
fs = fs[fs['feature_id'].isin(keep)].copy()
fs['n_libraries_detected'] = [(dm.loc[f] > 0).sum() for f in fs['feature_id']]
fs['prevalence_pct'] = fs['n_libraries_detected'] / N * 100
fs['mean_DPM'] = [dm.loc[f].mean() for f in fs['feature_id']]
fs = fs.sort_values(['n_libraries_detected', 'mean_DPM'], ascending=False)
fs[['feature_id', 'gene_name', 'ARO_ID', 'Primary_Drug_Class', 'CARD_Drug_Class', 'CARD_Resistance_Mechanism',
    'CARD_AMR_Gene_Family', 'n_libraries_detected', 'prevalence_pct', 'mean_DPM']].to_csv(
    out_file(MOD, 'arg_feature_summary.tsv'), sep='\t', index=False)

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
top.to_csv(out_file(MOD, 'arg_top35.csv'), index=False)
print(f'libraries={N}  retained ARGs={len(keep)}  high-frequency ARGs shown={len(top)}')
