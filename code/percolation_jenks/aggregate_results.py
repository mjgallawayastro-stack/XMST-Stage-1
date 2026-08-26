#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import json, math, zipfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent
ROOT=BASE/'output'
RUNS=ROOT/'runs'
N=1000

def wilson(k,n,z=1.959963984540054):
    if n==0:return (np.nan,np.nan)
    p=k/n; den=1+z*z/n
    ctr=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return ctr-half,ctr+half

summary_files=sorted(RUNS.glob('mc_*/run_summary.json'))
if len(summary_files)!=N:
    raise RuntimeError(f'Expected {N} completed runs, found {len(summary_files)}')
S=pd.DataFrame([json.loads(p.read_text()) for p in summary_files]).sort_values('realisation')
S.to_csv(ROOT/'mc1000_realisation_summary.csv',index=False)
rec=[]; splits=[]; place=[]
for d in sorted(RUNS.glob('mc_*')):
    rp=d/'template_recovery.csv'
    if rp.exists(): rec.append(pd.read_csv(rp))
    sp=d/'accepted_splits.csv'
    if sp.exists() and sp.stat().st_size>1:
        try:
            x=pd.read_csv(sp)
            if len(x): splits.append(x)
        except pd.errors.EmptyDataError: pass
    pp=d/'placements.csv'
    if pp.exists(): place.append(pd.read_csv(pp))
R=pd.concat(rec,ignore_index=True)
P=pd.concat(place,ignore_index=True)
B=pd.concat(splits,ignore_index=True) if splits else pd.DataFrame()
R.to_csv(ROOT/'mc1000_template_recovery.csv',index=False)
P.to_csv(ROOT/'mc1000_placements.csv',index=False)
B.to_csv(ROOT/'mc1000_accepted_branch_splits.csv',index=False)
if len(B):
    B['split_class']=np.select([
        B['correct_injected_separation'].fillna(False) & (B['dominant_fraction_a']>=0.5) & (B['dominant_fraction_b']>=0.5),
        B['correct_injected_separation'].fillna(False),
        (B['dominant_truth_a']=='background') | (B['dominant_truth_b']=='background'),
        B['dominant_truth_a']==B['dominant_truth_b'],
    ],['clean_injected_separation','mixed_injected_separation','background_involved','same_template_fragmentation'],default='other')
    B.to_csv(ROOT/'mc1000_accepted_branch_splits.csv',index=False)
T=(R.groupby(['stage','template_label'],as_index=False)
   .agg(true_size=('true_size','first'),detection_rate=('detected_ge50pct','mean'),
        strict_rate=('strict_recovery_c80_p50','mean'),completeness_mean=('completeness','mean'),
        completeness_median=('completeness','median'),purity_mean=('purity','mean'),
        purity_median=('purity','median'),fragment_count_mean=('fragment_count','mean'),
        merger_rate=('n_merged_templates',lambda x: np.mean(np.asarray(x)>0))))
T.to_csv(ROOT/'mc1000_template_summary.csv',index=False)
metrics=[]
def add(name,value,den=None):
    row={'metric':name,'value':value}
    if den is not None:
        lo,hi=wilson(int(value),int(den)); row.update({'denominator':den,'rate':value/den,'wilson95_low':lo,'wilson95_high':hi})
    metrics.append(row)
for prefix in ('global','final'):
    add(f'{prefix}_all_eight_detected_count',int(S[f'{prefix}_all_eight_detected'].sum()),N)
    add(f'{prefix}_all_eight_strict_count',int(S[f'{prefix}_all_eight_strict'].sum()),N)
    add(f'{prefix}_overlap_pair_resolved_count',int(S[f'{prefix}_overlap_pair_resolved'].sum()),N)
merged=~S.global_overlap_pair_resolved.astype(bool)
resolved=merged & S.final_overlap_pair_resolved.astype(bool)
add('overlap_resolved_by_pruning_given_global_merge_count',int(resolved.sum()),int(merged.sum()))
add('accepted_split_realisations',int((S.accepted_branch_splits>0).sum()),N)
add('no_split_realisations',int((S.accepted_branch_splits==0).sum()),N)
add('realisations_with_noninjected_or_mixed_split',int((S.noninjected_or_mixed_branch_splits>0).sum()),N)
if len(B):
    for cls,c in B.split_class.value_counts().items(): add(f'accepted_split_class_{cls}',int(c),len(B))
pd.DataFrame(metrics).to_csv(ROOT/'mc1000_headline_metrics.csv',index=False)
cols=['percolation_operational_cut_pc','jenks_threshold_pc','jenks_bootstrap_sd_pc','global_groups_ge10','final_groups_ge10','accepted_branch_splits','global_mean_completeness','final_mean_completeness','global_mean_purity','final_mean_purity']
Q=[]
for c in cols:
    x=pd.to_numeric(S[c],errors='coerce').dropna()
    Q.append({'metric':c,'mean':x.mean(),'sd':x.std(ddof=1),'p2_5':x.quantile(.025),'p16':x.quantile(.16),'median':x.median(),'p84':x.quantile(.84),'p97_5':x.quantile(.975),'min':x.min(),'max':x.max()})
pd.DataFrame(Q).to_csv(ROOT/'mc1000_distribution_summary.csv',index=False)
pd.crosstab(S.global_overlap_pair_resolved,S.final_overlap_pair_resolved,margins=True).to_csv(ROOT/'mc1000_overlap_transition_table.csv')

# Publication/support plots
fig,ax=plt.subplots(figsize=(8,5)); ax.hist(S.jenks_threshold_pc,bins=30); ax.axvline(S.jenks_threshold_pc.median(),linestyle='--',label=f"Median {S.jenks_threshold_pc.median():.3f} pc"); ax.set(xlabel='Jenks fracture threshold (pc)',ylabel='Realisations',title='XMST fracture-scale distribution: 1,000 realisations'); ax.legend(); fig.tight_layout(); fig.savefig(ROOT/'mc1000_jenks_threshold_distribution.png',dpi=220); plt.close(fig)
fig,ax=plt.subplots(figsize=(7,6)); ax.scatter(S.global_mean_completeness,S.global_mean_purity,s=10,label='Global Percolation–Jenks'); ax.scatter(S.final_mean_completeness,S.final_mean_purity,s=10,label='After reddening pruning'); ax.set(xlabel='Mean template completeness',ylabel='Mean template purity',title='Purity–completeness across 1,000 realisations'); ax.legend(); fig.tight_layout(); fig.savefig(ROOT/'mc1000_purity_completeness.png',dpi=220); plt.close(fig)
# convergence
for col,label,fname in [('jenks_threshold_pc','Mean Jenks fracture scale (pc)','mc1000_convergence_jenks.png'),('global_mean_completeness','Mean global completeness','mc1000_convergence_completeness.png'),('global_mean_purity','Mean global purity','mc1000_convergence_purity.png')]:
    x=pd.to_numeric(S[col],errors='coerce').to_numpy(float); cm=np.cumsum(x)/np.arange(1,len(x)+1)
    fig,ax=plt.subplots(figsize=(8,5)); ax.plot(np.arange(1,len(x)+1),cm); ax.set(xlabel='Number of Monte Carlo realisations',ylabel=label,title=f'Convergence of {label.lower()}'); fig.tight_layout(); fig.savefig(ROOT/fname,dpi=220); plt.close(fig)

n_global=int(S.global_all_eight_detected.sum()); n_final=int(S.final_all_eight_detected.sum())
glob_res=int(S.global_overlap_pair_resolved.sum()); fin_res=int(S.final_overlap_pair_resolved.sum())
report=f'''# XMST 1,000-realisation controlled Monte Carlo report

## Headline results

- All eight injected templates recovered globally: **{n_global}/{N}**
- All eight retained after reddening pruning: **{n_final}/{N}**
- Mean global completeness: **{S.global_mean_completeness.mean():.4f}**
- Mean global purity: **{S.global_mean_purity.mean():.4f}**
- Mean final completeness: **{S.final_mean_completeness.mean():.4f}**
- Mean final purity: **{S.final_mean_purity.mean():.4f}**
- Global overlap resolution: **{glob_res}/{N}**
- Final overlap resolution: **{fin_res}/{N}**
- Median Jenks fracture scale: **{S.jenks_threshold_pc.median():.3f} pc**
- Between-realisation Jenks scatter: **{S.jenks_threshold_pc.std(ddof=1):.3f} pc**
- Median within-realisation bootstrap uncertainty: **{S.jenks_bootstrap_sd_pc.median():.3f} pc**

## Runtime

- Median field runtime: **{S.elapsed_seconds.median():.1f} s**
- Total summed CPU field time: **{S.elapsed_seconds.sum()/3600:.2f} h**

The CSV tables in this directory contain the full per-field and per-template results. The convergence plots show whether the original 100-field validation had already stabilised.
'''
(ROOT/'MC1000_REPORT.md').write_text(report,encoding='utf-8')

package=ROOT/'XMST_MC1000_FINAL_RESULTS.zip'
include=[p for p in ROOT.iterdir() if p.is_file() and p.name != package.name]
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for p in include:z.write(p,arcname=p.name)
print(report)
print(f'Final package: {package}')
