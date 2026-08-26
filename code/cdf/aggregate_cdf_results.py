#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import json, math, zipfile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent; ROOT=BASE/'output'; RUNS=ROOT/'runs'; N=1000

def wilson(k,n,z=1.959963984540054):
    if n==0:return (np.nan,np.nan)
    p=k/n; den=1+z*z/n
    ctr=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return ctr-half,ctr+half

files=sorted(RUNS.glob('mc_*/run_summary.json'))
if len(files)!=N: raise RuntimeError(f'Expected {N} runs, found {len(files)}')
S=pd.DataFrame([json.loads(p.read_text()) for p in files]).sort_values('realisation')
S.to_csv(ROOT/'cdf_mc1000_realisation_summary.csv',index=False)
rec=[]; place=[]
for d in sorted(RUNS.glob('mc_*')):
    rec.append(pd.read_csv(d/'template_recovery.csv'))
    place.append(pd.read_csv(d/'placements.csv'))
R=pd.concat(rec,ignore_index=True); P=pd.concat(place,ignore_index=True)
R.to_csv(ROOT/'cdf_mc1000_template_recovery.csv',index=False)
P.to_csv(ROOT/'cdf_mc1000_placements.csv',index=False)
T=(R.groupby('template_label',as_index=False)
   .agg(true_size=('true_size','first'),detection_rate=('detected_ge50pct','mean'),
        strict_rate=('strict_recovery_c80_p50','mean'),completeness_mean=('completeness','mean'),
        completeness_median=('completeness','median'),purity_mean=('purity','mean'),
        purity_median=('purity','median'),fragment_count_mean=('fragment_count','mean'),
        merger_rate=('n_merged_templates',lambda x: np.mean(np.asarray(x)>0))))
T.to_csv(ROOT/'cdf_mc1000_template_summary.csv',index=False)
metrics=[]
def add(name,k,n=None):
    row={'metric':name,'value':k}
    if n is not None:
        lo,hi=wilson(int(k),int(n)); row.update({'denominator':n,'rate':k/n,'wilson95_low':lo,'wilson95_high':hi})
    metrics.append(row)
add('all_eight_detected_count',int(S.global_all_eight_detected.sum()),N)
add('all_eight_strict_count',int(S.global_all_eight_strict.sum()),N)
add('overlap_pair_resolved_count',int(S.global_overlap_pair_resolved.sum()),N)
pd.DataFrame(metrics).to_csv(ROOT/'cdf_mc1000_headline_metrics.csv',index=False)
cols=['cdf_threshold_pc','cdf_bootstrap_sd_pc','cdf_fit_improvement_fraction','global_groups_ge10','global_noninjected_dominant_groups','global_mean_completeness','global_mean_purity']
rows=[]
for c in cols:
    x=pd.to_numeric(S[c],errors='coerce').dropna()
    rows.append({'metric':c,'mean':x.mean(),'sd':x.std(ddof=1),'p2_5':x.quantile(.025),'p16':x.quantile(.16),'median':x.median(),'p84':x.quantile(.84),'p97_5':x.quantile(.975),'min':x.min(),'max':x.max()})
pd.DataFrame(rows).to_csv(ROOT/'cdf_mc1000_distribution_summary.csv',index=False)

fig,ax=plt.subplots(figsize=(8,5)); ax.hist(S.cdf_threshold_pc,bins=30); ax.axvline(S.cdf_threshold_pc.median(),linestyle='--',label=f"Median {S.cdf_threshold_pc.median():.3f} pc"); ax.set(xlabel='CDF fracture threshold (pc)',ylabel='Realisations',title='CDF fracture-scale distribution: 1,000 realisations'); ax.legend(); fig.tight_layout(); fig.savefig(ROOT/'cdf_mc1000_threshold_distribution.png',dpi=220); plt.close(fig)
for col,label,fname in [('cdf_threshold_pc','Mean CDF fracture scale (pc)','cdf_mc1000_convergence_threshold.png'),('global_mean_completeness','Mean global completeness','cdf_mc1000_convergence_completeness.png'),('global_mean_purity','Mean global purity','cdf_mc1000_convergence_purity.png')]:
    x=pd.to_numeric(S[col],errors='coerce').to_numpy(float); cm=np.cumsum(x)/np.arange(1,len(x)+1)
    fig,ax=plt.subplots(figsize=(8,5)); ax.plot(np.arange(1,len(x)+1),cm); ax.set(xlabel='Number of Monte Carlo realisations',ylabel=label,title=f'Convergence of {label.lower()}'); fig.tight_layout(); fig.savefig(ROOT/fname,dpi=220); plt.close(fig)
fig,ax=plt.subplots(figsize=(7,6)); ax.scatter(S.global_mean_completeness,S.global_mean_purity,s=10); ax.set(xlabel='Mean template completeness',ylabel='Mean template purity',title='CDF purity–completeness across 1,000 realisations'); fig.tight_layout(); fig.savefig(ROOT/'cdf_mc1000_purity_completeness.png',dpi=220); plt.close(fig)

report=f'''# CDF fracture comparison: 1,000 controlled Monte Carlo fields

The same deterministic 1,000 synthetic fields and injected templates used in the Percolation–Jenks campaign were regenerated from the same seeds. The only methodological change was the global fracture rule: a two-line least-squares breakpoint fitted to the full empirical MST edge-length CDF.

## Headline results

- All eight injected templates detected: **{int(S.global_all_eight_detected.sum())}/{N}**
- All eight meeting strict recovery criterion: **{int(S.global_all_eight_strict.sum())}/{N}**
- Mean completeness: **{S.global_mean_completeness.mean():.4f}**
- Mean purity: **{S.global_mean_purity.mean():.4f}**
- Overlap pair resolved: **{int(S.global_overlap_pair_resolved.sum())}/{N}**
- Median CDF fracture scale: **{S.cdf_threshold_pc.median():.3f} pc**
- Between-realisation threshold scatter: **{S.cdf_threshold_pc.std(ddof=1):.3f} pc**
- Median within-realisation bootstrap uncertainty: **{S.cdf_bootstrap_sd_pc.median():.3f} pc**
- Median number of retained groups (N >= 10): **{S.global_groups_ge10.median():.0f}**
- Median field runtime: **{S.elapsed_seconds.median():.1f} s**

The generated CSV files are comparison-ready against the Percolation–Jenks 1,000-field outputs.
'''
(ROOT/'CDF_MC1000_REPORT.md').write_text(report)
package=ROOT/'CDF_MC1000_FINAL_RESULTS.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for p in ROOT.iterdir():
        if p.is_file() and p.name!=package.name:z.write(p,arcname=p.name)
print(report)
