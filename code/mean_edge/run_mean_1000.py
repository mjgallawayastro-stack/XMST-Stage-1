#!/usr/bin/env python3
from __future__ import annotations
import argparse, os, subprocess, sys, time
from pathlib import Path

BASE=Path(__file__).resolve().parent
OUT=BASE/'output'; RUNS=OUT/'runs'; N=1000

def completed():
    return len(list(RUNS.glob('mc_*/run_summary.json')))

def check_dependencies():
    missing=[]
    for mod in ('numpy','pandas','scipy','sklearn','matplotlib'):
        try: __import__(mod)
        except ImportError: missing.append(mod)
    if missing:
        print('Missing Python packages: '+', '.join(missing), file=sys.stderr)
        print('Install with:', file=sys.stderr)
        print('  python3 -m pip install numpy pandas scipy scikit-learn matplotlib', file=sys.stderr)
        raise SystemExit(2)

def main():
    check_dependencies()
    ap=argparse.ArgumentParser(description='Run mean-edge fracturing on the same 1,000 controlled Monte Carlo fields.')
    ap.add_argument('--workers',type=int,default=max(1,min(8,(os.cpu_count() or 4)//2)))
    ap.add_argument('--overwrite',action='store_true')
    args=ap.parse_args()
    workers=max(1,min(args.workers,N))
    OUT.mkdir(exist_ok=True); RUNS.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env[k]='1'
    print('Mean-edge 1,000-field controlled Monte Carlo comparison')
    print('Fracture rule: arithmetic mean of the full global MST edge-length distribution')
    print(f'Workers: {workers}')
    print(f'Already complete: {completed()}/{N}')
    if completed()<N or args.overwrite:
        chunk=(N+workers-1)//workers
        procs=[]; logs=[]
        for i in range(workers):
            a=i*chunk; b=min(N,(i+1)*chunk)
            if a>=b: continue
            log=(OUT/f'worker_{i+1:02d}_{a:04d}_{b:04d}.log').open('a',buffering=1)
            cmd=[sys.executable,str(BASE/'mean_worker.py'),'--start',str(a),'--end',str(b)]
            if args.overwrite: cmd.append('--overwrite')
            procs.append(subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env,cwd=BASE))
            logs.append(log)
        last=-1
        try:
            while any(p.poll() is None for p in procs):
                n=completed()
                if n!=last:
                    print(f'Completed: {n}/{N} ({100*n/N:.1f}%)',flush=True)
                    last=n
                time.sleep(5)
        except KeyboardInterrupt:
            print('\nStopped safely. Rerun the same command to resume.')
            for p in procs:p.terminate()
            for p in procs:p.wait()
            return 130
        finally:
            for f in logs:f.close()
        if any(p.returncode for p in procs):
            print('One or more workers failed. Check output/worker_*.log',file=sys.stderr)
            return 1
    if completed()!=N:
        print(f'Incomplete: {completed()}/{N}. Rerun to resume.',file=sys.stderr)
        return 1
    subprocess.run([sys.executable,str(BASE/'aggregate_mean_results.py')],check=True,cwd=BASE,env=env)
    print('\nDONE')
    print(f'Report: {OUT/"MEAN_MC1000_REPORT.md"}')
    print(f'Results ZIP: {OUT/"MEAN_MC1000_FINAL_RESULTS.zip"}')
    return 0

if __name__=='__main__': raise SystemExit(main())
