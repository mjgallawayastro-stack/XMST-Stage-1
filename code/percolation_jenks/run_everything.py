#!/usr/bin/env python3
from __future__ import annotations
import argparse, os, subprocess, sys, time
from pathlib import Path

BASE=Path(__file__).resolve().parent
OUT=BASE/'output'; RUNS=OUT/'runs'
N=1000

def completed():
    return len(list(RUNS.glob('mc_*/run_summary.json')))

def main():
    ap=argparse.ArgumentParser(description='Run and aggregate the complete 1,000-field XMST controlled Monte Carlo campaign.')
    ap.add_argument('--workers',type=int,default=max(1,min(8,(os.cpu_count() or 4)//2)),help='Parallel worker processes (default: half CPU count, maximum 8)')
    ap.add_argument('--overwrite',action='store_true',help='Rerun completed fields')
    args=ap.parse_args()
    workers=max(1,min(args.workers,N))
    OUT.mkdir(exist_ok=True); RUNS.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env[k]='1'
    print(f'XMST 1,000-field production run')
    print(f'Package: {BASE}')
    print(f'Workers: {workers}')
    print(f'Already complete: {completed()}/{N}')
    if completed()<N or args.overwrite:
        bounds=[]
        chunk=(N+workers-1)//workers
        for i in range(workers):
            a=i*chunk; b=min(N,(i+1)*chunk)
            if a<b: bounds.append((a,b))
        procs=[]; logs=[]
        for i,(a,b) in enumerate(bounds,1):
            log=(OUT/f'worker_{i:02d}_{a:04d}_{b:04d}.log').open('a',buffering=1)
            cmd=[sys.executable,str(BASE/'xmst_worker.py'),'--start',str(a),'--end',str(b)]
            if args.overwrite: cmd.append('--overwrite')
            procs.append(subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env,cwd=BASE))
            logs.append(log)
        start=time.time(); last=-1
        try:
            while any(p.poll() is None for p in procs):
                n=completed()
                if n!=last:
                    elapsed=time.time()-start
                    rate=max(0,(n-last if last>=0 else 0))
                    print(f'Completed: {n}/{N} ({100*n/N:.1f}%)',flush=True)
                    last=n
                time.sleep(5)
        except KeyboardInterrupt:
            print('\nStopping workers. Completed fields are checkpointed; rerun the same command to resume.')
            for p in procs:p.terminate()
            for p in procs:p.wait()
            return 130
        finally:
            for f in logs:f.close()
        failed=[p.returncode for p in procs if p.returncode]
        if failed:
            print(f'One or more workers failed: {failed}. Check output/worker_*.log',file=sys.stderr)
            return 1
    n=completed()
    if n!=N:
        print(f'Run incomplete: {n}/{N}. Rerun this command to resume.',file=sys.stderr)
        return 1
    print('All 1,000 fields complete. Aggregating results...')
    subprocess.run([sys.executable,str(BASE/'aggregate_results.py')],check=True,cwd=BASE,env=env)
    print('\nDONE')
    print(f'Report: {OUT / "MC1000_REPORT.md"}')
    print(f'Final results ZIP: {OUT / "XMST_MC1000_FINAL_RESULTS.zip"}')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
