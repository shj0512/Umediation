"""Reproduce Umediation's simulation study; see validation/README.md."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
import pandas as pd
from scipy.stats import norm
from threadpoolctl import threadpool_limits
from inference import analyze_one, generate

HERE=Path(__file__).resolve().parent
PHASES=('core','screen','extensions','precision')

def run_job(task,design,output):
    st=next(s for s in design['strata'] if s['name']==task['stratum'])
    cases=[c for c in design['cases'] if c['stratum']==st['name'] and
           (not task.get('case_names') or c['case'] in task['case_names'])]
    folder=Path(output)/task['phase']
    folder.mkdir(parents=True,exist_ok=True)
    target=folder/f"{st['name']}_{task['start']}_{task['stop']}_b{task['B']}.csv.gz"
    stamp=target.with_suffix('').with_suffix('.complete.json')
    task_hash=hashlib.sha256(json.dumps(task,sort_keys=True).encode()).hexdigest()
    if target.exists() and stamp.exists():
        saved=json.loads(stamp.read_text())
        if saved.get('task_sha256')==task_hash and saved.get('result_sha256')==hashlib.sha256(target.read_bytes()).hexdigest():
            return str(target)
        raise RuntimeError(f'Output does not match the requested task: {target}')
    rows=[]
    start=time.perf_counter()
    with threadpool_limits(limits=1):
        for rep in range(task['start'],task['stop']+1):
            data=generate(st,rep)
            data['rep']=rep
            result=analyze_one(data,st,cases,task['B'],task['wild'],task['pairs'])
            expected=len(cases)*(1+int(task['pairs'])+2*int(task['wild']))
            if len(result)!=expected or any(row.get('outer_ok')!=1 for row in result):
                raise RuntimeError(f'Failed estimation in {st["name"]}, replication {rep}; no failed fits were discarded')
            for row in result:
                row.update(rep=rep,stratum=st['name'],phase=task['phase'])
            rows.extend(result)
    pd.DataFrame(rows).to_csv(target,index=False,float_format='%.17g',compression='gzip')
    stamp.write_text(json.dumps(dict(task=task,task_sha256=task_hash,
        result_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        rows=len(rows),seconds=time.perf_counter()-start),indent=2),encoding='utf-8')
    return str(target)

def summarize(paths,design):
    frames=[pd.read_csv(path,low_memory=False) for path in paths]
    df=pd.concat(frames,ignore_index=True).copy()
    if df.duplicated(['phase','stratum','scenario','method','rep']).any():
        raise RuntimeError('Duplicate simulation records')
    lookup={c['id']:c for c in design['cases']}
    L,U=design['range']
    left=df.scenario.map(lambda name:lookup[name]['b']*(lookup[name]['a1']+2*lookup[name]['a2']*L))
    right=df.scenario.map(lambda name:lookup[name]['b']*(lookup[name]['a1']+2*lookup[name]['a2']*U))
    pos=(left>1e-12)&(right< -1e-12)
    neg=(left< -1e-12)&(right>1e-12)
    df['reversal']=df.either_posneg.eq(1)|df.either_negpos.eq(1)
    df['false_reversal']=df.reversal&~(pos|neg)
    df['correct_direction']=(df.either_posneg.eq(1)&pos)|(df.either_negpos.eq(1)&neg)
    df['wrong_direction']=(df.either_posneg.eq(1)&neg)|(df.either_negpos.eq(1)&pos)
    names=['a1','a2','bpath','c1','c2','eta','kappa_m','kappa_d','kappa_t','ie_l','ie_u','delta_ie','ie1','ie2','ie3']
    import ols
    for j,name in enumerate(names):
        truth=df.scenario.map(lambda s:ols.evaluate(np.zeros((1,5)),lookup[s])[0,j])
        if name+'_lo' in df:
            available=df[name+'_lo'].notna()&df[name+'_hi'].notna()
            df[name+'_cover']=np.where(available,(df[name+'_lo']<=truth)&(truth<=df[name+'_hi']),np.nan)
            df[name+'_width']=df[name+'_hi']-df[name+'_lo']
    metrics=['kappa_m_cover','kappa_d_cover','kappa_t_cover','ie_l_cover','ie_u_cover',
             'band_cover','component_joint_cover','component_slope_mid_cover',
             'component_curvature_xm_cover','component_bpath_cover','component_kappa_d_cover',
             'reversal','false_reversal','correct_direction','wrong_direction']
    z=norm.isf(.025)
    result=[]
    for key,group in df.groupby(['phase','stratum','scenario','method'],sort=True):
        row=dict(zip(['phase','stratum','scenario','method'],key))
        row.update(R=len(group),valid=int(group.outer_ok.eq(1).sum()),failed=int(group.outer_ok.ne(1).sum()),case=lookup[key[2]]['case'])
        for metric in metrics:
            if metric not in group:
                continue
            values=group.loc[group.outer_ok.eq(1),metric].dropna()
            if not len(values):
                continue
            n=len(values)
            p=float(values.mean())
            den=1+z*z/n
            center=(p+z*z/(2*n))/den
            half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
            row[metric]=p
            row[metric+'_lo']=max(0.,center-half)
            row[metric+'_hi']=min(1.,center+half)
        for metric in ['kappa_m_width','ie_l_width','ie_u_width','band_mean_width','max_influence','wild_invalid']:
            if metric in group:
                values=group[metric].dropna()
                if len(values):
                    row[metric+'_mean']=float(values.mean())
                    row[metric+'_median']=float(values.median())
        result.append(row)
    return pd.DataFrame(result)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase',choices=(*PHASES,'all'),default='all')
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--output',type=Path,default=HERE/'results')
    parser.add_argument('--smoke',action='store_true',help='Run three samples of the G=5,N=400 design with 199 draws')
    args=parser.parse_args()
    if args.workers<1:
        parser.error('--workers must be positive')
    design=json.loads((HERE/'design.json').read_text(encoding='utf-8'))
    if args.smoke:
        tasks=[dict(phase='smoke',stratum='g5_n400',start=1,stop=3,B=199,
                    wild=True,pairs=True,case_names=['regular','zero_b','boundary'])]
    else:
        tasks=[t for t in design['jobs'] if args.phase=='all' or t['phase']==args.phase]
    args.output.mkdir(parents=True,exist_ok=True)
    paths=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_job,task,design,str(args.output)) for task in tasks]
        for future in as_completed(futures):
            paths.append(future.result())
            print(f'Completed {len(paths)}/{len(tasks)} tasks',flush=True)
    target=args.output/('smoke_summary.csv' if args.smoke else f'{args.phase}_summary.csv')
    summarize(paths,design).to_csv(target,index=False,float_format='%.12g')
    print(target)

if __name__=='__main__':
    main()
