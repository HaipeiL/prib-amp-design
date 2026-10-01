import argparse,csv,json,statistics
from pathlib import Path
import torch
from evaluator import load_ensemble,score_detailed_loaded


def load_rows(p):
    with open(p,newline="") as f:return list(csv.DictReader(f))

def wmean(vals,weights):
    num=sum(float(vals[k])*float(v) for k,v in weights.items()); den=sum(float(v) for v in weights.values())
    return num/den

def weighted_values(vals,weights):
    out=[]
    for k,w in weights.items(): out.extend([float(vals[k])]*int(w))
    return out

def quantile(xs,q):
    xs=sorted(float(x) for x in xs)
    if not xs:return float('nan')
    if len(xs)==1:return xs[0]
    p=(len(xs)-1)*q; lo=int(p); hi=min(lo+1,len(xs)-1); a=p-lo
    return xs[lo]*(1-a)+xs[hi]*a

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config',required=True); ap.add_argument('--lane',type=int,required=True); ap.add_argument('--run-dir',required=True)
    a=ap.parse_args(); cfg=json.loads(Path(a.config).read_text()); run=Path(a.run_dir); r02=Path(cfg['r02_run']); r01=Path(cfg['r01_run'])
    rows=load_rows(r02/'results/frontier_audit.csv')
    rows=[r for i,r in enumerate(rows) if i%2==a.lane]
    night=json.loads((r01/'configs/night_config.json').read_text()); ckpts=night['evaluator_checkpoints']
    device='cuda:0'; models,smap=load_ensemble(ckpts,device)
    missing=[sp for sp in cfg['official_species'] if sp not in smap]
    if missing: raise RuntimeError('official species missing from frozen evaluator map: '+repr(missing))
    outp=run/f'work/lane{a.lane}/panel_scores.csv'; outp.parent.mkdir(parents=True,exist_ok=True)
    fields=['candidate_id','sequence','r02_status','r02_reliability_floor','r02_reliability_geomean','physics_caution','physics_severe',
            'broad_mean','broad_p10','gram_negative_mean','gram_positive_mean','mdr_mean','panel_species_floor','panel_crossmodel_floor','panel_max_model_sd']
    fields += [f"pmic_{sp.replace(' ','_')}" for sp in cfg['official_species']]
    with open(outp,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for st in range(0,len(rows),int(cfg['eval_batch'])):
            chunk=rows[st:st+int(cfg['eval_batch'])]; seqs=[r['sequence'] for r in chunk]
            pred=score_detailed_loaded(seqs,models,smap,cfg['official_species'],batch=int(cfg['eval_batch']),device=device)
            for r,x in zip(chunk,pred):
                means={}; all_preds=[]; maxsd=0.0
                for sp,ens in zip(cfg['official_species'],x):
                    ens=[float(v) for v in ens]; means[sp]=statistics.mean(ens); all_preds.extend(ens)
                    maxsd=max(maxsd,statistics.pstdev(ens))
                broadvals=weighted_values(means,cfg['broad_weights'])
                rec={
                    'candidate_id':r['candidate_id'],'sequence':r['sequence'],'r02_status':r['audit_status'],
                    'r02_reliability_floor':r['reliability_floor'],'r02_reliability_geomean':r['reliability_geomean'],
                    'physics_caution':r.get('physics_caution',''),'physics_severe':r.get('physics_severe',''),
                    'broad_mean':wmean(means,cfg['broad_weights']),'broad_p10':quantile(broadvals,0.10),
                    'gram_negative_mean':wmean(means,cfg['gram_negative_weights']),'gram_positive_mean':wmean(means,cfg['gram_positive_weights']),
                    'mdr_mean':wmean(means,cfg['mdr_weights']),'panel_species_floor':min(means.values()),
                    'panel_crossmodel_floor':min(all_preds),'panel_max_model_sd':maxsd,
                }
                for sp in cfg['official_species']: rec[f"pmic_{sp.replace(' ','_')}"]=means[sp]
                w.writerow(rec)
            f.flush(); print(f'[lane{a.lane}] {min(st+len(chunk),len(rows))}/{len(rows)}',flush=True)
    print(f'PANEL_LANE{a.lane}=PASS rows={len(rows)}',flush=True)
if __name__=='__main__': main()
