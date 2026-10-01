import argparse,csv,json,math,statistics
from pathlib import Path
from evaluator import load_ensemble,score_detailed_loaded
from common import valid_seq,atomic_json

SPECIES=["Escherichia coli","Pseudomonas aeruginosa","Staphylococcus aureus"]
HYDRO=set("AVILMFWY"); POS=set("KRH"); NEG=set("DE"); POLAR=set("STNQ"); SPECIAL=set("GPC")

def contrast(aa):
    if aa in HYDRO: return "S"
    if aa in POS: return "Q"
    if aa in NEG: return "N"
    if aa in POLAR: return "L"
    return "K"

def perturbations(s):
    # Deterministic exhaustive position-wise two-probe scan + terminal truncation.
    out=[]; seen={s}
    for i,a in enumerate(s):
        aa="A" if a!="A" else "G"
        v=s[:i]+aa+s[i+1:]
        if v not in seen: out.append(("alanine_scan",i,v)); seen.add(v)
        aa2=contrast(a)
        if aa2==a: aa2="A" if a!="A" else "G"
        v=s[:i]+aa2+s[i+1:]
        if v not in seen: out.append(("contrast_scan",i,v)); seen.add(v)
    for k in (1,2):
        if len(s)-k>=8:
            v=s[k:]
            if v not in seen: out.append((f"n_trunc_{k}",-1,v)); seen.add(v)
            v=s[:-k]
            if v not in seen: out.append((f"c_trunc_{k}",-1,v)); seen.add(v)
    return out

def pred_min3(x): return min(statistics.mean(v) for v in x)
def pred_model_min3(x,e): return min(x[s][e] for s in range(len(x)))
def q(vals,p):
    vals=sorted(vals)
    if not vals: return float("nan")
    x=(len(vals)-1)*p; lo=int(math.floor(x)); hi=int(math.ceil(x))
    if lo==hi:return vals[lo]
    return vals[lo]*(hi-x)+vals[hi]*(x-lo)

def summarize_candidate(r,perts,preds,delta,catdelta):
    base=float(r["pmic_min3"]); vals=[pred_min3(x) for x in preds]; drops=[base-v for v in vals]
    nmodels=len(preds[0][0]) if preds else 0
    model_p10=[q([pred_model_min3(x,e) for x in preds],0.10) for e in range(nmodels)]
    by={}
    for (kind,pos,seq),v in zip(perts,vals): by.setdefault(kind,[]).append(v)
    worst_i=max(range(len(drops)),key=lambda i:drops[i]) if drops else 0
    return {
      "candidate_id":r["candidate_id"],"n_perturb":len(vals),"perturb_p10_min3":q(vals,0.10),"perturb_p25_min3":q(vals,0.25),
      "perturb_median_min3":q(vals,0.50),"drop_p90":q(drops,0.90),"drop_worst":max(drops) if drops else 0.0,
      "retention_fraction_0p5":sum(v>=base-delta for v in vals)/max(1,len(vals)),
      "catastrophic_fraction_1p0":sum(v<base-catdelta for v in vals)/max(1,len(vals)),
      "perturb_crossmodel_p10_floor":min(model_p10) if model_p10 else float("nan"),
      "alanine_p10_min3":q(by.get("alanine_scan",[]),0.10),"contrast_p10_min3":q(by.get("contrast_scan",[]),0.10),
      "trunc_p10_min3":q([v for k,vs in by.items() if "trunc" in k for v in vs],0.10),
      "worst_perturb_kind":perts[worst_i][0] if vals else "","worst_perturb_position":perts[worst_i][1] if vals else "",
      "worst_perturb_sequence":perts[worst_i][2] if vals else "","worst_perturb_pmic_min3":vals[worst_i] if vals else ""
    }

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); ap.add_argument("--run-dir",required=True); ap.add_argument("--lane",type=int,required=True)
    a=ap.parse_args(); cfg=json.loads(Path(a.config).read_text()); run=Path(a.run_dir); prev=Path(cfg["previous_run"])
    ncfg=json.loads((prev/"configs/night_config.json").read_text()); lane=a.lane
    with open(run/"work/frontier.csv",newline="") as f: rows=[r for r in csv.DictReader(f) if int(r["candidate_id"])%2==lane]
    outdir=run/"work"/f"lane{lane}"; outdir.mkdir(parents=True,exist_ok=True); shard=int(cfg["lane_candidate_shard"])
    models,smap=load_ensemble(ncfg["evaluator_checkpoints"],device="cuda:0")
    fields=["candidate_id","n_perturb","perturb_p10_min3","perturb_p25_min3","perturb_median_min3","drop_p90","drop_worst","retention_fraction_0p5",
      "catastrophic_fraction_1p0","perturb_crossmodel_p10_floor","alanine_p10_min3","contrast_p10_min3","trunc_p10_min3","worst_perturb_kind","worst_perturb_position","worst_perturb_sequence","worst_perturb_pmic_min3"]
    done=0
    for st in range(0,len(rows),shard):
        rr=rows[st:st+shard]; fp=outdir/f"shard_{st:06d}_{st+len(rr):06d}.csv"
        if fp.exists(): done+=len(rr); continue
        allp=[]; bounds=[]
        for r in rr:
            p=perturbations(r["sequence"]); bounds.append((len(allp),len(allp)+len(p),p)); allp.extend([x[2] for x in p])
        preds=score_detailed_loaded(allp,models,smap,SPECIES,batch=int(cfg["eval_batch"]),device="cuda:0")
        tmp=fp.with_suffix(".csv.tmp")
        with open(tmp,"w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
            for r,(lo,hi,p) in zip(rr,bounds):
                w.writerow(summarize_candidate(r,p,preds[lo:hi],float(cfg["retention_delta_pmic"]),float(cfg["catastrophic_delta_pmic"])))
        tmp.replace(fp); done+=len(rr); atomic_json(outdir/"progress.json",{"lane":lane,"completed":done,"target":len(rows),"last_shard":fp.name})
        print(f"[prib lane{lane}] {done}/{len(rows)}",flush=True)
    atomic_json(outdir/"complete.json",{"status":"PASS","lane":lane,"rows":len(rows)})
    print(f"PRIB_LANE{lane}=PASS rows={len(rows)}",flush=True)
if __name__=="__main__": main()
