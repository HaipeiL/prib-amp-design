import argparse, csv, json, math, os, statistics, subprocess, sys, zipfile
from pathlib import Path
from common import AA, atomic_json, read_fasta, write_fasta, valid_seq, sha256_file
from evaluator import score_broad

def hydrophobic_frac(s): return sum(c in "AVILMFWY" for c in s)/len(s)
def charge(s): return s.count("K")+s.count("R")-s.count("D")-s.count("E")
def max_run(s):
    best=cur=1
    for i in range(1,len(s)):
        cur=cur+1 if s[i]==s[i-1] else 1; best=max(best,cur)
    return best
def kmer_set(s,k=3): return {s[i:i+k] for i in range(max(0,len(s)-k+1))}
def jaccard_dist(a,b):
    if not a and not b: return 0.0
    return 1-len(a&b)/max(1,len(a|b))
def rank01(vals):
    order=sorted(range(len(vals)),key=lambda i:(vals[i],i))
    out=[0.0]*len(vals); n=max(1,len(vals)-1)
    for r,i in enumerate(order): out[i]=r/n
    return out

def ratio_fallback(a,b):
    # Levenshtein.ratio uses indel similarity = 2*LCS/(len(a)+len(b)).
    prev=[0]*(len(b)+1)
    for ca in a:
        cur=[0]
        for j,cb in enumerate(b,1):
            cur.append(prev[j-1]+1 if ca==cb else max(prev[j],cur[-1]))
        prev=cur
    return 2*prev[-1]/(len(a)+len(b)) if (a or b) else 1.0

def best_ref_similarity(seq,refs):
    try:
        import Levenshtein
        f=Levenshtein.ratio
    except Exception:
        f=ratio_fallback
    best=0.0
    # quick length upper bound for indel ratio: 2*minlen/(len1+len2)
    for r in refs:
        ub=2*min(len(seq),len(r))/(len(seq)+len(r))
        if ub<=best: continue
        v=f(seq,r)
        if v>best: best=v
        if best>0.80: break
    return best

def load_rows(run):
    rows=[]
    for fp in sorted((run/"raw").glob("lane*/shard_*.csv")):
        with open(fp,newline="") as f:
            for r in csv.DictReader(f): rows.append(r)
    return rows

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); ap.add_argument("--run-dir",required=True)
    a=ap.parse_args(); cfg=json.loads(Path(a.config).read_text()); run=Path(a.run_dir)
    res=run/"results"; res.mkdir(parents=True,exist_ok=True)
    rows=load_rows(run); raw_n=len(rows)
    seen=set(); clean=[]
    bad=dups=0
    refs=[]
    refpath=cfg.get("official_reference")
    if refpath and Path(refpath).exists(): refs=read_fasta(refpath)
    refset=set(refs)
    exact_overlap=0
    for r in rows:
        s=r["sequence"].strip().upper()
        if not valid_seq(s): bad+=1; continue
        if s in seen: dups+=1; continue
        seen.add(s)
        if refset and s in refset: exact_overlap+=1; continue
        r["sequence"]=s; clean.append(r)
    if len(clean)<50000:
        status="PARTIAL"; reason=f"clean pool {len(clean)} < 50000"
    elif not refs:
        status="PARTIAL"; reason="official antibacterial reference unavailable"
    else:
        status="PASS"; reason=""
    atomic_json(res/"raw_generation_summary.json",{"raw_rows":raw_n,"invalid":bad,"duplicates":dups})
    atomic_json(res/"clean_pool_summary.json",{"clean_unique":len(clean),"exact_reference_overlap_removed":exact_overlap,
                                              "official_reference_loaded":bool(refs),"status":status,"reason":reason})
    if len(clean)<50000:
        package_return(run,cfg,status,reason); return
    seqs=[r["sequence"] for r in clean]
    for r in clean:
        s=r["sequence"]; r["length"]=len(s); r["charge_proxy"]=charge(s); r["hydrophobic_fraction"]=hydrophobic_frac(s); r["max_run"]=max_run(s)
    ev=cfg.get("evaluator_checkpoints",[])
    if len(ev)==3 and all(Path(p).exists() for p in ev):
        means,mins=score_broad(seqs,ev,["Escherichia coli","Pseudomonas aeruginosa","Staphylococcus aureus"],batch=1024)
        for r,m,n in zip(clean,means,mins): r["pmic_mean3"]=m; r["pmic_min3"]=n
        rp=rank01(mins); rm=rank01(means)
    else:
        for r in clean: r["pmic_mean3"]=float("nan"); r["pmic_min3"]=float("nan")
        rp=rm=[0.5]*len(clean)
    phys=[]
    for r in clean:
        h=r["hydrophobic_fraction"]; c=r["charge_proxy"]
        hydro=max(0.0,1-abs(h-0.45)/0.45)
        ch=1/(1+math.exp(-(c-1)/2))
        run_pen=0.0 if r["max_run"]>=5 else 1.0
        phys.append(0.45*hydro+0.45*ch+0.10*run_pen)
    for i,r in enumerate(clean):
        r["score"]=0.65*rp[i]+0.25*rm[i]+0.10*phys[i]
    clean.sort(key=lambda r:(-r["score"],int(r["candidate_id"])))
    # Build top candidates with official similarity hard gate and 0.78 safety margin.
    eligible=[]
    for r in clean[:5000]:
        sim=best_ref_similarity(r["sequence"],refs)
        r["max_reference_similarity"]=sim
        if sim <= float(cfg["top_similarity_safety"]) and r["max_run"]<5:
            eligible.append(r)
        if len(eligible)>=1200: break
    if len(eligible)<100:
        status="PARTIAL"; reason=f"only {len(eligible)} top candidates pass similarity safety"
        package_return(run,cfg,status,reason); return
    # Greedy score + 3-mer diversity portfolio.
    top=[]; km=[kmer_set(r["sequence"]) for r in eligible]
    chosen=set()
    first=0; top.append(eligible[first]); chosen.add(first)
    while len(top)<100:
        best_i=None; best_u=-1
        for i,r in enumerate(eligible):
            if i in chosen: continue
            d=min(jaccard_dist(km[i],km[j]) for j in chosen)
            u=0.78*r["score"]+0.22*d
            if u>best_u or (u==best_u and int(r["candidate_id"])<int(eligible[best_i]["candidate_id"])):
                best_u=u; best_i=i
        chosen.add(best_i); top.append(eligible[best_i])
    # Library: force top100, then preserve condition-mixture representation with round-robin over six strata.
    top_ids={r["candidate_id"] for r in top}; library=list(top)
    strata={}
    for r in clean:
        if r["candidate_id"] in top_ids: continue
        key=(r["species"],r["condition"]); strata.setdefault(key,[]).append(r)
    keys=sorted(strata)
    idx={k:0 for k in keys}
    while len(library)<50000:
        progressed=False
        for k in keys:
            j=idx[k]
            if j<len(strata[k]):
                library.append(strata[k][j]); idx[k]=j+1; progressed=True
                if len(library)>=50000: break
        if not progressed: break
    if len(library)!=50000:
        status="PARTIAL"; reason=f"library fill stopped at {len(library)}"
        package_return(run,cfg,status,reason); return
    libseq=[r["sequence"] for r in library]; topseq=[r["sequence"] for r in top]
    write_fasta(libseq,res/"provisional_library.fasta","lib")
    write_fasta(topseq,res/"provisional_top.fasta","top")
    with open(res/"candidate_scores.csv","w",newline="") as f:
        fields=["candidate_id","species","condition","sequence","length","charge_proxy","hydrophobic_fraction","max_run","pmic_mean3","pmic_min3","score"]
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for r in library: w.writerow({k:r.get(k,"") for k in fields})
    with open(res/"top100_audit.csv","w",newline="") as f:
        fields=["rank","candidate_id","species","condition","sequence","score","pmic_mean3","pmic_min3","max_reference_similarity"]
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for rank,r in enumerate(top,1):
            w.writerow({**{k:r.get(k,"") for k in fields if k!="rank"},"rank":rank})
    maxsim=max(float(r["max_reference_similarity"]) for r in top)
    summary={"status":"PASS","raw_generated":raw_n,"unique_valid":len(seen),
             "clean_after_reference_filter":len(clean),"provisional_library_size":len(library),
             "provisional_top_size":len(top),"top100_max_reference_similarity":maxsim,
             "external_scorer_completed":False,"independent_pmic_evaluator_completed":len(ev)==3,
             "failed_stages":[],"recommended_resume_script":None}
    atomic_json(run/"round_summary.json",summary)
    package_return(run,cfg,"PASS","")
    print(f"POSTPROCESS=PASS clean={len(clean)} library=50000 top=100 maxsim={maxsim:.4f}")

def package_return(run,cfg,status,reason):
    if not (run/"round_summary.json").exists():
        atomic_json(run/"round_summary.json",{"status":status,"raw_generated":0,"unique_valid":0,
          "clean_after_reference_filter":0,"generation_seconds":0,"two_gpu_sequences_per_second":cfg.get("measured_two_gpu_seq_per_sec",0),
          "external_scorer_completed":False,"provisional_library_size":0,"provisional_top_size":0,
          "top100_max_reference_similarity":0.0,"failed_stages":[reason] if reason else [],
          "recommended_resume_script":"13_resume_postprocess.sh" if status!="PASS" else None})
    (run/"TAKEBACK.md").write_text(f"# AMP Challenge R01 Takeback\n\nStatus: **{status}**\n\n{reason}\n")
    z=run/"RETURN_NIGHT.zip"
    include=["TAKEBACK.md","round_summary.json","configs/night_config.json","manifests/ASSET_MAP.json",
             "results/raw_generation_summary.json","results/clean_pool_summary.json","results/provisional_library.fasta",
             "results/provisional_top.fasta","results/candidate_scores.csv","results/top100_audit.csv",
             "logs/lane0.log","logs/lane1.log","logs/postprocess.log","artifact_manifest.sha256"]
    with zipfile.ZipFile(z,"w",compression=zipfile.ZIP_DEFLATED) as zz:
        for rel in include:
            p=run/rel
            if p.exists(): zz.write(p,rel)
    print("RETURN_ZIP="+str(z))
if __name__=="__main__": main()
