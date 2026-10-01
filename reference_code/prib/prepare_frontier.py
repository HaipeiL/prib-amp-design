import argparse,csv,json,math,statistics
from pathlib import Path
from common import read_fasta,valid_seq,descriptors,physics_flags,rank01,atomic_json
from evaluator import score_detailed

SPECIES=["Escherichia coli","Pseudomonas aeruginosa","Staphylococcus aureus"]

def load_raw(prev):
    rows=[]
    for fp in sorted((prev/"raw").glob("lane*/shard_*.csv")):
        with open(fp,newline="") as f: rows.extend(csv.DictReader(f))
    return rows

def summarize_pred(x):
    # x species x ensemble
    sp_means=[statistics.mean(v) for v in x]
    sp_sds=[statistics.pstdev(v) for v in x]
    flat=[z for v in x for z in v]
    return {
      "pmic_mean3":statistics.mean(sp_means),"pmic_min3":min(sp_means),"pmic_max3":max(sp_means),
      "species_spread":max(sp_means)-min(sp_means),"ensemble_std_max":max(sp_sds),
      "ensemble_std_mean":statistics.mean(sp_sds),"ensemble_range":max(flat)-min(flat),
      "pmic_model_floor":min(min(x[s][e] for s in range(3)) for e in range(len(x[0]))),
      **{f"pmic_{i}":sp_means[i] for i in range(3)}
    }

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); ap.add_argument("--run-dir",required=True)
    a=ap.parse_args(); cfg=json.loads(Path(a.config).read_text()); run=Path(a.run_dir); prev=Path(cfg["previous_run"])
    ncfg=json.loads((prev/"configs/night_config.json").read_text()); refs=set(read_fasta(ncfg["official_reference"]))
    raw=load_raw(prev); seen=set(); clean=[]; bad=dups=exact=0
    for r in raw:
        s=r["sequence"].strip().upper()
        if not valid_seq(s): bad+=1; continue
        if s in seen: dups+=1; continue
        seen.add(s)
        if s in refs: exact+=1; continue
        clean.append(r)
    seqs=[r["sequence"] for r in clean]
    preds=score_detailed(seqs,ncfg["evaluator_checkpoints"],SPECIES,batch=int(cfg["eval_batch"]))
    for r,x in zip(clean,preds):
        r.update(summarize_pred(x)); d=descriptors(r["sequence"]); r.update(d)
        caution,severe=physics_flags(d,cfg); r["physics_caution"]=";".join(caution); r["physics_severe"]=";".join(severe)
    # Reproduce activity-oriented base ranking, but preserve evaluator uncertainty separately.
    rp=rank01([r["pmic_min3"] for r in clean]); rm=rank01([r["pmic_mean3"] for r in clean])
    phys=[]
    for r in clean:
        h=r["hydrophobic_fraction"]; c=r["charge_proxy"]
        hydro=max(0.0,1-abs(h-0.45)/0.45); ch=1/(1+math.exp(-(c-1)/2)); run_pen=0.0 if r["max_residue_run"]>=5 else 1.0
        phys.append(0.45*hydro+0.45*ch+0.10*run_pen)
    for i,r in enumerate(clean): r["base_score"]=0.65*rp[i]+0.25*rm[i]+0.10*phys[i]
    clean.sort(key=lambda r:(-r["base_score"],int(r["candidate_id"])))
    # 6-stratum representation + global high-score fill. This is screening breadth, not a final quota.
    picked={}; per=int(cfg["frontier_per_stratum"])
    for r in clean:
        k=(r["species"],r["condition"]); picked.setdefault(k,[])
        if len(picked[k])<per: picked[k].append(r)
    frontier=[]; ids=set()
    for k in sorted(picked):
        for r in picked[k]:
            if r["candidate_id"] not in ids: frontier.append(r); ids.add(r["candidate_id"])
    # guarantee prior Top100 is stress-tested
    oldtop=read_fasta(prev/"results/provisional_top.fasta") if (prev/"results/provisional_top.fasta").exists() else []
    oldset=set(oldtop)
    byseq={r["sequence"]:r for r in clean}
    for s in oldtop:
        if s in byseq and byseq[s]["candidate_id"] not in ids:
            frontier.append(byseq[s]); ids.add(byseq[s]["candidate_id"])
    target=int(cfg["frontier_total"])
    for r in clean:
        if len(frontier)>=target: break
        if r["candidate_id"] not in ids: frontier.append(r); ids.add(r["candidate_id"])
    frontier.sort(key=lambda r:int(r["candidate_id"]))
    out=run/"work"; out.mkdir(parents=True,exist_ok=True)
    fields=["candidate_id","species","condition","sequence","base_score","pmic_mean3","pmic_min3","pmic_max3","species_spread",
      "ensemble_std_max","ensemble_std_mean","ensemble_range","pmic_model_floor","pmic_0","pmic_1","pmic_2","length","charge_proxy","charge_density",
      "hydrophobic_fraction","max_hydrophobic_run","max_residue_run","sequence_entropy","hydrophobic_moment_alpha","cysteine_count","methionine_count",
      "physics_caution","physics_severe","was_r01_top100"]
    with open(out/"frontier.csv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for r in frontier:
            rr={k:r.get(k,"") for k in fields}; rr["was_r01_top100"]=int(r["sequence"] in oldset); w.writerow(rr)
    atomic_json(run/"results/prepare_summary.json",{
      "raw":len(raw),"clean":len(clean),"duplicates":dups,"invalid":bad,"exact_reference_removed":exact,
      "frontier":len(frontier),"old_top100_in_frontier":sum(r["sequence"] in oldset for r in frontier),"status":"PASS"})
    print(f"PREPARE=PASS clean={len(clean)} frontier={len(frontier)}",flush=True)
if __name__=="__main__": main()
