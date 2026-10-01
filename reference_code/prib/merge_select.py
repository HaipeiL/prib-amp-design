import argparse,csv,json,math,statistics
from pathlib import Path
from common import read_fasta,write_fasta,rank01,physics_flags,descriptors,kmer_set,jaccard_similarity,atomic_json,sha256_file

try:
    import Levenshtein
except Exception as e:
    raise SystemExit("Levenshtein package required for R02 reference scan: "+str(e))

def best_ref_similarity(seq,refs,stop=0.80):
    best=0.0
    for r in refs:
        ub=2*min(len(seq),len(r))/(len(seq)+len(r))
        if ub<=best: continue
        v=Levenshtein.ratio(seq,r)
        if v>best: best=v
        if best>stop: break
    return best

def gmean(vals):
    vals=[max(1e-9,float(v)) for v in vals]
    return math.exp(sum(math.log(v) for v in vals)/len(vals))

def load_csv(p):
    with open(p,newline="") as f:return list(csv.DictReader(f))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); ap.add_argument("--run-dir",required=True)
    a=ap.parse_args(); cfg=json.loads(Path(a.config).read_text()); run=Path(a.run_dir); prev=Path(cfg["previous_run"])
    rows=load_csv(run/"work/frontier.csv"); pert={}
    for fp in sorted((run/"work").glob("lane*/shard_*.csv")):
        for r in load_csv(fp): pert[r["candidate_id"]]=r
    if len(pert)!=len(rows): raise RuntimeError(f"perturb summaries {len(pert)} != frontier {len(rows)}")
    numeric=["base_score","pmic_mean3","pmic_min3","pmic_max3","species_spread","ensemble_std_max","ensemble_std_mean","ensemble_range",
      "pmic_model_floor","charge_density","hydrophobic_fraction","sequence_entropy","hydrophobic_moment_alpha"]
    pn=["perturb_p10_min3","perturb_p25_min3","perturb_median_min3","drop_p90","drop_worst","retention_fraction_0p5","catastrophic_fraction_1p0",
      "perturb_crossmodel_p10_floor","alanine_p10_min3","contrast_p10_min3","trunc_p10_min3"]
    for r in rows:
        for k in numeric:
            try:r[k]=float(r[k])
            except:pass
        r.update(pert[r["candidate_id"]])
        for k in pn:
            try:r[k]=float(r[k])
            except:r[k]=float("nan")
        d=descriptors(r["sequence"]); caution,severe=physics_flags(d,cfg)
        r["physics_caution"]=";".join(caution); r["physics_severe"]=";".join(severe)
    # Percentile axes: potency, perturbed lower-tail potency, local retention, model consensus.
    axes={
      "p_potency":rank01([r["pmic_min3"] for r in rows]),
      "p_perturb":rank01([r["perturb_crossmodel_p10_floor"] for r in rows]),
      "p_retention":rank01([r["retention_fraction_0p5"] for r in rows]),
      "p_consensus":rank01([-r["ensemble_std_max"] for r in rows]),
    }
    for i,r in enumerate(rows):
        vals=[axes[k][i] for k in ("p_potency","p_perturb","p_retention","p_consensus")]
        for k in axes:r[k]=axes[k][i]
        r["reliability_floor"]=min(vals); r["reliability_geomean"]=gmean(vals)
        # Status is relative audit triage, not a biological safety label.
        if r["physics_severe"]: r["audit_status"]="REJECT"
        elif r["reliability_floor"]>=0.50: r["audit_status"]="PASS"
        elif r["reliability_floor"]>=0.25: r["audit_status"]="CAUTION"
        else:r["audit_status"]="REJECT"
    rows.sort(key=lambda r:(-r["reliability_floor"],-r["reliability_geomean"],-r["pmic_min3"],int(r["candidate_id"])))
    # Only expensive official-reference scan on the strongest audit pool.
    nscan=min(int(cfg["reference_scan_pool"]),len(rows)); refs=read_fasta(json.loads((prev/"configs/night_config.json").read_text())["official_reference"])
    eligible=[]
    for r in rows[:nscan]:
        if r["audit_status"]=="REJECT": continue
        sim=best_ref_similarity(r["sequence"],refs)
        r["max_reference_similarity"]=sim
        if sim<=float(cfg["top_similarity_safety"]): eligible.append(r)
    if len(eligible)<int(cfg["minimum_pass_pool"]):
        raise RuntimeError(f"only {len(eligible)} eligible after PRIB+reference gates")
    # Diversity is a portfolio gate, not part of the reliability score.
    top=[]; km=[]; gate=float(cfg["pairwise_kmer_similarity_gate"]); relax=float(cfg["pairwise_kmer_similarity_relaxed"])
    def fill(threshold):
        for r in eligible:
            if len(top)>=int(cfg["top_size"]): break
            if any(r["candidate_id"]==x["candidate_id"] for x in top): continue
            k=kmer_set(r["sequence"])
            if km and max(jaccard_similarity(k,z) for z in km)>threshold: continue
            top.append(r); km.append(k)
    fill(gate)
    if len(top)<int(cfg["top_size"]): fill(relax)
    if len(top)<int(cfg["top_size"]):
        for r in eligible:
            if len(top)>=int(cfg["top_size"]):break
            if all(r["candidate_id"]!=x["candidate_id"] for x in top): top.append(r); km.append(kmer_set(r["sequence"]))
    if len(top)!=int(cfg["top_size"]): raise RuntimeError(f"top selection stopped at {len(top)}")
    # Final 50k = PRIB top100 + frozen R01 50k reservoir, preserving official constraints.
    oldlib=read_fasta(prev/"results/provisional_library.fasta"); topseq=[r["sequence"] for r in top]; seen=set(); lib=[]
    for s in topseq+oldlib:
        if s not in seen: seen.add(s); lib.append(s)
        if len(lib)>=int(cfg["library_size"]):break
    if len(lib)!=int(cfg["library_size"]): raise RuntimeError(f"library fill {len(lib)}")
    out=run/"results"; out.mkdir(parents=True,exist_ok=True)
    write_fasta(lib,out/"prib_library.fasta","lib"); write_fasta(topseq,out/"prib_top.fasta","top")
    oldtop=set(read_fasta(prev/"results/provisional_top.fasta")); overlap=sum(s in oldtop for s in topseq)
    fields=["rank","candidate_id","species","condition","sequence","audit_status","pmic_min3","pmic_mean3","ensemble_std_max","perturb_p10_min3","perturb_median_min3",
      "retention_fraction_0p5","catastrophic_fraction_1p0","perturb_crossmodel_p10_floor","drop_p90","drop_worst","p_potency","p_perturb","p_retention","p_consensus","reliability_floor","reliability_geomean",
      "physics_caution","physics_severe","max_reference_similarity","worst_perturb_kind","worst_perturb_position","worst_perturb_sequence","worst_perturb_pmic_min3","was_r01_top100"]
    with open(out/"prib_top100_audit.csv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for i,r in enumerate(top,1):
            rr={k:r.get(k,"") for k in fields}; rr["rank"]=i; w.writerow(rr)
    allfields=["candidate_id","species","condition","sequence","audit_status","pmic_min3","pmic_mean3","ensemble_std_max","perturb_p10_min3","perturb_median_min3","retention_fraction_0p5",
      "catastrophic_fraction_1p0","perturb_crossmodel_p10_floor","drop_p90","drop_worst","p_potency","p_perturb","p_retention","p_consensus","reliability_floor","reliability_geomean","physics_caution","physics_severe","was_r01_top100"]
    with open(out/"frontier_audit.csv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=allfields); w.writeheader()
        for r in rows:w.writerow({k:r.get(k,"") for k in allfields})
    summary={
      "status":"PASS","frontier_size":len(rows),"eligible_after_prib_and_reference":len(eligible),"top_size":len(top),"library_size":len(lib),
      "r01_top100_overlap":overlap,"r01_top100_replaced":100-overlap,
      "top100_max_reference_similarity":max(float(r["max_reference_similarity"]) for r in top),
      "top100_median_pmic_min3":statistics.median(float(r["pmic_min3"]) for r in top),
      "top100_median_perturb_p10_min3":statistics.median(float(r["perturb_p10_min3"]) for r in top),
      "top100_median_perturb_crossmodel_p10_floor":statistics.median(float(r["perturb_crossmodel_p10_floor"]) for r in top),
      "top100_median_retention_0p5":statistics.median(float(r["retention_fraction_0p5"]) for r in top),
      "top100_median_ensemble_std_max":statistics.median(float(r["ensemble_std_max"]) for r in top),
      "selection_principle":"official gates -> PRIB relative audit -> maximin reliability floor -> diversity gate",
      "safety_note":"Physics descriptors and perturbation robustness are pre-experimental risk proxies, not HC50 or biological safety predictions."
    }
    atomic_json(out/"prib_summary.json",summary)
    report=f'''# PRIB Experimental-Readiness Audit R02\n\nStatus: **PASS**\n\n- Frontier stress-tested: {len(rows)}\n- Eligible after PRIB + official reference gate: {len(eligible)}\n- Final library: {len(lib)}\n- Final top list: {len(top)}\n- R01 Top100 retained: {overlap}/100\n- R01 Top100 replaced by the audit: {100-overlap}/100\n- Max reference similarity in final Top100: {summary['top100_max_reference_similarity']:.6f}\n- Median broad-spectrum baseline pMIC proxy (min over 3 species): {summary['top100_median_pmic_min3']:.4f}\n- Median perturbation lower-tail pMIC proxy (10th percentile): {summary['top100_median_perturb_p10_min3']:.4f}\n- Median cross-model perturbation lower-tail floor: {summary['top100_median_perturb_crossmodel_p10_floor']:.4f}\n- Median fraction retaining within 0.5 pMIC under perturbation: {summary['top100_median_retention_0p5']:.3f}\n- Median maximum evaluator disagreement SD: {summary['top100_median_ensemble_std_max']:.4f}\n\n## Method\nGenerator output was frozen before this audit. R02 does not retrain M2. A 12k high-potential, six-stratum frontier is stress-tested using deterministic position-wise alanine-like substitutions, physicochemical-contrast substitutions, and N/C-terminal truncations. The same frozen independent pMIC ensemble is used only as a measurement instrument for local perturbation response. Selection is failure-oriented: potency, perturbation lower-tail performance, local retention, and evaluator consensus are converted to within-frontier percentiles; the candidate's audit score is the minimum (maximin) across those axes. Sequence-only physicochemical descriptors are used as caution/severe risk flags, not as toxicity claims. Official novelty/reference constraints are applied afterward, followed by a diversity gate.\n\n## Interpretation boundary\nThis is a computational pre-experimental audit. It does not establish MIC, HC50, hemolysis, pharmacokinetics, or safety. Its purpose is to reduce model-fragile and physically extreme candidates before wet-lab selection.\n'''
    (out/"PRIB_AUDIT_REPORT.md").write_text(report)
    print("PRIB_SELECT=PASS",json.dumps(summary,sort_keys=True),flush=True)
if __name__=="__main__": main()
