import argparse,csv,json,math,statistics,sys
from pathlib import Path
from common import read_fasta,write_fasta,rank01,kmer_set,jaccard_similarity,atomic_json,sha256_file
try:
    import Levenshtein
except Exception as e:
    raise SystemExit('Levenshtein required; include R02 src on PYTHONPATH: '+str(e))

def load_csv(p):
    with open(p,newline='') as f:return list(csv.DictReader(f))

def best_ref_similarity(seq,refs,stop=0.80):
    best=0.0
    for r in refs:
        ub=2*min(len(seq),len(r))/(len(seq)+len(r))
        if ub<=best: continue
        v=Levenshtein.ratio(seq,r)
        if v>best:best=v
        if best>stop:break
    return best

def gmean(vals):
    vals=[max(1e-9,float(v)) for v in vals]
    return math.exp(sum(math.log(v) for v in vals)/len(vals))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--config',required=True); ap.add_argument('--run-dir',required=True)
    a=ap.parse_args(); cfg=json.loads(Path(a.config).read_text()); run=Path(a.run_dir); r02=Path(cfg['r02_run']); r01=Path(cfg['r01_run'])
    rows=load_csv(run/'work/lane0/panel_scores.csv')+load_csv(run/'work/lane1/panel_scores.csv')
    if len(rows)!=12000: raise RuntimeError(f'panel rows {len(rows)} != 12000')
    nums=['r02_reliability_floor','r02_reliability_geomean','broad_mean','broad_p10','gram_negative_mean','gram_positive_mean','mdr_mean','panel_species_floor','panel_crossmodel_floor','panel_max_model_sd']
    for r in rows:
        for k in nums:r[k]=float(r[k])
    axes={
      'p_broad':rank01([r['broad_p10'] for r in rows]),
      'p_gneg':rank01([r['gram_negative_mean'] for r in rows]),
      'p_gpos':rank01([r['gram_positive_mean'] for r in rows]),
      'p_mdr':rank01([r['mdr_mean'] for r in rows]),
      'p_panel_floor':rank01([r['panel_crossmodel_floor'] for r in rows]),
      'p_consensus':rank01([-r['panel_max_model_sd'] for r in rows]),
    }
    for i,r in enumerate(rows):
        for k,v in axes.items():r[k]=v[i]
        core=[r['r02_reliability_floor'],r['p_broad'],r['p_gneg'],r['p_gpos'],r['p_mdr'],r['p_panel_floor'],r['p_consensus']]
        r['panel_reliability_floor']=min(core); r['panel_reliability_geomean']=gmean(core)
        r['caution_count']=0 if not r.get('physics_caution') else len([x for x in r['physics_caution'].split(';') if x])
    rows.sort(key=lambda r:(r['r02_status']!='PASS', bool(r.get('physics_severe')), -r['panel_reliability_floor'], -r['panel_reliability_geomean'], r['caution_count'], -r['broad_p10'], int(r['candidate_id'])))
    night=json.loads((r01/'configs/night_config.json').read_text()); refs=read_fasta(night['official_reference'])
    eligible=[]
    for r in rows[:min(int(cfg['reference_scan_pool']),len(rows))]:
        if r['r02_status']!='PASS' or r.get('physics_severe'):continue
        sim=best_ref_similarity(r['sequence'],refs); r['max_reference_similarity']=sim
        if sim<=float(cfg['top_similarity_safety']):eligible.append(r)
    if len(eligible)<100:raise RuntimeError(f'eligible {len(eligible)} < 100')
    top=[]; km=[]
    def fill(gate):
        for r in eligible:
            if len(top)>=100:return
            if any(r['candidate_id']==x['candidate_id'] for x in top):continue
            k=kmer_set(r['sequence'])
            if km and max(jaccard_similarity(k,z) for z in km)>gate:continue
            top.append(r);km.append(k)
    fill(float(cfg['pairwise_kmer_similarity_gate'])); fill(float(cfg['pairwise_kmer_similarity_relaxed']))
    if len(top)<100:
        for r in eligible:
            if len(top)>=100:break
            if all(r['candidate_id']!=x['candidate_id'] for x in top):top.append(r)
    if len(top)!=100:raise RuntimeError('top100 fill failed')
    topseq=[r['sequence'] for r in top]; oldlib=read_fasta(r02/'results/prib_library.fasta'); seen=set();lib=[]
    for s in topseq+oldlib:
        if s not in seen:seen.add(s);lib.append(s)
        if len(lib)>=50000:break
    if len(lib)!=50000:raise RuntimeError(f'library size {len(lib)}')
    out=run/'results';out.mkdir(parents=True,exist_ok=True);write_fasta(lib,out/'panel_library.fasta','lib');write_fasta(topseq,out/'panel_top.fasta','top')
    oldtop=set(read_fasta(r02/'results/prib_top.fasta'));overlap=sum(s in oldtop for s in topseq)
    basecols=['rank','candidate_id','sequence','panel_reliability_floor','panel_reliability_geomean','r02_reliability_floor','broad_mean','broad_p10','gram_negative_mean','gram_positive_mean','mdr_mean','panel_species_floor','panel_crossmodel_floor','panel_max_model_sd','physics_caution','max_reference_similarity']
    spcols=[f"pmic_{sp.replace(' ','_')}" for sp in cfg['official_species']]
    with open(out/'panel_top100_audit.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=basecols+spcols);w.writeheader()
        for i,r in enumerate(top,1):
            q={k:r.get(k,'') for k in basecols+spcols};q['rank']=i;w.writerow(q)
    with open(out/'panel_frontier_scores.csv','w',newline='') as f:
        cols=list(rows[0].keys());w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
    sm={
      'status':'PASS','frontier_size':len(rows),'top_size':100,'library_size':50000,'r02_top100_overlap':overlap,'r02_top100_replaced':100-overlap,
      'top100_max_reference_similarity':max(r['max_reference_similarity'] for r in top),
      'top100_median_broad_p10':statistics.median(r['broad_p10'] for r in top),
      'top100_median_gram_negative_mean':statistics.median(r['gram_negative_mean'] for r in top),
      'top100_median_gram_positive_mean':statistics.median(r['gram_positive_mean'] for r in top),
      'top100_median_mdr_mean':statistics.median(r['mdr_mean'] for r in top),
      'top100_median_panel_crossmodel_floor':statistics.median(r['panel_crossmodel_floor'] for r in top),
      'top100_median_r02_reliability_floor':statistics.median(r['r02_reliability_floor'] for r in top),
      'interpretation':'Species-level frozen-evaluator proxy aligned to the 10 unique species underlying the official 20-strain panel. It is not strain-specific MIC prediction and does not model HC50.'
    }
    atomic_json(out/'panel_summary.json',sm)
    report=f'''# PRIB Pan-Panel Audit R04\n\nStatus: **PASS**\n\nR04 keeps the frozen M2 generator and R02 perturbation audit unchanged, then adds an independent species-level stress test aligned to the 10 unique bacterial species represented in the official 20-strain panel. The frozen three-model pMIC evaluator ensemble is queried for every candidate/species pair. Official strain multiplicities are used only as weights when summarizing Broad, Gram-negative, Gram-positive, and MDR activity proxies.\n\nThis is deliberately **not** presented as strain-specific prediction: the evaluator conditions on species labels, not ATCC/AIC/MDR strain identities. HC50/hemolysis is also not predicted. Physicochemical flags remain pre-experimental risk proxies only.\n\n- R02 Top100 retained: {overlap}/100\n- R02 Top100 replaced: {100-overlap}/100\n- Median broad-panel lower-tail proxy: {sm['top100_median_broad_p10']:.4f}\n- Median Gram-negative proxy: {sm['top100_median_gram_negative_mean']:.4f}\n- Median Gram-positive proxy: {sm['top100_median_gram_positive_mean']:.4f}\n- Median MDR proxy: {sm['top100_median_mdr_mean']:.4f}\n- Median cross-model/species floor: {sm['top100_median_panel_crossmodel_floor']:.4f}\n- Max official-reference similarity: {sm['top100_max_reference_similarity']:.6f}\n'''
    (out/'PANEL_AUDIT_REPORT.md').write_text(report)
    # manifest
    files=['panel_library.fasta','panel_top.fasta','panel_top100_audit.csv','panel_frontier_scores.csv','panel_summary.json','PANEL_AUDIT_REPORT.md']
    with open(out/'R04_manifest.sha256','w') as f:
        for name in files:f.write(f"{sha256_file(out/name)}  {name}\n")
    print('PANEL_SELECT=PASS',json.dumps(sm,sort_keys=True),flush=True)
if __name__=='__main__':main()
