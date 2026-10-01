from pathlib import Path
import hashlib, json, math
from collections import Counter

AA=set("ACDEFGHIKLMNPQRSTVWY")
HYDRO=set("AVILMFWY")
EISENBERG={
    "A":0.62,"C":0.29,"D":-0.90,"E":-0.74,"F":1.19,"G":0.48,"H":-0.40,
    "I":1.38,"K":-1.50,"L":1.06,"M":0.64,"N":-0.78,"P":0.12,"Q":-0.85,
    "R":-2.53,"S":-0.18,"T":-0.05,"V":1.08,"W":0.81,"Y":0.26,
}

def sha256_file(path, chunk=1024*1024):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        while True:
            b=f.read(chunk)
            if not b: break
            h.update(b)
    return h.hexdigest()

def atomic_json(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+".tmp")
    tmp.write_text(json.dumps(obj,indent=2,sort_keys=True)+"\n")
    tmp.replace(path)

def read_fasta(path):
    seqs=[]; cur=[]
    for line in Path(path).read_text().splitlines():
        line=line.strip()
        if not line: continue
        if line.startswith(">"):
            if cur: seqs.append("".join(cur).upper()); cur=[]
        else: cur.append(line)
    if cur: seqs.append("".join(cur).upper())
    return seqs

def write_fasta(seqs,path,prefix="seq"):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with open(path,"w") as f:
        for i,s in enumerate(seqs,1): f.write(f">{prefix}{i:05d}\n{s}\n")

def valid_seq(s,min_len=8,max_len=50):
    return min_len<=len(s)<=max_len and not (set(s)-AA)

def rank01(vals, higher_is_better=True):
    vals=list(vals); n=len(vals)
    if n==0: return []
    order=sorted(range(n),key=lambda i:(vals[i],i),reverse=not higher_is_better)
    out=[0.0]*n; den=max(1,n-1)
    for r,i in enumerate(order): out[i]=r/den
    return out

def net_charge_proxy(s):
    return s.count("K")+s.count("R")-s.count("D")-s.count("E")

def hydrophobic_fraction(s): return sum(c in HYDRO for c in s)/len(s)

def max_run_from_set(s,chars):
    best=cur=0
    for c in s:
        cur=cur+1 if c in chars else 0; best=max(best,cur)
    return best

def max_residue_run(s):
    if not s: return 0
    best=cur=1
    for i in range(1,len(s)):
        cur=cur+1 if s[i]==s[i-1] else 1; best=max(best,cur)
    return best

def shannon_entropy(s):
    n=len(s); c=Counter(s)
    return -sum((v/n)*math.log2(v/n) for v in c.values())

def hydrophobic_moment_alpha(s):
    # Ideal alpha helix: 100 degrees/residue. Sequence-only proxy, not a structure claim.
    theta=math.radians(100.0); x=y=0.0
    for i,a in enumerate(s):
        h=EISENBERG[a]; x+=h*math.cos(i*theta); y+=h*math.sin(i*theta)
    return math.sqrt(x*x+y*y)/len(s)

def descriptors(s):
    q=net_charge_proxy(s); L=len(s); hf=hydrophobic_fraction(s)
    return {
      "length":L,"charge_proxy":q,"charge_density":q/L,"hydrophobic_fraction":hf,
      "max_hydrophobic_run":max_run_from_set(s,HYDRO),"max_residue_run":max_residue_run(s),
      "sequence_entropy":shannon_entropy(s),"hydrophobic_moment_alpha":hydrophobic_moment_alpha(s),
      "cysteine_count":s.count("C"),"methionine_count":s.count("M")
    }

def physics_flags(d,cfg):
    p=cfg["physics"]; caution=[]; severe=[]
    if d["hydrophobic_fraction"]>p["hydrophobic_fraction_caution_high"]: caution.append("high_hydrophobic_fraction")
    if d["hydrophobic_fraction"]>p["hydrophobic_fraction_severe_high"]: severe.append("extreme_hydrophobic_fraction")
    if d["charge_density"]>p["charge_density_caution_high"]: caution.append("high_charge_density")
    if d["charge_density"]>p["charge_density_severe_high"]: severe.append("extreme_charge_density")
    if d["max_hydrophobic_run"]>=p["max_hydrophobic_run_caution"]: caution.append("hydrophobic_run")
    if d["max_hydrophobic_run"]>=p["max_hydrophobic_run_severe"]: severe.append("extreme_hydrophobic_run")
    if d["sequence_entropy"]<p["entropy_caution_low"]: caution.append("low_complexity")
    if d["sequence_entropy"]<p["entropy_severe_low"]: severe.append("extreme_low_complexity")
    if d["cysteine_count"]>=p["cysteine_caution_count"]: caution.append("multiple_cysteines")
    if d["cysteine_count"]>=p["cysteine_severe_count"]: severe.append("many_cysteines")
    return caution,severe

def kmer_set(s,k=3): return {s[i:i+k] for i in range(max(0,len(s)-k+1))}

def jaccard_similarity(a,b):
    if not a and not b: return 1.0
    return len(a&b)/max(1,len(a|b))
