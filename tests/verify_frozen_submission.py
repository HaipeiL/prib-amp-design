from pathlib import Path
import hashlib

AA=set("ACDEFGHIKLMNPQRSTVWY")
ROOT=Path(__file__).resolve().parents[1]

def read_fasta(p):
    seqs=[]; cur=[]
    for line in p.read_text().splitlines():
        line=line.strip()
        if not line: continue
        if line.startswith(">"):
            if cur: seqs.append("".join(cur)); cur=[]
        else: cur.append(line.upper())
    if cur: seqs.append("".join(cur))
    return seqs

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()

lib=read_fasta(ROOT/"generate/library.fasta")
top=read_fasta(ROOT/"generate/top.fasta")
assert len(lib)==50000, len(lib)
assert len(set(lib))==50000
assert len(top)==100 and len(set(top))==100
assert set(top)<=set(lib)
assert all(8<=len(s)<=50 and set(s)<=AA for s in lib)
assert (ROOT/"generate/library.fasta").read_bytes()==(ROOT/"assets/frozen_library.fasta").read_bytes()
assert (ROOT/"generate/top.fasta").read_bytes()==(ROOT/"assets/frozen_top.fasta").read_bytes()
print("FROZEN_SUBMISSION=PASS")
print("library_sha256",digest(ROOT/"generate/library.fasta"))
print("top_sha256",digest(ROOT/"generate/top.fasta"))
