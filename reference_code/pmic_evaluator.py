import torch, json, math
AA="ACDEFGHIKLMNPQRSTVWY"; AA_TO_ID={a:i for i,a in enumerate(AA)}; PAD=20; VOCAB=21

class IndependentPMICEvaluator(torch.nn.Module):
    def __init__(self,n_species,d_model=384,nhead=6,layers=8,ff=1536,max_len=100,dropout=0.15):
        super().__init__(); self.n_species=int(n_species); self.d_model=int(d_model); self.nhead=int(nhead)
        self.layers=int(layers); self.ff=int(ff); self.max_len=int(max_len); self.dropout=float(dropout)
        self.aa=torch.nn.Embedding(VOCAB,self.d_model,padding_idx=PAD); self.pos=torch.nn.Embedding(self.max_len,self.d_model)
        self.species=torch.nn.Embedding(self.n_species+1,self.d_model)
        block=torch.nn.TransformerEncoderLayer(d_model=self.d_model,nhead=self.nhead,dim_feedforward=self.ff,
                                               dropout=self.dropout,activation="gelu",batch_first=True,norm_first=True)
        self.encoder=torch.nn.TransformerEncoder(block,num_layers=self.layers); self.norm=torch.nn.LayerNorm(self.d_model)
        self.head=torch.nn.Sequential(torch.nn.Linear(self.d_model,self.d_model),torch.nn.GELU(),
                                      torch.nn.Dropout(self.dropout),torch.nn.Linear(self.d_model,1))
    def forward(self,tokens,mask,species_id):
        B,L=tokens.shape; pos=torch.arange(L,device=tokens.device)[None,:]
        h=self.aa(tokens)+self.pos(pos)+self.species(species_id)[:,None,:]
        h=self.encoder(h,src_key_padding_mask=~mask); h=self.norm(h); h=h.masked_fill(~mask.unsqueeze(-1),0.0)
        denom=mask.sum(dim=1,keepdim=True).clamp_min(1).to(h.dtype); pooled=h.sum(dim=1)/denom
        return self.head(pooled).squeeze(-1)

def _encode(seqs,species,map_,device):
    L=max(map(len,seqs)); tok=torch.full((len(seqs),L),PAD,dtype=torch.long,device=device)
    mask=torch.zeros((len(seqs),L),dtype=torch.bool,device=device); sid=torch.zeros(len(seqs),dtype=torch.long,device=device)
    for i,(s,sp) in enumerate(zip(seqs,species)):
        ids=[AA_TO_ID[x] for x in s]; tok[i,:len(ids)]=torch.tensor(ids,device=device); mask[i,:len(ids)]=True
        sid[i]=int(map_.get(sp,0))
    return tok,mask,sid

def load_ensemble(paths,device="cuda:0"):
    models=[]; refmap=None
    for p in paths:
        o=torch.load(p,map_location="cpu"); sm=o["species_map"]
        if refmap is None: refmap=sm
        elif sm!=refmap: raise RuntimeError("Evaluator species maps differ")
        m=IndependentPMICEvaluator(**o["model_config"]).to(device); m.load_state_dict(o["model"]); m.eval()
        models.append((m,float(o["target_mean"]),float(o["target_std"])))
    return models,refmap

@torch.no_grad()
def score_broad(seqs,paths,species3,batch=1024,device="cuda:0"):
    models,smap=load_ensemble(paths,device); out_mean=[]; out_min=[]
    for st in range(0,len(seqs),batch):
        ss=seqs[st:st+batch]; per=[]
        for sp in species3:
            tok,mask,sid=_encode(ss,[sp]*len(ss),smap,device); ens=[]
            for m,mu,sd in models:
                with torch.autocast("cuda",dtype=torch.bfloat16):
                    z=m(tok,mask,sid)
                ens.append(z.float()*sd+mu)
            per.append(torch.stack(ens).mean(0))
        mat=torch.stack(per,dim=1)
        out_mean += mat.mean(1).cpu().tolist(); out_min += mat.min(1).values.cpu().tolist()
    return out_mean,out_min
