import torch
from torch import nn

AA = "ACDEFGHIKLMNPQRSTVWY"
PAD = 20
MASK = 21
VOCAB = 22

class AdaptiveEditTransformer(nn.Module):
    def __init__(self, d_model=768, nhead=12, layers=12, ff=3072,
                 max_len=64, species_buckets=256, dropout=0.1,
                 condition_content=True, condition_cardinality=True):
        super().__init__()
        self.d_model=int(d_model); self.nhead=int(nhead); self.layers=int(layers)
        self.ff=int(ff); self.max_len=int(max_len); self.species_buckets=int(species_buckets)
        self.dropout=float(dropout)
        self.condition_content=bool(condition_content)
        self.condition_cardinality=bool(condition_cardinality)
        self.tok=nn.Embedding(VOCAB,self.d_model,padding_idx=PAD)
        self.pos=nn.Embedding(self.max_len+2,self.d_model)
        self.time=nn.Sequential(nn.Linear(1,self.d_model),nn.SiLU(),nn.Linear(self.d_model,self.d_model))
        self.pmic=nn.Sequential(nn.Linear(1,self.d_model),nn.SiLU(),nn.Linear(self.d_model,self.d_model))
        self.species=nn.Embedding(self.species_buckets,self.d_model)
        layer=nn.TransformerEncoderLayer(self.d_model,self.nhead,self.ff,self.dropout,
                                         batch_first=True,norm_first=True,activation="gelu")
        self.enc=nn.TransformerEncoder(layer,self.layers)
        self.norm=nn.LayerNorm(self.d_model)
        self.residue=nn.Linear(self.d_model,20)
        self.delete=nn.Linear(self.d_model,1)
        self.insert=nn.Linear(self.d_model,1)
        self.insert_aa=nn.Linear(self.d_model,20)

    def get_config(self):
        return {"d_model":self.d_model,"nhead":self.nhead,"layers":self.layers,"ff":self.ff,
                "max_len":self.max_len,"species_buckets":self.species_buckets,"dropout":self.dropout,
                "condition_content":self.condition_content,
                "condition_cardinality":self.condition_cardinality}

    def forward(self,tokens,mask,t,pmic,species_id):
        B,L=tokens.shape
        pos=torch.arange(L,device=tokens.device)[None,:]
        h=self.tok(tokens)+self.pos(pos)+self.time(t[:,None]).unsqueeze(1)
        h=self.enc(h,src_key_padding_mask=~mask)
        h=self.norm(h)
        h=h.masked_fill(~mask.unsqueeze(-1),0.0)
        cond=self.pmic(pmic[:,None]).unsqueeze(1)+self.species(species_id).unsqueeze(1)
        cond_content=cond if self.condition_content else cond*0.0
        cond_cardinality=cond if self.condition_cardinality else cond*0.0
        residue_logits=self.residue(h+cond_content)
        delete_logits=self.delete(h+cond_cardinality).squeeze(-1)
        zero_gap=torch.zeros(B,1,self.d_model,device=h.device,dtype=h.dtype)
        gap_base=torch.cat([zero_gap,h],dim=1)
        insert_logits=self.insert(gap_base+cond_cardinality.expand(B,L+1,self.d_model)).squeeze(-1)
        insert_aa_logits=self.insert_aa(gap_base+cond_content.expand(B,L+1,self.d_model))
        return {"residue":residue_logits,"delete":delete_logits,
                "insert":insert_logits,"insert_aa":insert_aa_logits}
