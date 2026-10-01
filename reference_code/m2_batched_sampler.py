import torch
from model import AA, PAD

def _pad(xs,device):
    B=len(xs); L=max(int(x.numel()) for x in xs)
    tok=torch.full((B,L),PAD,dtype=torch.long,device=device)
    mask=torch.zeros((B,L),dtype=torch.bool,device=device)
    for i,x in enumerate(xs):
        n=x.numel(); tok[i,:n]=x; mask[i,:n]=True
    return tok,mask

@torch.no_grad()
def generate_batch(model, pmic, species_id, seeds, steps=20, min_len=8, max_len=50,
                   device="cuda", delete_scale=0.35, insert_scale=0.30):
    B=len(seeds)
    assert len(pmic)==B and len(species_id)==B
    dev=torch.device(device)
    gens=[]
    xs=[]
    for seed in seeds:
        g=torch.Generator(device=dev); g.manual_seed(int(seed)); gens.append(g)
        xs.append(torch.randint(0,20,(min_len,),generator=g,device=dev))
    pmic_t=torch.tensor(pmic,dtype=torch.float32,device=dev)
    sid_t=torch.tensor(species_id,dtype=torch.long,device=dev)
    for step_index in range(steps):
        tval=1.0-step_index/max(steps-1,1)
        tok,mask=_pad(xs,dev)
        t=torch.full((B,),tval,dtype=torch.float32,device=dev)
        out=model(tok,mask,t,pmic_t,sid_t)
        res=[]
        for i,x in enumerate(xs):
            L=x.numel()
            y=torch.multinomial(torch.softmax(out["residue"][i,:L],-1),1,
                                generator=gens[i]).squeeze(-1)
            dp=torch.sigmoid(out["delete"][i,:L])*(float(delete_scale)*tval)
            keep=torch.rand(dp.shape,device=dev,generator=gens[i])>dp
            if int(keep.sum())<min_len:
                keep[:min_len]=True
            res.append(y[keep])
        xs=res
        tok,mask=_pad(xs,dev)
        out=model(tok,mask,t,pmic_t,sid_t)
        new=[]
        for i,x in enumerate(xs):
            L=x.numel()
            ip=torch.sigmoid(out["insert"][i,:L+1])*(float(insert_scale)*tval)
            ia=torch.argmax(out["insert_aa"][i,:L+1],-1)
            arr=[]
            for gap in range(L+1):
                if len(arr)<max_len and torch.rand((),device=dev,generator=gens[i])<ip[gap]:
                    arr.append(int(ia[gap]))
                if gap<L and len(arr)<max_len:
                    arr.append(int(x[gap]))
            while len(arr)<min_len:
                arr.append(int(torch.randint(0,20,(1,),generator=gens[i],device=dev)))
            new.append(torch.tensor(arr[:max_len],dtype=torch.long,device=dev))
        xs=new
    return ["".join(AA[i] for i in x.tolist()) for x in xs]
