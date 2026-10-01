import torch
from model import AA

def decode(x):
    return "".join(AA[i] for i in x if 0 <= i < 20)

@torch.no_grad()
def generate_one(model, pmic, species_id=0, steps=20, min_len=8, max_len=50,
                 seed=1, fixed_length=None, device="cuda", delete_scale=0.35,
                 insert_scale=0.30, return_trace=False):
    generator=torch.Generator(device=device); generator.manual_seed(seed)
    initial_length=int(fixed_length) if fixed_length is not None else int(min_len)
    x=torch.randint(0,20,(1,initial_length),generator=generator,device=device)
    trace_steps=[]; total_delete=0; total_insert=0; total_floor_fill=0
    for step_index in range(steps):
        t=1.0-step_index/max(steps-1,1)
        mask=torch.ones_like(x,dtype=torch.bool)
        out=model(x,mask,torch.tensor([t],device=device),
                  torch.tensor([pmic],device=device),
                  torch.tensor([species_id],device=device))
        x=torch.multinomial(torch.softmax(out["residue"][0],-1),1,
                            generator=generator).squeeze(-1)[None,:]
        length_before_edit=int(x.shape[1])
        deleted_this_step=inserted_this_step=floor_fill_this_step=0
        if fixed_length is None:
            delete_probability=torch.sigmoid(out["delete"][0])*(float(delete_scale)*t)
            keep=torch.rand(delete_probability.shape,device=device,generator=generator)>delete_probability
            if int(keep.sum()) < min_len:
                keep[:min_len]=True
            deleted_this_step=int(x.shape[1])-int(keep.sum())
            x=x[:,keep]
            mask=torch.ones_like(x,dtype=torch.bool)
            out=model(x,mask,torch.tensor([t],device=device),
                      torch.tensor([pmic],device=device),
                      torch.tensor([species_id],device=device))
            insert_probability=torch.sigmoid(out["insert"][0])*(float(insert_scale)*t)
            insert_aa=torch.argmax(out["insert_aa"][0],-1)
            new_tokens=[]
            for gap in range(x.shape[1]+1):
                if len(new_tokens)<max_len and torch.rand((),device=device,generator=generator)<insert_probability[gap]:
                    new_tokens.append(int(insert_aa[gap])); inserted_this_step+=1
                if gap<x.shape[1] and len(new_tokens)<max_len:
                    new_tokens.append(int(x[0,gap]))
            if len(new_tokens)<min_len:
                floor_fill_this_step=min_len-len(new_tokens)
                new_tokens += [int(torch.randint(0,20,(1,),generator=generator,device=device))
                               for _ in range(floor_fill_this_step)]
            x=torch.tensor(new_tokens[:max_len],device=device,dtype=torch.long)[None,:]
        total_delete+=deleted_this_step; total_insert+=inserted_this_step
        total_floor_fill+=floor_fill_this_step
        trace_steps.append({"step":int(step_index),"t":float(t),
                            "length_before_edit":int(length_before_edit),
                            "delete_events":int(deleted_this_step),
                            "insert_events":int(inserted_this_step),
                            "floor_fill_events":int(floor_fill_this_step),
                            "length_after_edit":int(x.shape[1])})
    sequence=decode(x[0].tolist())
    trace={"initial_length":int(initial_length),"final_length":len(sequence),
           "delete_count":int(total_delete),"insert_count":int(total_insert),
           "floor_fill_count":int(total_floor_fill),
           "net_cardinality_change":int(len(sequence)-initial_length),
           "delete_scale":float(delete_scale),"insert_scale":float(insert_scale),
           "steps":trace_steps}
    return (sequence,trace) if return_trace else sequence
