"""Similarity of population time courses on a shared temporal axis."""
import itertools
import numpy as np
import pandas as pd
from .common import UNIT,ANIMAL,BIN,save,correlation,rng
from .rates import profiles
from .statistics import compare,summary,adjust

def vectors(unit,c,window,metric="Response_Z",width=None):
    width=width or c["bin_width_s"]; rows=[]
    n=int(round((window[1]-window[0])/width))
    for condition,definition in c["rsa"]["conditions"].items():
        stage,onset=definition["stage"],definition["onset_s"]
        q=unit[unit.Stage.eq(stage)&unit[BIN[0]].ge(onset+window[0]-1e-8)&unit[BIN[1]].le(onset+window[1]+1e-8)].copy()
        if q.empty:continue
        q["Feature"]=np.floor((q[BIN[0]]-onset-window[0]+1e-8)/width).astype(int)
        # For the wider-bin sensitivity analysis, average already normalized bins.
        u=q.groupby(UNIT+["Feature"],as_index=False).agg(Value=(metric,"mean"),N_Source_Bins=(metric,"count"),N_Trials=("N_Trials","min"))
        needed=int(round(width/c["bin_width_s"]))
        u.loc[u.N_Source_Bins.ne(needed),"Value"]=np.nan
        a=u.groupby(ANIMAL+["Feature"],as_index=False).agg(Value=("Value","mean"),N_Valid_Units=("Value","count"),N_Total_Units=("Unit_ID","size"),Min_Trials=("N_Trials","min"),Max_Trials=("N_Trials","max"))
        a["Condition"]=condition;a["N_Features_Expected"]=n;a["Time_s"]=window[0]+(a.Feature+.5)*width
        rows.append(a)
    return pd.concat(rows,ignore_index=True) if rows else pd.DataFrame()

def matrix_edges(p):
    rows=[]
    for ident,q in p.groupby(["Group","Mouse_ID","Region"]):
        arrays={}
        for condition,v in q.groupby("Condition"):
            expected=int(v.N_Features_Expected.iloc[0]);arrays[condition]=v.set_index("Feature").Value.reindex(range(expected)).to_numpy()
        for a,b in itertools.product(arrays,repeat=2):
            r=correlation(arrays[a],arrays[b]);rows.append(dict(zip(["Group","Mouse_ID","Region"],ident))|dict(Condition_1=a,Condition_2=b,Similarity_r=r,N_Features=len(arrays[a]),Status="ok" if np.isfinite(r) else "undefined_profile"))
    return pd.DataFrame(rows)

def group_tests(edges,c,family):
    groups=list(c["groups"]);a,b=c["rsa"]["primary_edge"];rows=[]
    for region in c["regions"]:
        q=edges[edges.Region.eq(region)&edges.Condition_1.eq(a)&edges.Condition_2.eq(b)]
        av=q[q.Group.eq(groups[0])].Similarity_r.to_numpy();bv=q[q.Group.eq(groups[1])].Similarity_r.to_numpy()
        eps=c["rsa"]["fisher_clip_epsilon"]
        az=np.arctanh(np.clip(av,-1+eps,1-eps));bz=np.arctanh(np.clip(bv,-1+eps,1-eps))
        row=compare(az,bz,c["statistics"],family,region)
        row.update(Family=family,Region=region,Group_A=groups[0],Group_B=groups[1],Condition_1=a,Condition_2=b,Mean_r_A=np.nanmean(av) if np.isfinite(av).any() else np.nan,Mean_r_B=np.nanmean(bv) if np.isfinite(bv).any() else np.nan,Test_Scale="Fisher_z",N_Clipped_A=int(np.sum(np.abs(av)>1-eps)),N_Clipped_B=int(np.sum(np.abs(bv)>1-eps)))
        rows.append(row)
    return pd.DataFrame(rows)

def run(d,c,out,unit):
    allp=vectors(unit,c,c["rsa"]["window_s"]);allp["Split"]="All"
    splits=[allp]
    for label,parity in [("Odd",1),("Even",0)]:
        part=d[d.Trial_Order.mod(2).eq(parity)]
        if part.empty:continue
        u,_,_=profiles(part,c)
        p=vectors(u,c,c["rsa"]["window_s"]);p["Split"]=label;splits.append(p)
    p=pd.concat(splits,ignore_index=True);save(p,out,"plotdata/rsa_profiles.csv")
    edges=matrix_edges(allp);save(edges,out,"plotdata/rsa_matrices.csv")
    save(edges[edges.Status.ne("ok")],out,"audit/rsa_unavailable.csv")
    rel=[]
    for ident,q in p.groupby(["Group","Mouse_ID","Region","Condition"]):
        expected=int(q.N_Features_Expected.iloc[0]);arr={};counts={}
        for label in ["Odd","Even"]:
            s=q[q.Split.eq(label)];arr[label]=s.set_index("Feature").Value.reindex(range(expected)).to_numpy()
            counts[f"{label}_Min_Trials"]=s.Min_Trials.min();counts[f"{label}_Max_Trials"]=s.Max_Trials.max();counts[f"{label}_Min_Valid_Units"]=s.N_Valid_Units.min()
        rel.append(dict(zip(["Group","Mouse_ID","Region","Condition"],ident))|counts|dict(Reliability_r=correlation(arr["Odd"],arr["Even"])))
    save(pd.DataFrame(rel),out,"plotdata/rsa_reliability.csv")
    save(adjust(group_tests(edges,c,"rsa_primary")),out,"statistics/rsa_primary.csv")
    desc=[]
    for ident,q in edges.groupby(["Group","Region","Condition_1","Condition_2"]):
        desc.append(dict(zip(["Group","Region","Condition_1","Condition_2"],ident))|summary(q.Similarity_r,c["statistics"],"rsa_description",*ident))
    save(pd.DataFrame(desc),out,"statistics/rsa_edge_summaries.csv")
    sensitivity=[];tests=[]
    for variant in c["rsa"]["variants"]:
        vp=vectors(unit,c,variant["window_s"],variant["metric"],variant["bin_width_s"])
        ed=matrix_edges(vp);a,b=c["rsa"]["primary_edge"]
        ed=ed[ed.Condition_1.eq(a)&ed.Condition_2.eq(b)].copy();ed["Variant"]=variant["name"];sensitivity.append(ed)
        if c["study"]=="association":
            t=group_tests(ed,c,"rsa_sensitivity");t["Variant"]=variant["name"];tests.append(t)
    save(pd.concat(sensitivity,ignore_index=True),out,"plotdata/rsa_sensitivity.csv")
    if tests:save(adjust(pd.concat(tests,ignore_index=True)),out,"statistics/rsa_sensitivity.csv")
    return edges

def count_controls(d,c,out,animal_rates):
    """Descriptive unit subsampling and control-count thinning; no reliability matching."""
    opt=c["rsa"]["count_controls"]
    if not opt["enabled"]:return
    n=opt["units_per_stage"]; groups=list(c["groups"]); control,comparison=groups
    a,b=c["rsa"]["primary_edge"];ca,cb=c["rsa"]["conditions"][a],c["rsa"]["conditions"][b]
    probabilities={}
    for reg in c["regions"]:
        q=animal_rates[animal_rates.Region.eq(reg)&animal_rates.Window.eq(c["intervention_window"])]
        av=q[q.Group.eq(control)].Event_FR_Hz.mean();bv=q[q.Group.eq(comparison)].Event_FR_Hz.mean()
        probabilities[reg]=min(1.,bv/av) if np.isfinite(av) and np.isfinite(bv) and av>0 and bv>=0 else np.nan
    caches={}
    for ident,q in d[d.Stage.isin([ca["stage"],cb["stage"]])].groupby(ANIMAL):
        units=[]
        for uid,v in q.groupby("Unit_ID"):
            pivot=v.pivot(index="Trial_Order",columns=BIN[0],values="Spike_Count").sort_index(axis=1)
            units.append((str(uid),pivot.to_numpy(dtype=np.int64),pivot.columns.to_numpy(float)))
        caches[ident]=units
    def population(items,condition,keep,random):
        values=[];valid=0
        for _,counts,t in items:
            y=random.binomial(counts,keep) if keep<1 else counts
            rates=y/c["bin_width_s"]
            baseline=rates[:,(t>=c["baseline"][0])&(t<c["baseline"][1])]
            sd=baseline.std(ddof=1);mu=baseline.mean()
            if not np.isfinite(sd) or sd<=0:continue
            win=c["rsa"]["window_s"];on=condition["onset_s"]
            v=(rates[:,(t>=on+win[0])&(t<on+win[1])].mean(0)-mu)/sd
            values.append(v);valid+=1
        return (np.mean(values,axis=0) if values else np.full(int(round((c["rsa"]["window_s"][1]-c["rsa"]["window_s"][0])/c["bin_width_s"])),np.nan)),valid
    rows=[];skips=[]
    animals=d[["Group","Mouse_ID","Region"]].drop_duplicates()
    for group,mouse,reg in animals.itertuples(index=False,name=None):
        u=caches.get((group,mouse,reg,ca["stage"]),[]);v=caches.get((group,mouse,reg,cb["stage"]),[])
        if min(len(u),len(v))<n:
            skips.append(dict(Group=group,Mouse_ID=mouse,Region=reg,Reason="Too few units for configured subsampling",N_CL=len(u),N_TS=len(v)));continue
        for rep in range(opt["repeats"]):
            r=rng(c["statistics"]["seed"],"count_control",group,mouse,reg,rep)
            x=[u[i] for i in r.choice(len(u),n,replace=False)];y=[v[i] for i in r.choice(len(v),n,replace=False)]
            pv,nv=population(x,ca,1.,r)
            for variant,keep in [("unit_subsample",1.),("unit_subsample_and_thinning",probabilities[reg] if group==control else 1.)]:
                if not np.isfinite(keep):continue
                tv,nt=population(y,cb,keep,r)
                rows.append(dict(Group=group,Mouse_ID=mouse,Region=reg,Repeat=rep+1,Variant=variant,Keep_Probability=keep,N_Selected_Units=n,N_CL_Valid_Units=nv,N_TS_Valid_Units=nt,Similarity_r=correlation(pv,tv)))
    md=pd.DataFrame(rows)
    save(pd.DataFrame(skips,columns=["Group","Mouse_ID","Region","Reason","N_CL","N_TS"]),out,"audit/count_control_exclusions.csv")
    save(md if len(md) else pd.DataFrame(columns=["Group","Mouse_ID","Region","Repeat","Variant","Keep_Probability","Similarity_r"]),out,"plotdata/rsa_count_control_repeats.csv")
    if md.empty:return
    am=md.groupby(["Group","Mouse_ID","Region","Variant","Keep_Probability"],as_index=False).agg(Median_r=("Similarity_r","median"),Monte_Carlo_Lower=("Similarity_r",lambda x:x.quantile(.025)),Monte_Carlo_Upper=("Similarity_r",lambda x:x.quantile(.975)),N_Repeats=("Repeat","nunique"),N_Finite=("Similarity_r","count"))
    save(am,out,"plotdata/rsa_count_controls.csv")
    rows=[]
    for ident,q in am.groupby(["Group","Region","Variant"]):rows.append(dict(zip(["Group","Region","Variant"],ident))|summary(q.Median_r,c["statistics"],"matching",*ident))
    save(pd.DataFrame(rows),out,"statistics/rsa_count_controls_descriptive.csv")
