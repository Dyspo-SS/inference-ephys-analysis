"""PSTHs and prespecified response windows, with equal animal weighting."""
import numpy as np
import pandas as pd
from .common import UNIT,ANIMAL,BIN,TRIAL,save
from .statistics import signflip,compare,summary,adjust

def profiles(d,c):
    base=d[(d[BIN[0]]>=c["baseline"][0])&(d[BIN[1]]<=c["baseline"][1])]
    b=base.groupby(UNIT).Rate_Hz.agg(["mean","std"]).reset_index().rename(columns={"mean":"Baseline_FR_Hz","std":"Baseline_SD_Hz"})
    u=d.groupby(UNIT+BIN,as_index=False).agg(FR_Hz=("Rate_Hz","mean"),N_Trials=("Trial_ID","nunique"))
    u=u.merge(b,on=UNIT,how="left",validate="many_to_one")
    u["Response_Z"]=(u.FR_Hz-u.Baseline_FR_Hz)/u.Baseline_SD_Hz.where(u.Baseline_SD_Hz>0)
    p=u.groupby(ANIMAL+BIN,as_index=False).agg(FR_Hz=("FR_Hz","mean"),Response_Z=("Response_Z","mean"),N_Units=("Unit_ID","size"),N_Z_Units=("Response_Z","count"))
    return u,p,b

def run(d,c,out):
    u,p,b=profiles(d,c)
    save(u,out,"plotdata/unit_psth.csv");save(p,out,"plotdata/animal_psth.csv");save(b,out,"tables/baseline_parameters.csv")
    save(b[b.Baseline_SD_Hz.isna()|b.Baseline_SD_Hz.le(0)],out,"audit/normalization_unavailable.csv")
    gg=p.groupby(["Group","Region","Stage",*BIN],as_index=False).agg(Mean_FR_Hz=("FR_Hz","mean"),SEM_FR_Hz=("FR_Hz","sem"),Mean_Z=("Response_Z","mean"),SEM_Z=("Response_Z","sem"),N_Animals=("Mouse_ID","size"),N_Z_Animals=("Response_Z","count"))
    save(gg,out,"plotdata/group_psth.csv")
    unit_rows=[]
    for stage,spec in c["stages"].items():
        q=d[d.Stage.eq(stage)]
        if q.empty: continue
        base=q[q[BIN[0]].ge(c["baseline"][0])&q[BIN[1]].le(c["baseline"][1])].groupby(UNIT+TRIAL).Rate_Hz.mean().rename("Baseline_FR_Hz").reset_index()
        for name,w in spec["windows"].items():
            lo,hi=w["interval"]
            x=q[q[BIN[0]].ge(lo)&q[BIN[1]].le(hi)].groupby(UNIT+TRIAL).Rate_Hz.mean().rename("Event_FR_Hz").reset_index().merge(base,on=UNIT+TRIAL,validate="one_to_one")
            x["Delta_FR_Hz"]=x.Event_FR_Hz-x.Baseline_FR_Hz
            v=x.groupby(UNIT,as_index=False).agg(Event_FR_Hz=("Event_FR_Hz","mean"),Baseline_FR_Hz=("Baseline_FR_Hz","mean"),Delta_FR_Hz=("Delta_FR_Hz","mean"),N_Trials=("Trial_ID","nunique"))
            v["Window"]=name;v["Window_Start_s"]=lo;v["Window_End_s"]=hi;unit_rows.append(v)
    ur=pd.concat(unit_rows,ignore_index=True)
    ar=ur.groupby(ANIMAL+["Window","Window_Start_s","Window_End_s"],as_index=False).agg(Event_FR_Hz=("Event_FR_Hz","mean"),Baseline_FR_Hz=("Baseline_FR_Hz","mean"),Delta_FR_Hz=("Delta_FR_Hz","mean"),N_Units=("Unit_ID","size"),Min_Trials=("N_Trials","min"),Max_Trials=("N_Trials","max"))
    ar["Log_Rate_Ratio"]=np.nan
    ok=ar.Event_FR_Hz.gt(0)&ar.Baseline_FR_Hz.gt(0)
    ar.loc[ok,"Log_Rate_Ratio"]=np.log(ar.loc[ok,"Event_FR_Hz"]/ar.loc[ok,"Baseline_FR_Hz"])
    save(ur,out,"plotdata/unit_window_rates.csv");save(ar,out,"plotdata/animal_window_rates.csv")
    if c["study"]=="intervention":
        rows=[];groups=list(c["groups"]);window=c["intervention_window"]
        for family,metrics in [("rates_primary",["Event_FR_Hz","Log_Rate_Ratio"]),("rates_secondary",["Baseline_FR_Hz","Delta_FR_Hz"])]:
            for region in c["regions"]:
                q=ar[ar.Region.eq(region)&ar.Window.eq(window)]
                for metric in metrics:
                    row=dict(Family=family,Region=region,Window=window,Metric=metric,Group_A=groups[0],Group_B=groups[1])
                    row.update(compare(q[q.Group.eq(groups[0])][metric],q[q.Group.eq(groups[1])][metric],c["statistics"],family,region,metric))
                    rows.append(row)
        st=pd.DataFrame(rows);st=pd.concat([adjust(q) for _,q in st.groupby("Family",sort=False)],ignore_index=True)
    else:
        rows=[]
        for group,gs in c["groups"].items():
            for region in c["regions"]:
                for stage in gs["stages"]:
                    for name,w in c["stages"][stage]["windows"].items():
                        q=ar[ar.Group.eq(group)&ar.Region.eq(region)&ar.Window.eq(name)&ar.Stage.eq(stage)]
                        if w["family"]=="descriptive":
                            row=summary(q.Delta_FR_Hz,c["statistics"],group,region,name);row["P_Raw"]=np.nan
                        else: row=signflip(q.Delta_FR_Hz,c["statistics"],group,region,name)
                        row.update(Group=group,Region=region,Stage=stage,Window=name,Metric="Delta_FR_Hz",Family=w["family"],Mean_Baseline_Hz=q.Baseline_FR_Hz.mean(),Mean_Event_Hz=q.Event_FR_Hz.mean())
                        rows.append(row)
        st=pd.DataFrame(rows)
        st=pd.concat([adjust(q,"None" if fam in ["v1_primary","descriptive"] else "Holm") for fam,q in st.groupby("Family",sort=False)],ignore_index=True)
    save(st,out,"statistics/response_tests.csv")
    summaries=[]
    for ident,q in ar.groupby(["Group","Region","Stage","Window"]):
        for metric in ["Baseline_FR_Hz","Event_FR_Hz","Delta_FR_Hz","Log_Rate_Ratio"]:
            row=dict(zip(["Group","Region","Stage","Window"],ident));row["Metric"]=metric
            row.update(summary(q[metric],c["statistics"],*ident,metric));summaries.append(row)
    save(pd.DataFrame(summaries),out,"statistics/response_summaries.csv")
    return u,p,ar
