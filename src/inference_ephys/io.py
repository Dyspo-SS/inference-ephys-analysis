"""Read an explicit study configuration and validate binned counts."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .common import UNIT, BIN, TRIAL, save, write_json, checksum

def read_config(path):
    path=Path(path).resolve()
    c=json.loads(path.read_text(encoding="utf-8"))
    for k in ["study","input","groups","regions","stages","baseline","statistics","glm","rsa"]:
        if k not in c: raise ValueError(f"Missing configuration field: {k}")
    if c["study"] not in ["association","intervention","v1"]:
        raise ValueError("study must be association, intervention, or v1")
    if c["study"]!="v1" and len(c["groups"])!=2: raise ValueError("Exactly two independent groups required")
    if c["study"]=="v1" and len(c["groups"])!=1: raise ValueError("V1 workflow requires one group")
    if c["invalid_bin_policy"] not in ["error","drop_unit_trial"]: raise ValueError("Invalid missing-bin policy")
    if c["bin_width_s"]<=0:raise ValueError("bin_width_s must be positive")
    if c["glm"].get("failure_policy") not in ["error","exclude"]:raise ValueError("GLM failure policy must be error or exclude")
    if not 0<c["statistics"]["ci_level"]<1:raise ValueError("CI level must be between zero and one")
    for k in ["bootstrap_draws","permutation_draws","exact_limit"]:
        if not isinstance(c["statistics"][k],int) or c["statistics"][k]<1:raise ValueError(f"Invalid {k}")
    if c["glm"]["enabled"] and c["glm"]["shift_draws"]<1:raise ValueError("At least one GLM shift is required")
    if c["glm"]["ridge_lambda"]<0:raise ValueError("Negative ridge penalty")
    if not 0<c["rsa"]["fisher_clip_epsilon"]<.01:raise ValueError("Invalid Fisher clipping epsilon")
    def on_grid(v,origin):return np.isclose((v-origin)/c["bin_width_s"],round((v-origin)/c["bin_width_s"]),atol=1e-7,rtol=0)
    for name,s in c["stages"].items():
        lo,hi=s["epoch"]
        if lo>=hi or c["baseline"][0]<lo or c["baseline"][1]>hi: raise ValueError(f"Invalid epoch/baseline for {name}")
        if not all(on_grid(v,lo) for v in [hi,*c["baseline"]]):raise ValueError(f"Epoch/baseline not aligned to bin grid: {name}")
        for event,window in s["events"].items():
            if event not in ["Sound","Light","Shock"] or not lo<=window[0]<window[1]<=hi:
                raise ValueError(f"Invalid event {name}/{event}")
            if not all(on_grid(v,lo) for v in window):raise ValueError(f"Event boundary inside bin: {name}/{event}")
        for w in s["windows"].values():
            a,b=w["interval"]
            if not lo<=a<b<=hi or not all(on_grid(v,lo) for v in [a,b]):raise ValueError(f"Invalid response window: {name}")
    for group,spec in c["groups"].items():
        if not spec["stages"] or set(spec["stages"])-set(c["stages"]):raise ValueError(f"Invalid stage selection for {group}")
    if c["rsa"]["enabled"]:
        if set(c["rsa"]["primary_edge"])-set(c["rsa"]["conditions"]):raise ValueError("Unknown primary RSA condition")
        for variant in [dict(window_s=c["rsa"]["window_s"],bin_width_s=c["bin_width_s"]),*c["rsa"]["variants"]]:
            lo,hi=variant["window_s"];width=variant["bin_width_s"]
            if width<=0 or hi<=lo or (hi-lo)/width<3 or not np.isclose(width/c["bin_width_s"],round(width/c["bin_width_s"])) or not np.isclose((hi-lo)/width,round((hi-lo)/width)):raise ValueError("RSA windows need at least three complete features aligned to the count grid")
            for condition,spec in c["rsa"]["conditions"].items():
                s=c["stages"][spec["stage"]];left,right=s["epoch"];a,b=lo+spec["onset_s"],hi+spec["onset_s"]
                if not left<=a<b<=right or not all(on_grid(v,left) for v in [a,b]):raise ValueError(f"RSA window outside/aligned incorrectly: {condition}")
    return c,path

def load_counts(c, config_path, out):
    path=(config_path.parent/c["input"]).resolve()
    ids=UNIT+["Trial_ID"]
    required=ids+["Trial_Order",*BIN,"Exposure_s","Spike_Count","Sound","Light","Shock","Valid_Bin"]
    d=pd.read_csv(path,dtype={k:"string" for k in ids})
    missing=set(required)-set(d)
    if missing: raise ValueError(f"Missing columns: {sorted(missing)}")
    if d.empty: raise ValueError("The input contains no observations")
    if d[ids].isna().any().any() or d[ids].apply(lambda v:v.str.strip().eq("").any()).any():
        raise ValueError("Empty identifiers")
    for k in set(required)-set(ids): d[k]=pd.to_numeric(d[k],errors="raise")
    if d.duplicated(UNIT+["Trial_ID",BIN[0]]).any(): raise ValueError("Duplicate unit/trial/bin keys")
    if not d.Valid_Bin.isin([0,1]).all(): raise ValueError("Valid_Bin must be 0 or 1")
    if not np.isfinite(d[["Trial_Order",*BIN,"Sound","Light","Shock"]]).all().all():
        raise ValueError("Nonfinite event/time/order values")
    if not ((d.Trial_Order>=1)&(d.Trial_Order%1==0)).all(): raise ValueError("Trial_Order must be a positive integer")
    if not d[["Sound","Light","Shock"]].isin([0,1]).all().all(): raise ValueError("Event indicators must be binary")
    if set(d.Group)-set(c["groups"]): raise ValueError("Unexpected groups; map labels in config or input")
    if set(d.Region)-set(c["regions"]): raise ValueError("Unexpected regions")
    if set(d.Stage)-set(c["stages"]): raise ValueError("Unexpected stage")
    original_rows=len(d);exclusions=[]
    def exclude(q,reason):
        rec=q[UNIT+TRIAL].drop_duplicates().copy();rec["Reason"]=reason;exclusions.append(rec)
    keep=pd.Series(False,index=d.index)
    for g,spec in c["groups"].items(): keep|=(d.Group==g)&d.Stage.isin(spec["stages"])
    exclude(d[~keep],"Stage outside configured recording selection")
    d=d[keep].copy()
    if d.empty: raise ValueError("No rows in the selected design")
    # A recording session is assigned one stage; units are keyed by session.
    if (d.groupby(["Group","Mouse_ID","Session_ID"]).Stage.nunique()>1).any():
        raise ValueError("A Session_ID cannot name more than one stage for an animal")
    if (d.groupby(["Group","Mouse_ID","Region","Stage"]).Session_ID.nunique()>1).any():
        raise ValueError("Select one session per mouse/region/stage before analysis")
    events=d[["Group","Mouse_ID","Session_ID",*TRIAL,*BIN,"Sound","Light","Shock"]].drop_duplicates()
    if events.duplicated(["Group","Mouse_ID","Session_ID","Trial_ID",BIN[0]]).any():
        raise ValueError("Event labels, timing, or trial order differ across units/regions")
    tr=d[["Group","Mouse_ID","Session_ID",*TRIAL]].drop_duplicates()
    if tr.duplicated(["Group","Mouse_ID","Session_ID","Trial_ID"]).any() or tr.duplicated(["Group","Mouse_ID","Session_ID","Trial_Order"]).any():
        raise ValueError("Trial_ID and Trial_Order must map one to one within session")
    bad_indices=[];bad_records=[];bw=c["bin_width_s"]
    for ident,q in d.groupby(UNIT+["Trial_ID"],sort=True):
        q=q.sort_values(BIN[0]);stage=ident[-2];spec=c["stages"][stage]
        expected=np.arange(spec["epoch"][0],spec["epoch"][1],bw)
        reason=None
        if (q.Valid_Bin!=1).any(): reason="Invalid bin in unit-trial epoch"
        elif len(q)!=len(expected) or not np.allclose(q[BIN[0]],expected,atol=1e-8,rtol=0): reason="Incomplete or unexpected time grid"
        elif not np.allclose(q[BIN[1]]-q[BIN[0]],bw,rtol=0,atol=1e-8): reason="Unequal bin widths"
        elif not np.isfinite(q[["Exposure_s","Spike_Count"]]).all().all(): reason="Missing valid-bin exposure/count"
        elif not np.allclose(q.Exposure_s,bw,atol=1e-8,rtol=0): reason="Partial or invalid exposure"
        elif (q.Spike_Count<0).any() or (q.Spike_Count%1!=0).any(): reason="Noninteger or negative counts"
        # Never drop corrupt counts or mislabeled events as artifact handling.
        if reason in ["Noninteger or negative counts"]: raise ValueError(f"{ident}: {reason}")
        for event in ["Sound","Light","Shock"]:
            window=spec["events"].get(event)
            target=((q[BIN[0]]>=window[0]-1e-8)&(q[BIN[1]]<=window[1]+1e-8)).astype(int) if window else np.zeros(len(q),int)
            if not np.array_equal(q[event],target): raise ValueError(f"{ident}: {event} labels disagree with configured event intervals")
        if reason:
            bad_indices.extend(q.index); rec=dict(zip(UNIT+["Trial_ID"],ident));rec["Trial_Order"]=int(q.Trial_Order.iloc[0]);rec["Reason"]=reason;bad_records.append(rec)
    if bad_records: exclusions.append(pd.DataFrame(bad_records))
    ex=pd.concat(exclusions,ignore_index=True) if exclusions else pd.DataFrame(columns=UNIT+TRIAL+["Reason"])
    save(ex,out,"audit/exclusions.csv")
    if bad_records and c["invalid_bin_policy"]=="error":
        raise ValueError("Invalid/incomplete epochs. See audit/exclusions.csv; no zero imputation is performed.")
    d=d.drop(index=bad_indices).copy()
    if d.empty: raise ValueError("No valid complete unit-trials remain")
    observed=set(d[["Group","Region","Stage"]].itertuples(index=False,name=None))
    required={(g,r,s) for g,v in c["groups"].items() for r in c["regions"] for s in v["stages"]}
    if required-observed:raise ValueError(f"Entire configured group/region/stage absent after selection: {sorted(required-observed)}")
    d["Spike_Count"]=d.Spike_Count.astype(np.int64)
    d["Rate_Hz"]=d.Spike_Count/d.Exposure_s
    d=d.sort_values(UNIT+TRIAL+[BIN[0]]).reset_index(drop=True)
    selected=d[UNIT].drop_duplicates()
    save(selected,out,"tables/unit_manifest.csv")
    save(d[["Group","Mouse_ID","Session_ID","Stage",*TRIAL,*BIN,"Sound","Light","Shock"]].drop_duplicates(),out,"tables/event_bins.csv")
    qc=d.groupby(UNIT,as_index=False).agg(N_Trials=("Trial_ID","nunique"),N_Bins=("Spike_Count","size"),Total_Spikes=("Spike_Count","sum"),Exposure_Total_s=("Exposure_s","sum"),Zero_Bin_Fraction=("Spike_Count",lambda v:(v==0).mean()),Max_Count=("Spike_Count","max"))
    qc["Mean_FR_Hz"]=qc.Total_Spikes/qc.Exposure_Total_s
    save(qc,out,"tables/unit_count_qc.csv")
    alloc=qc.groupby(["Group","Mouse_ID","Session_ID","Region","Stage"],as_index=False).agg(N_Unit_Sessions=("Unit_ID","size"),Min_Trials=("N_Trials","min"),Max_Trials=("N_Trials","max"),N_Bins=("N_Bins","sum"),Total_Spikes=("Total_Spikes","sum"))
    save(alloc,out,"tables/sample_allocation.csv")
    total=dict(Input_SHA256=checksum(path),Input_Rows=original_rows,Retained_Rows=len(d),
               N_Animals=len(d[["Group","Mouse_ID"]].drop_duplicates()),N_Sessions=len(d[["Group","Mouse_ID","Session_ID"]].drop_duplicates()),
               N_Unit_Sessions=len(selected),N_Trials=len(d[["Group","Mouse_ID","Session_ID","Trial_ID"]].drop_duplicates()),
               Note="Regional session counts overlap when both regions were recorded together; channel and waveform QC are not recoverable from binned counts.")
    write_json(total,Path(out)/"audit/input_summary.json")
    return d,path
