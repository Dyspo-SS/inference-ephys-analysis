"""Optional session-level behavioral panel; animals are not linked to electrophysiology."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from .common import save
from .statistics import adjust

def run(c,config_path,out):
    opt=c.get("behavior",{})
    if not opt.get("input"):return
    d=pd.read_csv((config_path.parent/opt["input"]).resolve(),dtype={"Mouse_ID":"string","Group":"string"})
    needed=["Mouse_ID","Group","Day","Freezing_Percent"]
    if not set(needed)<=set(d):raise ValueError("Behavior columns required: "+", ".join(needed))
    if d[needed].isna().any().any() or d.duplicated(["Group","Mouse_ID","Day"]).any():raise ValueError("Missing or duplicate animal/day behavioral records")
    if not np.isfinite(d[["Day","Freezing_Percent"]]).all().all() or not d.Freezing_Percent.between(0,100).all():raise ValueError("Invalid freezing percentage/day")
    if set(d.Group)!=set(opt["groups"]) or set(d.Day)!=set(opt["days"]):raise ValueError("Behavior groups/days differ from planned comparison family")
    save(d[needed],out,"plotdata/behavior_animals.csv")
    desc=d.groupby(["Group","Day"],as_index=False).agg(Mean=("Freezing_Percent","mean"),SEM=("Freezing_Percent","sem"),N_Animals=("Mouse_ID","size"));save(desc,out,"plotdata/behavior_groups.csv")
    rows=[]
    for day in opt["days"]:
        for control in opt["controls"]:
            a=d[d.Day.eq(day)&d.Group.eq(opt["reference"])].Freezing_Percent.to_numpy()
            b=d[d.Day.eq(day)&d.Group.eq(control)].Freezing_Percent.to_numpy()
            result=ttest_ind(a,b,equal_var=False) if min(len(a),len(b))>=2 else None
            rows.append(dict(Day=day,Group_A=opt["reference"],Group_B=control,N_A=len(a),N_B=len(b),Difference_A_Minus_B=a.mean()-b.mean(),T=result.statistic if result else np.nan,DF=result.df if result else np.nan,P_Raw=result.pvalue if result else np.nan,Test="two_sided_Welch",Family="behavior_day_control"))
    save(adjust(pd.DataFrame(rows)),out,"statistics/behavior_tests.csv")
