"""Publication panels from exported plotting tables, without embedded result values."""
from pathlib import Path
import string
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from .common import BIN,rng,write_json,save

def style(c):
    plt.rcParams.update({"font.family":c["plot"]["font"],"font.size":8,"axes.labelsize":8,
        "xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8,"legend.frameon":False,
        "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":.7,
        "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})

def read(out,name):
    path=Path(out)/name
    columns=pd.read_csv(path,nrows=0).columns
    ids=["Group","Mouse_ID","Session_ID","Unit_ID","Region","Stage","Trial_ID"]
    return pd.read_csv(path,dtype={k:str for k in ids if k in columns})
def color(c,group):return c["groups"][group]["color"]
def label(c,group):return c["groups"][group].get("label",group.replace("_"," "))
def tag(ax,letter):ax.text(-.14,1.12,letter,transform=ax.transAxes,weight="bold",fontsize=11,va="top")

def export(fig,out,name,c):
    dest=Path(out)/"figures";dest.mkdir(exist_ok=True)
    for ext in c["plot"]["formats"]:fig.savefig(dest/f"{name}.{ext}",dpi=c["plot"]["dpi"],facecolor="white")
    plt.close(fig)

def points(ax,q,x,metric,group,c):
    shapes=["o","s","^","v","D","P","X","<",">","h"]
    for i,row in enumerate(q.sort_values("Mouse_ID").itertuples(index=False)):
        v=getattr(row,metric)
        if not np.isfinite(v):continue
        jitter=rng(c["statistics"]["seed"],"jitter",group,row.Mouse_ID,x).uniform(-.12,.12)
        marker=c.get("_animal_markers",{}).get((group,str(row.Mouse_ID)),shapes[i%len(shapes)])
        ax.scatter(x+jitter,v,s=15,marker=marker,color=color(c,group),alpha=.8,zorder=3)
    vals=q[metric].dropna().to_numpy()
    if len(vals):
        r=rng(c["statistics"]["seed"],"plot_ci",group,metric,x)
        boot=vals[r.integers(0,len(vals),(c["statistics"]["bootstrap_draws"],len(vals)))].mean(1)
        lo,hi=np.quantile(boot,[(1-c["statistics"]["ci_level"])/2,(1+c["statistics"]["ci_level"])/2]);mean=vals.mean()
        ax.errorbar(x+.19,mean,yerr=[[max(0,mean-lo)],[max(0,hi-mean)]],fmt="D",color="black",ms=3,capsize=2,lw=.8)

def annotate_p(ax,p,txt="adjusted p"):
    ax.text(.98,.97,f"{txt} = {p:.3g}" if np.isfinite(p) else "not estimable",transform=ax.transAxes,ha="right",va="top",fontsize=8)

def event_shading(ax,c,stage):
    palette={"Sound":"#4C78A8","Light":"#D6B95D","Shock":"#B66A80"}
    for name,(a,b) in c["stages"][stage]["events"].items():
        ax.axvspan(a,b,color=palette[name],alpha=.09,lw=0)
    ax.axvline(0,color=".5",ls="--",lw=.6)

def response_panel(fig,spec,out,c,group,letter="B"):
    unit=read(out,"plotdata/unit_psth.csv");gp=read(out,"plotdata/group_psth.csv")
    stages=c["groups"][group]["stages"];regions=c["regions"]
    gs=spec.subgridspec(len(regions),len(stages),hspace=.7,wspace=.38)
    first=None
    for r,region in enumerate(regions):
        for s,stage in enumerate(stages):
            inner=gs[r,s].subgridspec(2,1,height_ratios=[1.2,1],hspace=.05)
            ax=fig.add_subplot(inner[0]);bottom=fig.add_subplot(inner[1],sharex=ax)
            q=unit[unit.Group.eq(group)&unit.Region.eq(region)&unit.Stage.eq(stage)]
            if q.empty:
                ax.text(.5,.5,"No recording",ha="center",transform=ax.transAxes);bottom.set_axis_off();continue
            heat=q.pivot(index=["Mouse_ID","Session_ID","Unit_ID"],columns=BIN[0],values="Response_Z").sort_index()
            lo,hi=c["stages"][stage]["epoch"];lim=c["plot"]["heatmap_z_limit"]
            cmap=plt.colormaps["RdBu_r"].copy();cmap.set_bad(".8")
            image=ax.imshow(heat.to_numpy(),aspect="auto",interpolation="nearest",extent=[lo,hi,len(heat),0],cmap=cmap,vmin=-lim,vmax=lim,rasterized=True)
            sizes=heat.groupby(level="Mouse_ID",sort=False).size();bounds=sizes.cumsum()
            for b in bounds.iloc[:-1]:ax.axhline(b,color=".4",lw=.45)
            ax.set_yticks([.5,len(heat)-.5],["1",str(len(heat))]);ax.tick_params(axis="x",bottom=False,labelbottom=False)
            if s==0:ax.set_ylabel(f"{region}\nUnits")
            if r==0:ax.text(.5,1.10,stage,transform=ax.transAxes,ha="center")
            qg=gp[gp.Group.eq(group)&gp.Region.eq(region)&gp.Stage.eq(stage)]
            t=qg[BIN[0]].to_numpy()+c["bin_width_s"]/2;y=qg.Mean_Z.to_numpy();se=qg.SEM_Z.to_numpy()
            event_shading(bottom,c,stage);bottom.plot(t,y,color="black",lw=1);bottom.fill_between(t,y-se,y+se,color=".6",alpha=.25,lw=0)
            bottom.axhline(0,color=".6",lw=.5,ls=":");bottom.set_xlim(lo,hi);bottom.set_xlabel("Time from first cue (s)")
            if s==0:bottom.set_ylabel("Response (z)")
            if s==len(stages)-1:
                cb=fig.colorbar(image,ax=ax,fraction=.045,pad=.025);cb.set_ticks([-lim,0,lim]);cb.set_label("Baseline z")
            if first is None:first=ax
    if first is not None:
        tag(first,letter);first.text(0,1.3,label(c,group),transform=first.transAxes,color=color(c,group))

def glm_panel(fig,spec,out,c,letter="C"):
    a=read(out,"plotdata/glm_animal_gains.csv");st=read(out,"statistics/glm_event_locking.csv")
    regions=c["regions"];groups=list(c["groups"]);gs=spec.subgridspec(1,len(regions),wspace=.45)
    for j,reg in enumerate(regions):
        ax=fig.add_subplot(gs[j]);stage_list=[s for s in c["stages"] if s in c["glm"]["stages"]]
        for i,stage in enumerate(stage_list):
            for k,g in enumerate(groups):
                x=i+(k-(len(groups)-1)/2)*.32
                q=a[a.Group.eq(g)&a.Region.eq(reg)&a.Stage.eq(stage)&a.Term.eq("All_Events")]
                if q.empty:continue
                points(ax,q,x,"Gain_Bits_Per_Spike",g,c)
        ax.axhline(0,color=".6",lw=.6,ls="--");ax.set_xticks(range(len(stage_list)),stage_list,rotation=15);ax.set_ylabel("Gain (bits/spike)");ax.text(.02,.97,reg,transform=ax.transAxes,va="top")
        # Exact corrected probabilities remain in the source statistical table.
        if j==0:
            tag(ax,letter)
            handles=[Line2D([],[],marker="o",ls="",color=color(c,g),label=label(c,g)) for g in groups]
            ax.legend(handles=handles,loc="upper right",fontsize=7)
    return st

def rsa_panel(fig,spec,out,c,letter="D"):
    edges=read(out,"plotdata/rsa_matrices.csv");stats=read(out,"statistics/rsa_primary.csv")
    groups=list(c["groups"]);conds=c["rsa"]["primary_edge"] if c["study"]=="association" else list(c["rsa"]["conditions"])
    gs=spec.subgridspec(len(c["regions"]),len(groups)+1,width_ratios=[1]*len(groups)+[1.5],wspace=.5,hspace=.65)
    first=None
    for i,reg in enumerate(c["regions"]):
        for j,g in enumerate(groups):
            ax=fig.add_subplot(gs[i,j]);q=edges[edges.Group.eq(g)&edges.Region.eq(reg)]
            mat=q.groupby(["Condition_1","Condition_2"]).Similarity_r.mean().unstack().reindex(index=conds,columns=conds) if len(q) else pd.DataFrame(np.nan,index=conds,columns=conds)
            image=ax.imshow(mat,cmap="RdBu_r",vmin=-1,vmax=1,interpolation="nearest")
            ax.set_xticks(range(len(conds)),conds);ax.set_yticks(range(len(conds)),conds)
            for ii in range(len(conds)):
                for jj in range(len(conds)):
                    val=mat.iloc[ii,jj]
                    ax.text(jj,ii,f"{val:.2f}" if np.isfinite(val) else "NA",ha="center",va="center",color="white" if abs(val)>.6 else "black")
            if i==0:ax.text(.5,1.15,label(c,g),ha="center",transform=ax.transAxes,color=color(c,g))
            if j==0:ax.set_ylabel(reg)
            if first is None:first=ax
        ax=fig.add_subplot(gs[i,-1]);a,b=c["rsa"]["primary_edge"]
        for j,g in enumerate(groups):
            q=edges[edges.Group.eq(g)&edges.Region.eq(reg)&edges.Condition_1.eq(a)&edges.Condition_2.eq(b)]
            points(ax,q,j,"Similarity_r",g,c)
        ax.set_xticks(range(len(groups)),[label(c,g) for g in groups]);ax.set_ylabel(f"{a}–{b} similarity (r)");ax.set_ylim(-1.05,1.18);ax.axhline(0,color=".6",ls="--",lw=.6)
        s=stats[stats.Region.eq(reg)]
        if len(s):annotate_p(ax,s.P_Adjusted.iloc[0])
    if first is not None:tag(first,letter)

def behavior_panel(fig,spec,out,c,letter="A"):
    ax=fig.add_subplot(spec);t=read(out,"plotdata/behavior_groups.csv")
    for g,q in t.groupby("Group"):
        ax.errorbar(q.Day,q.Mean,yerr=q.SEM,marker="o",ms=3,lw=1,label=g.replace("_"," "))
    ax.set_xlabel("Day");ax.set_ylabel("Freezing (%)");ax.legend();tag(ax,letter)

def response_comparison(fig,spec,out,c,metric,letter):
    ar=read(out,"plotdata/animal_window_rates.csv");st=read(out,"statistics/response_tests.csv")
    gs=spec.subgridspec(1,len(c["regions"]),wspace=.55);groups=list(c["groups"])
    names={"Event_FR_Hz":"Sound firing rate (Hz)","Log_Rate_Ratio":"log(sound / baseline)","Baseline_FR_Hz":"Pre-sound rate (Hz)","Delta_FR_Hz":"Sound minus baseline (Hz)"}
    for i,reg in enumerate(c["regions"]):
        ax=fig.add_subplot(gs[i]);q=ar[ar.Region.eq(reg)&ar.Window.eq(c["intervention_window"])]
        for j,g in enumerate(groups):points(ax,q[q.Group.eq(g)],j,metric,g,c)
        ax.set_xticks(range(len(groups)),[label(c,g) for g in groups]);ax.set_ylabel(names[metric]);ax.text(.02,.97,reg,transform=ax.transAxes,va="top")
        q=st[st.Region.eq(reg)&st.Metric.eq(metric)]
        if len(q):annotate_p(ax,q.P_Adjusted.iloc[0])
        if i==0:tag(ax,letter)

def simple_grid(out,c,filename,table,metric,facets,xcol,ylim=None):
    df=read(out,table)
    if df.empty:return
    keys=list(df.groupby(facets,sort=False));cols=min(3,len(keys));rows=int(np.ceil(len(keys)/cols))
    fig,axs=plt.subplots(rows,cols,figsize=(190/25.4,max(2.5,2.3*rows)),squeeze=False,layout="constrained")
    for i,(ident,q) in enumerate(keys):
        ax=axs.flat[i];xvals=list(q[xcol].unique())
        if not isinstance(ident,tuple):ident=(ident,)
        for j,x in enumerate(xvals):
            for k,g in enumerate(c["groups"]):
                s=q[q[xcol].eq(x)&q.Group.eq(g)]
                if len(s):points(ax,s,j+(k-(len(c["groups"])-1)/2)*.32,metric,g,c)
        ax.set_xticks(range(len(xvals)),xvals,rotation=25,ha="right");ax.set_ylabel(metric.replace("_"," "));ax.text(.02,.97," · ".join(map(str,ident)),transform=ax.transAxes,va="top");ax.axhline(0,color=".6",ls="--",lw=.5)
        if ylim:ax.set_ylim(*ylim)
        tag(ax,string.ascii_uppercase[i%26])
    handles=[Line2D([],[],marker="o",ls="",color=color(c,g),label=label(c,g)) for g in c["groups"]]
    axs.flat[0].legend(handles=handles,loc="lower right",fontsize=7)
    for ax in list(axs.flat)[len(keys):]:ax.set_axis_off()
    export(fig,out,filename,c)

def quality(out,c):
    q=read(out,"tables/unit_count_qc.csv");b=read(out,"tables/baseline_parameters.csv");q=q.merge(b,on=["Group","Mouse_ID","Session_ID","Region","Unit_ID","Stage"])
    alloc=read(out,"tables/sample_allocation.csv")
    fig,axes=plt.subplots(len(c["regions"]),3,figsize=(190/25.4,2.6*len(c["regions"])),squeeze=False,layout="constrained")
    for i,reg in enumerate(c["regions"]):
        for col,metric in enumerate(["N_Unit_Sessions","Zero_Bin_Fraction","Baseline_FR_Hz"]):
            ax=axes[i,col];data=alloc if col==0 else q;data=data[data.Region.eq(reg)];stages=[s for s in c["stages"] if s in set(data.Stage)]
            for j,stage in enumerate(stages):
                for k,g in enumerate(c["groups"]):
                    z=data[data.Stage.eq(stage)&data.Group.eq(g)];x=j+(k-(len(c["groups"])-1)/2)*.3
                    if col:
                        # Unit distributions are descriptive; black diamonds show mouse means.
                        for mouse,m in z.groupby("Mouse_ID"):
                            jitter=rng(c["statistics"]["seed"],"qc",reg,metric,stage,g,mouse).uniform(-.10,.10,len(m));ax.scatter(x+jitter,m[metric],s=3,color=color(c,g),alpha=.28)
                        z=z.groupby(["Group","Mouse_ID"],as_index=False)[metric].mean()
                    points(ax,z,x,metric,g,c)
            ax.set_xticks(range(len(stages)),stages,rotation=20);ax.set_ylabel(f"{reg}\n"+metric.replace("_"," "));tag(ax,string.ascii_uppercase[i*3+col])
    handles=[Line2D([],[],marker="o",ls="",color=color(c,g),label=label(c,g)) for g in c["groups"]]
    axes[0,0].legend(handles=handles,loc="lower right",fontsize=7)
    export(fig,out,"Supplement_Count_QC",c)

def glm_supplements(out,c):
    null=read(out,"plotdata/glm_group_null.csv");st=read(out,"statistics/glm_event_locking.csv")
    keys=list(null.groupby(["Group","Region","Stage"]));cols=3;rows=int(np.ceil(len(keys)/cols))
    fig,axs=plt.subplots(rows,cols,figsize=(190/25.4,2.1*rows),squeeze=False,layout="constrained")
    for i,(ident,q) in enumerate(keys):
        ax=axs.flat[i];g,r,s=ident;v=st[st.Group.eq(g)&st.Region.eq(r)&st.Stage.eq(s)]
        ax.hist(q.Gain_Bits_Per_Spike,bins=30,color=".65");ax.axvline(v.Mean.iloc[0],color=color(c,g),lw=1.3);ax.set_xlabel("Gain (bits/spike)");ax.set_ylabel("Shift draws");ax.text(.02,.97,f"{label(c,g)} · {r} · {s}",transform=ax.transAxes,va="top",fontsize=8);tag(ax,string.ascii_uppercase[i%26])
    for ax in list(axs.flat)[len(keys):]:ax.set_axis_off()
    export(fig,out,"Supplement_GLM_Nulls",c)
    simple_grid(out,c,"Supplement_GLM_DropOne","plotdata/glm_animal_gains.csv","Gain_Bits_Per_Spike",["Region","Group"],"Term")
    pr=read(out,"plotdata/glm_animal_predictions.csv");keys=list(pr.groupby(["Region","Stage"]));cols=3;rows=int(np.ceil(len(keys)/cols))
    fig,axs=plt.subplots(rows,cols,figsize=(190/25.4,2.4*rows),squeeze=False,layout="constrained")
    for i,((reg,stage),q) in enumerate(keys):
        ax=axs.flat[i];event_shading(ax,c,stage)
        for group,v in q.groupby("Group"):
            p=v.groupby(BIN[0])[["Observed_FR_Hz","Predicted_Full_FR_Hz","Predicted_Reduced_FR_Hz"]].mean()
            for metric,ls in [("Observed_FR_Hz","-"),("Predicted_Full_FR_Hz","--"),("Predicted_Reduced_FR_Hz",":")]:ax.plot(p.index+c["bin_width_s"]/2,p[metric],ls=ls,color=color(c,group),lw=1,label=f"{label(c,group)} {metric.split('_')[0]} {metric.split('_')[1]}")
        ax.set_xlabel("Time (s)");ax.set_ylabel("Rate (Hz)");ax.text(.02,.97,f"{reg} · {stage}",transform=ax.transAxes,va="top");tag(ax,string.ascii_uppercase[i%26])
    axs.flat[0].legend(fontsize=7,loc="upper right")
    for ax in list(axs.flat)[len(keys):]:ax.set_axis_off()
    export(fig,out,"Supplement_GLM_Predictions",c)
    dg=read(out,"tables/glm_diagnostics.csv");fig,axs=plt.subplots(1,2,figsize=(190/25.4,2.8),layout="constrained")
    for group,q in dg.groupby("Group"):
        axs[0].hist(q.Train_Pearson_Dispersion.dropna(),bins=30,alpha=.4,label=label(c,group),color=color(c,group));axs[1].hist(q.Full_Iterations,bins=20,alpha=.4,color=color(c,group))
    for ax,name,t in zip(axs,["Training Pearson dispersion","Optimizer iterations"],["A","B"]):ax.set_xlabel(name);ax.set_ylabel("Unit-fold count");tag(ax,t)
    axs[0].legend();export(fig,out,"Supplement_GLM_Diagnostics",c)

def psth_comparison(out,c):
    p=read(out,"plotdata/group_psth.csv");stages=[s for s in c["stages"] if s in set(p.Stage)]
    fig,axs=plt.subplots(len(c["regions"]),len(stages),figsize=(190/25.4,2.4*len(c["regions"])),squeeze=False,layout="constrained")
    for i,reg in enumerate(c["regions"]):
        for j,stage in enumerate(stages):
            ax=axs[i,j];event_shading(ax,c,stage)
            for group,q in p[p.Region.eq(reg)&p.Stage.eq(stage)].groupby("Group"):
                t=q[BIN[0]].to_numpy()+c["bin_width_s"]/2;y=q.Mean_FR_Hz.to_numpy();se=q.SEM_FR_Hz.to_numpy();ax.plot(t,y,color=color(c,group),label=label(c,group));ax.fill_between(t,y-se,y+se,color=color(c,group),alpha=.15,lw=0)
            ax.set_xlabel("Time (s)");ax.set_ylabel(f"{reg} rate (Hz)");ax.text(.98,.97,stage,ha="right",va="top",transform=ax.transAxes);tag(ax,string.ascii_uppercase[i*len(stages)+j])
    axs[0,0].legend();export(fig,out,"Supplement_Raw_PSTH",c)

def rsa_profiles_plot(out,c):
    p=read(out,"plotdata/rsa_profiles.csv");p=p[p.Split.eq("All")]
    fig,axs=plt.subplots(1,len(c["regions"]),figsize=(190/25.4,3),squeeze=False,layout="constrained")
    palettes={"PS":"#4C78A8","CL":"#DBA242","TS":"#A44F81"}
    for ax,reg in zip(axs.flat,c["regions"]):
        for gi,group in enumerate(c["groups"]):
            for cond,q in p[p.Region.eq(reg)&p.Group.eq(group)].groupby("Condition"):
                s=q.groupby("Time_s").Value.agg(["mean","sem"]);ax.plot(s.index,s["mean"],ls="-" if gi==0 else "--",color=palettes.get(cond,".4"),label=f"{cond} {label(c,group)}")
        ax.set_xlabel("Time from cue (s)");ax.set_ylabel(f"{reg} response (z)")
    axs.flat[0].legend(fontsize=7);export(fig,out,"Supplement_RSA_Profiles",c)

def additional_matrices(out,c):
    edges=read(out,"plotdata/rsa_matrices.csv")
    groups=[g for g in c["groups"] if "Pairing" in c["groups"][g]["stages"]]
    if not groups:return
    conditions=list(c["rsa"]["conditions"])
    fig,axes=plt.subplots(len(c["regions"]),len(groups),figsize=(max(3,3*len(groups)),3*len(c["regions"])),squeeze=False,layout="constrained")
    for i,reg in enumerate(c["regions"]):
        for j,group in enumerate(groups):
            ax=axes[i,j];q=edges[edges.Region.eq(reg)&edges.Group.eq(group)]
            if q.empty:ax.set_axis_off();continue
            mat=q.groupby(["Condition_1","Condition_2"]).Similarity_r.mean().unstack().reindex(index=conditions,columns=conditions)
            ax.imshow(mat,vmin=-1,vmax=1,cmap="RdBu_r",interpolation="nearest")
            for ii in range(len(conditions)):
                for jj in range(len(conditions)):
                    v=mat.iloc[ii,jj];ax.text(jj,ii,f"{v:.2f}" if np.isfinite(v) else "NA",ha="center",va="center",color="white" if abs(v)>.6 else "black")
            ax.set_xticks(range(len(conditions)),conditions);ax.set_yticks(range(len(conditions)),conditions);ax.set_ylabel(reg);ax.text(.5,1.07,label(c,group),ha="center",transform=ax.transAxes);tag(ax,string.ascii_uppercase[(i*len(groups)+j)%26])
    export(fig,out,"Supplement_RSA_Available_Matrices",c)

def v1_main(out,c):
    u=read(out,"plotdata/unit_psth.csv");g=read(out,"plotdata/group_psth.csv");a=read(out,"plotdata/animal_window_rates.csv");st=read(out,"statistics/response_tests.csv")
    group=next(iter(c["groups"]));stage=c["groups"][group]["stages"][0]
    fig=plt.figure(figsize=(190/25.4,5.7));gs=fig.add_gridspec(2,2,left=.10,right=.96,bottom=.12,top=.93,hspace=.55,wspace=.5)
    ax=fig.add_subplot(gs[0,:]);heat=u.pivot(index=["Mouse_ID","Session_ID","Unit_ID"],columns=BIN[0],values="Response_Z").sort_index();lo,hi=c["stages"][stage]["epoch"];lim=c["plot"]["heatmap_z_limit"]
    im=ax.imshow(heat,aspect="auto",interpolation="nearest",extent=[lo,hi,len(heat),0],vmin=-lim,vmax=lim,cmap="RdBu_r",rasterized=True)
    sizes=heat.groupby(level="Mouse_ID",sort=False).size();bounds=sizes.cumsum()
    for b in bounds.iloc[:-1]:ax.axhline(b,color=".5",lw=.4)
    ax.set_yticks(bounds-sizes/2,sizes.index);ax.set_ylabel("V1 units by mouse");ax.set_xlabel("Time (s)");fig.colorbar(im,ax=ax,fraction=.025,pad=.02,label="Baseline z");tag(ax,"A")
    ax=fig.add_subplot(gs[1,0]);event_shading(ax,c,stage);t=g[BIN[0]].to_numpy()+c["bin_width_s"]/2;y=g.Mean_FR_Hz.to_numpy();e=g.SEM_FR_Hz.to_numpy();ax.plot(t,y,color="black");ax.fill_between(t,y-e,y+e,color=".6",alpha=.3);ax.set_xlabel("Time (s)");ax.set_ylabel("Firing rate (Hz)");tag(ax,"B")
    ax=fig.add_subplot(gs[1,1]);q=a[a.Window.eq(c["v1_primary_window"])]
    for _,row in q.iterrows():ax.plot([0,1],[row.Baseline_FR_Hz,row.Event_FR_Hz],color=".65",lw=.7)
    points(ax,q,0,"Baseline_FR_Hz",group,c);points(ax,q,1,"Event_FR_Hz",group,c);ax.set_xticks([0,1],["Pre-sound","Sound"]);ax.set_ylabel("Firing rate (Hz)");s=st[st.Window.eq(c["v1_primary_window"])];annotate_p(ax,s.P_Raw.iloc[0],"p");tag(ax,"C");export(fig,out,"V1_Main",c)

def run(out,c):
    # Local-only style metadata must not mutate the configuration saved in the run manifest.
    c=dict(c)
    animals=read(out,"plotdata/animal_psth.csv")[["Group","Mouse_ID"]].drop_duplicates()
    shapes=["o","s","^","v","D","P","X","<",">","h"]
    mapping=[];c["_animal_markers"]={}
    for group,q in animals.groupby("Group"):
        for i,mouse in enumerate(sorted(q.Mouse_ID.astype(str))):
            marker=shapes[i%len(shapes)];c["_animal_markers"][(group,mouse)]=marker
            mapping.append(dict(Group=group,Mouse_ID=mouse,Marker=marker,Color=color(c,group)))
    save(pd.DataFrame(mapping),out,"plotdata/animal_marker_key.csv")
    style(c)
    if c["study"]=="association":
        has_behavior=(Path(out)/"plotdata/behavior_groups.csv").exists()
        fig=plt.figure(figsize=(190/25.4,10.8 if has_behavior else 9.4));ratios=[1.1,3.2,1.4,2.5] if has_behavior else [3.2,1.4,2.5]
        gs=fig.add_gridspec(len(ratios),1,height_ratios=ratios,left=.10,right=.94,top=.96,bottom=.055,hspace=.7)
        offset=int(has_behavior)
        if has_behavior:behavior_panel(fig,gs[0],out,c)
        response_panel(fig,gs[offset],out,c,next(iter(c["groups"])),"B");glm_panel(fig,gs[offset+1],out,c,"C");rsa_panel(fig,gs[offset+2],out,c,"D");export(fig,out,"Fig2" if has_behavior else "Fig2_Electrophysiology_BCD",c)
        for group in list(c["groups"])[1:]:
            fig=plt.figure(figsize=(190/25.4,5));gs=fig.add_gridspec(1,1,left=.10,right=.94,bottom=.12,top=.9);response_panel(fig,gs[0],out,c,group,"A");export(fig,out,"Supplement_Control_Population",c)
    elif c["study"]=="intervention":
        fig=plt.figure(figsize=(190/25.4,9));gs=fig.add_gridspec(3,1,height_ratios=[1,1,2.4],left=.12,right=.94,top=.96,bottom=.065,hspace=.65)
        response_comparison(fig,gs[0],out,c,"Event_FR_Hz","J");response_comparison(fig,gs[1],out,c,"Log_Rate_Ratio","K");rsa_panel(fig,gs[2],out,c,"L");export(fig,out,"Fig4_Electrophysiology_JKL",c)
        fig=plt.figure(figsize=(190/25.4,5));gs=fig.add_gridspec(2,1,left=.13,right=.94,top=.94,bottom=.11,hspace=.6);response_comparison(fig,gs[0],out,c,"Baseline_FR_Hz","A");response_comparison(fig,gs[1],out,c,"Delta_FR_Hz","B");export(fig,out,"Supplement_Intervention_Rates",c)
    else:v1_main(out,c)
    quality(out,c);psth_comparison(out,c)
    simple_grid(out,c,"Supplement_Response_Windows","plotdata/animal_window_rates.csv","Delta_FR_Hz",["Region"],"Window")
    if c["glm"]["enabled"]:glm_supplements(out,c)
    if c["rsa"]["enabled"]:
        rsa_profiles_plot(out,c)
        if c["study"]=="association":additional_matrices(out,c)
        simple_grid(out,c,"Supplement_RSA_Reliability","plotdata/rsa_reliability.csv","Reliability_r",["Region"],"Condition",(-1.05,1.08))
        simple_grid(out,c,"Supplement_RSA_Sensitivity","plotdata/rsa_sensitivity.csv","Similarity_r",["Region"],"Variant",(-1.05,1.08))
        if (Path(out)/"plotdata/rsa_count_controls.csv").exists():simple_grid(out,c,"Supplement_RSA_CountControls","plotdata/rsa_count_controls.csv","Median_r",["Region"],"Variant",(-1.05,1.08))
    write_json(dict(Width_mm=190,Font=c["plot"]["font"],Heatmap_Row_Order="Mouse_ID, Session_ID, Unit_ID",Heatmap_Interpolation="nearest",PSTH_Error="animal SEM",Point_Error="animal percentile bootstrap CI",Formats=c["plot"]["formats"]),Path(out)/"audit/figure_settings.json")
