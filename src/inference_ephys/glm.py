"""Trial-held-out ridge Poisson models and shared within-trial shift nulls."""
import numpy as np
import pandas as pd
from .common import UNIT,ANIMAL,BIN,save,rng
from .statistics import summary,adjust

def fit_batch(X,Y,exposure,opt):
    Y=np.atleast_2d(Y)
    if np.linalg.matrix_rank(X)!=X.shape[1]:raise ValueError("Rank-deficient design")
    if (Y.sum(1)<=0).any():raise ValueError("Zero training spikes; finite intercept MLE unavailable")
    penalty=np.full(X.shape[1],opt["ridge_lambda"]);penalty[0]=0
    beta=np.zeros((len(Y),X.shape[1]));beta[:,0]=np.log(Y.sum(1)/exposure.sum())
    for iteration in range(1,opt["max_iter"]+1):
        eta=beta@X.T
        if np.any(np.abs(eta)>100):raise ValueError("Extreme linear predictor")
        mu=np.exp(eta)*exposure
        grad=(Y-mu)@X-beta*penalty
        info=np.einsum("ni,bn,nj->bij",X,mu,X,optimize=True)+np.diag(penalty)
        step=np.linalg.solve(info,grad[...,None])[...,0]
        objective=(Y*eta-mu).sum(1)-.5*(beta*beta*penalty).sum(1)
        scale=np.ones(len(Y));accepted=False
        for _ in range(40):
            proposal=beta+step*scale[:,None];ep=proposal@X.T
            with np.errstate(over="ignore",invalid="ignore"):
                value=(Y*ep-np.exp(ep)*exposure).sum(1)-.5*(proposal*proposal*penalty).sum(1)
            bad=~np.isfinite(value)|(value<objective-1e-9)
            if not bad.any():accepted=True;break
            scale[bad]*=.5
        if not accepted:raise ValueError("Optimizer line search failed")
        change=step*scale[:,None];beta+=change
        # Both the parameter step and score must be small; tiny line-search steps alone do not imply convergence.
        eta=beta@X.T;gradient=(Y-np.exp(eta)*exposure)@X-beta*penalty
        if np.max(np.abs(change))<opt["tolerance"] and np.max(np.abs(gradient))<opt["gradient_tolerance"]:
            return beta,iteration
    raise ValueError("Poisson optimizer did not converge")

def run(d,c,out):
    opt=c["glm"];M=opt["shift_draws"]
    chosen=d[d.Stage.isin(opt["stages"])].copy()
    unit_rows=[];fold_rows=[];coef_rows=[];pred_rows=[];diagnostics=[];perm_rows=[];failures=[];split_rows=[]
    units=list(chosen.groupby(UNIT,sort=True))
    for number,(ident,g) in enumerate(units,1):
        meta=dict(zip(UNIT,ident));local={k:[] for k in ["units","folds","coefs","preds","diagnostics","perms","splits"]}
        try:
            terms=list(c["stages"][meta["Stage"]]["events"])
            if not terms:raise ValueError("No configured event terms")
            y=[];orders=[];trial_ids=[];reference=None
            for (order,tid),q in g.groupby(["Trial_Order","Trial_ID"],sort=True):
                q=q.sort_values(BIN[0]);a=q[BIN+["Exposure_s"]+terms].to_numpy(float)
                if reference is None:reference=a
                elif a.shape!=reference.shape or not np.allclose(a,reference,rtol=0,atol=1e-8):raise ValueError("Unequal trial timing for shift model")
                y.append(q.Spike_Count.to_numpy(float));orders.append(order);trial_ids.append(tid)
            y=np.asarray(y);orders=np.asarray(orders,float);T,B=y.shape
            if y.sum()<=0:raise ValueError("No spikes; bits/spike undefined")
            if meta["Stage"] in opt["leave_one_trial_out_stages"]:
                if T<3:raise ValueError("Leave-one-trial-out needs at least three trials")
                splits=[(f"Leave_{tid}_out",orders!=order) for tid,order in zip(trial_ids,orders)]
            else:
                if min(np.sum(orders%2==0),np.sum(orders%2==1))<2:raise ValueError("Need at least two trials in each odd/even partition")
                splits=[("Odd_to_Even",orders%2==1),("Even_to_Odd",orders%2==0)]
            categories,inv=np.unique(reference[:,3:],axis=0,return_inverse=True);C=len(categories)
            mask=np.eye(C)[inv];exposure=(reference[:,2,None]*mask).sum(0)
            batch=np.empty((M+1,T,C));batch[0]=y@mask
            for tr,tid in enumerate(trial_ids):
                # Region and Unit_ID intentionally absent: simultaneously recorded signals share shifts.
                random=rng(c["statistics"]["seed"],"glm_shifts",meta["Group"],meta["Mouse_ID"],meta["Session_ID"],tid)
                shifts=random.integers(0,B,size=M)
                batch[1:,tr,:]=y[tr][(np.arange(B)[None,:]-shifts[:,None])%B]@mask
            gain=np.zeros(M+1);term_gain={t:0. for t in terms};total_spikes=0.
            for fold,train in splits:
                test=~train;mean=orders[train].mean();sd=orders[train].std(ddof=0)
                if sd<=0:raise ValueError("Constant training trial order")
                z=(orders-mean)/sd
                X=np.column_stack([np.ones(T*C),np.repeat(z,C),np.tile(categories,(T,1))]);e=np.tile(exposure,T)
                Y=batch.reshape(M+1,-1);tr=np.repeat(train,C);te=np.repeat(test,C)
                bet,niter=fit_batch(X[tr],Y[:,tr],e[tr],opt)
                # Reduced model depends only on trial totals, invariant under within-trial shifts.
                reduced,nred=fit_batch(X[tr,:2],Y[0,tr],e[tr],opt)
                eta=bet@X[te].T;eta0=reduced@X[te,:2].T
                ll=(Y[:,te]*eta-np.exp(eta)*e[te]).sum(1)
                ll0=(Y[:,te]*eta0-np.exp(eta0)*e[te]).sum(1)
                delta=ll-ll0;gain+=delta;spikes=Y[0,te].sum();total_spikes+=spikes
                local["folds"].append(meta|dict(Fold=fold,N_Train_Trials=int(train.sum()),N_Test_Trials=int(test.sum()),Heldout_Spikes=spikes,LL_Full_Constant_Omitted=ll[0],LL_Reduced_Constant_Omitted=ll0[0],Delta_LL_Nats=delta[0],Gain_Bits_Per_Spike=delta[0]/(spikes*np.log(2)) if spikes else np.nan))
                for model,names,values in [("Full",["Intercept","Trial_Order_Z"]+terms,bet[0]),("Reduced",["Intercept","Trial_Order_Z"],reduced[0])]:
                    for name,value in zip(names,values):local["coefs"].append(meta|dict(Fold=fold,Model=model,Term=name,Beta=value,Exp_Beta=np.exp(value),Training_Order_Mean=mean,Training_Order_SD=sd))
                for j,term in enumerate(terms):
                    keep=[k for k in range(X.shape[1]) if k!=j+2]
                    br,_=fit_batch(X[tr][:,keep],Y[0,tr],e[tr],opt)
                    er=br@X[te][:,keep].T
                    lr=(Y[0,te]*er[0]-np.exp(er[0])*e[te]).sum()
                    term_gain[term]+=ll[0]-lr
                # Predictions and Pearson dispersion are evaluated on original bins.
                dispersion_sum=0.;training_bins=0
                for k,tid in enumerate(trial_ids):
                    xb=np.column_stack([np.ones(B),np.full(B,z[k]),reference[:,3:]])
                    pred=np.exp(xb@bet[0]);pred0=np.exp(xb[:,:2]@reduced[0]);mu=pred*reference[:,2]
                    local["splits"].append(meta|dict(Fold=fold,Trial_ID=tid,Trial_Order=int(orders[k]),Role="Train" if train[k] else "Validation"))
                    if train[k]:
                        dispersion_sum+=np.sum((y[k]-mu)**2/mu);training_bins+=B
                    else:
                        for j in range(B):local["preds"].append(meta|dict(Trial_ID=tid,Time_Bin_Start_s=reference[j,0],Observed_FR_Hz=y[k,j]/reference[j,2],Predicted_Full_FR_Hz=pred[j],Predicted_Reduced_FR_Hz=pred0[j]))
                local["diagnostics"].append(meta|dict(Fold=fold,Full_Iterations=niter,Reduced_Iterations=nred,Train_Pearson_Dispersion=dispersion_sum/(training_bins-X.shape[1]) if training_bins>X.shape[1] else np.nan,N_Train_Bins=training_bins,Status="converged"))
            scores=gain/(total_spikes*np.log(2))
            for term,value in [("All_Events",gain[0]),*term_gain.items()]:
                local["units"].append(meta|dict(Term=term,Gain_Bits_Per_Spike=value/(total_spikes*np.log(2)),Delta_LL_Nats=value,Heldout_Spikes=total_spikes,N_Trials=T))
            for k,val in enumerate(scores[1:],1):local["perms"].append(meta|dict(Permutation_ID=k,Gain_Bits_Per_Spike=val))
            for target,key in [(unit_rows,"units"),(fold_rows,"folds"),(coef_rows,"coefs"),(pred_rows,"preds"),(diagnostics,"diagnostics"),(perm_rows,"perms"),(split_rows,"splits")]:target.extend(local[key])
        except (ValueError,np.linalg.LinAlgError,FloatingPointError) as e:
            failures.append(meta|dict(Reason=str(e)))
        if number%25==0:print(f"GLM {number}/{len(units)} unit-sessions",flush=True)
    save(pd.DataFrame(failures,columns=UNIT+["Reason"]),out,"audit/glm_failures.csv")
    if failures and opt["failure_policy"]=="error":raise ValueError("GLM failures: inspect audit/glm_failures.csv. No GLM results finalized.")
    if not unit_rows:raise ValueError("No estimable GLM units")
    for rows,name in [(fold_rows,"folds"),(coef_rows,"coefficients"),(diagnostics,"diagnostics"),(split_rows,"trial_splits")]:save(pd.DataFrame(rows),out,f"tables/glm_{name}.csv")
    u=save(pd.DataFrame(unit_rows),out,"plotdata/glm_unit_gains.csv")
    a=u.groupby(ANIMAL+["Term"],as_index=False).agg(Gain_Bits_Per_Spike=("Gain_Bits_Per_Spike","mean"),N_Units=("Unit_ID","size"))
    save(a,out,"plotdata/glm_animal_gains.csv")
    pr=pd.DataFrame(pred_rows).groupby(UNIT+[BIN[0]],as_index=False)[["Observed_FR_Hz","Predicted_Full_FR_Hz","Predicted_Reduced_FR_Hz"]].mean()
    save(pr,out,"plotdata/glm_unit_predictions.csv")
    pa=pr.groupby(ANIMAL+[BIN[0]],as_index=False)[["Observed_FR_Hz","Predicted_Full_FR_Hz","Predicted_Reduced_FR_Hz"]].mean();save(pa,out,"plotdata/glm_animal_predictions.csv")
    perm=pd.DataFrame(perm_rows).groupby(ANIMAL+["Permutation_ID"],as_index=False).Gain_Bits_Per_Spike.mean()
    save(perm,out,"plotdata/glm_animal_null.csv")
    null=perm.groupby(["Group","Region","Stage","Permutation_ID"],as_index=False).Gain_Bits_Per_Spike.mean();save(null,out,"plotdata/glm_group_null.csv")
    rows=[]
    for group,gs in c["groups"].items():
        for reg in c["regions"]:
            for stage in gs["stages"]:
                if stage not in opt["stages"]:continue
                q=a[a.Group.eq(group)&a.Region.eq(reg)&a.Stage.eq(stage)&a.Term.eq("All_Events")]
                row=summary(q.Gain_Bits_Per_Spike,c["statistics"],"glm",group,reg,stage)
                nv=null[null.Group.eq(group)&null.Region.eq(reg)&null.Stage.eq(stage)].Gain_Bits_Per_Spike.to_numpy()
                row.update(Group=group,Region=reg,Stage=stage,Term="All_Events",P_Raw=(1+np.sum(nv>=row["Mean"]-1e-12))/(len(nv)+1) if len(nv) and np.isfinite(row["Mean"]) else np.nan,N_Shifts=len(nv),N_Units=int(q.N_Units.sum()),Family="glm_event_locking",Test="upper_tail_shared_circular_shift")
                rows.append(row)
    save(adjust(pd.DataFrame(rows),opt["correction"]),out,"statistics/glm_event_locking.csv")
    desc=[]
    for ident,q in a[a.Term.ne("All_Events")].groupby(["Group","Region","Stage","Term"]):desc.append(dict(zip(["Group","Region","Stage","Term"],ident))|summary(q.Gain_Bits_Per_Spike,c["statistics"],"glm_drop_one",*ident))
    save(pd.DataFrame(desc),out,"statistics/glm_drop_one_descriptive.csv")
