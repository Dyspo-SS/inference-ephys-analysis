"""Animal-level resampling; no units, bins, or resamples count as animals."""
import itertools
import math
import numpy as np
from .common import rng

def finite(values):
    a = np.asarray(values, float)
    return a[np.isfinite(a)]

def summary(values, cfg, *key):
    a = finite(values)
    result = dict(N_Animals=len(a), Mean=float(a.mean()) if len(a) else np.nan,
                  CI_Lower=np.nan, CI_Upper=np.nan)
    if len(a) >= 2:
        r = rng(cfg["seed"], "bootstrap", *key)
        boots = a[r.integers(0, len(a), (cfg["bootstrap_draws"], len(a)))].mean(1)
        alpha = (1 - cfg["ci_level"]) / 2
        result.update(zip(["CI_Lower", "CI_Upper"], np.quantile(boots, [alpha, 1-alpha])))
    return result

def signflip(values, cfg, *key):
    a = finite(values); n = len(a)
    out = summary(a, cfg, *key)
    out.update(P_Raw=np.nan, Effect_dz=np.nan, Assignments=0, Test="unavailable")
    if n < 2:
        return out
    out["Effect_dz"] = a.mean()/a.std(ddof=1) if a.std(ddof=1) > 0 else np.nan
    threshold = abs(a.mean()) - 1e-12
    total = 2**n
    exceed = 0
    if total <= cfg["exact_limit"]:
        # Chunked enumeration prevents exponential memory allocation.
        for start in range(0, total, 8192):
            ids = np.arange(start, min(start+8192, total), dtype=np.uint64)
            signs = 2*((ids[:, None] >> np.arange(n, dtype=np.uint64)) & 1).astype(float)-1
            exceed += np.count_nonzero(np.abs(signs@a/n) >= threshold)
        out.update(P_Raw=exceed/total, Assignments=total, Test="exact_two_sided_sign_flip")
    else:
        r = rng(cfg["seed"], "signflip", *key); total = cfg["permutation_draws"]
        for start in range(0, total, 8192):
            signs = r.choice([-1., 1.], size=(min(8192, total-start), n))
            exceed += np.count_nonzero(np.abs(signs@a/n) >= threshold)
        out.update(P_Raw=(exceed+1)/(total+1), Assignments=total, Test="monte_carlo_two_sided_sign_flip")
    return out

def compare(a, b, cfg, *key):
    a,b = finite(a),finite(b); na,nb = len(a),len(b)
    out = dict(N_A=na, N_B=nb, Mean_A=a.mean() if na else np.nan,
               Mean_B=b.mean() if nb else np.nan, Difference_A_Minus_B=np.nan,
               CI_Lower=np.nan, CI_Upper=np.nan, P_Raw=np.nan, Hedges_g=np.nan,
               Assignments=0, Test="unavailable")
    if min(na,nb) < 2:
        return out
    obs = a.mean()-b.mean(); pool=np.r_[a,b]; total=math.comb(na+nb,na)
    exceed=0; threshold=abs(obs)-1e-12
    if total <= cfg["exact_limit"]:
        iterator=itertools.combinations(range(na+nb),na)
        while True:
            batch=list(itertools.islice(iterator,8192))
            if not batch: break
            s=pool[np.asarray(batch)].sum(1)
            exceed+=np.count_nonzero(np.abs(s/na-(pool.sum()-s)/nb)>=threshold)
        p=exceed/total; test="exact_two_sided_group_permutation"
    else:
        total=cfg["permutation_draws"]; r=rng(cfg["seed"],"labels",*key)
        for _ in range(total):
            idx=r.permutation(len(pool)); v=pool[idx[:na]].mean()-pool[idx[na:]].mean()
            exceed+=abs(v)>=threshold
        p=(exceed+1)/(total+1);test="monte_carlo_two_sided_group_permutation"
    r=rng(cfg["seed"],"difference_bootstrap",*key); m=cfg["bootstrap_draws"]
    boots=a[r.integers(0,na,(m,na))].mean(1)-b[r.integers(0,nb,(m,nb))].mean(1)
    alpha=(1-cfg["ci_level"])/2;lo,hi=np.quantile(boots,[alpha,1-alpha])
    pooled=np.sqrt(((na-1)*a.var(ddof=1)+(nb-1)*b.var(ddof=1))/(na+nb-2))
    out.update(Difference_A_Minus_B=obs, CI_Lower=lo, CI_Upper=hi, P_Raw=p,
               Hedges_g=(1-3/(4*(na+nb)-9))*obs/pooled if pooled>0 else np.nan,
               Assignments=total, Test=test)
    return out

def adjust(table, method="Holm"):
    """Keep every planned hypothesis; unavailable p values consume a slot as p=1."""
    t=table.copy()
    t["P_Adjusted"]=np.nan;t["Correction"]=method;t["Family_Size"]=len(t)
    if not len(t): return t
    p=t.P_Raw.to_numpy(float);good=np.isfinite(p);work=np.where(good,p,1.)
    order=np.argsort(work);ranked=work[order];n=len(p)
    if method=="Holm": adj=np.minimum(1,np.maximum.accumulate(ranked*np.arange(n,0,-1)))
    elif method=="BH": adj=np.minimum(1,np.minimum.accumulate((ranked*n/np.arange(1,n+1))[::-1])[::-1])
    elif method=="None": adj=ranked
    else: raise ValueError("Unknown correction method")
    restored=np.empty(n);restored[order]=adj;restored[~good]=np.nan
    t["P_Adjusted"]=restored
    return t
