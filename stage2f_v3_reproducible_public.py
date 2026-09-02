from pathlib import Path
import csv,gzip,math
import numpy as np

ROOT=Path.cwd()
DATA=ROOT/"MotorUnit_files_for_analysis"/"Stage2E_model_ready"/"Stage2E_primary_model_ready_public.csv.gz"
OUT=ROOT/"MotorUnit_files_for_analysis"/"Stage2F_v3_final_crossed_hierarchical"; OUT.mkdir(parents=True,exist_ok=True)
N_CHAINS,N_ITER,BURN_IN,THIN=4,8000,3000,5
BASE_SEED,N_SWEEPS=260817,3
BETA_PRIOR_SD,IG_SHAPE,IG_RATE=2.5,2.5,0.5
PRIMARY="plateau_train_rate_hz"; ALT="plateau_mean_instantaneous_rate_hz"

def read_rows(p):
    op=gzip.open if str(p).endswith(".gz") else open
    with op(p,"rt",encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def fnum(x):
    try:
        v=float(x); return v if math.isfinite(v) else None
    except:return None
def fint(x):
    try:return int(float(x))
    except:return None
def prep(rows,outcome):
    pre=[]
    for r in rows:
        p=r.get("participant_folder","").strip(); rec=r.get("recording_id","").strip()
        b=r.get("stage2e_bout_uid","").strip(); mu=fint(r.get("original_mu_number"))
        y=fnum(r.get(outcome)); a=fnum(r.get("muap_peak_to_peak_mean_raw"))
        if not p or not rec or not b or mu is None or y is None or a is None or a<=0:continue
        q=dict(r); q["_y"]=y; q["_la"]=math.log(a); q["_mu"]=f"{p}||{rec}||{mu}"; pre.append(q)
    by={}
    for r in pre:by.setdefault(r["stage2e_bout_uid"],[]).append(r)
    out=[]
    for rs in by.values():
        if len(rs)<2:continue
        z=np.array([r["_la"] for r in rs]); sd=float(z.std(ddof=1))
        if not math.isfinite(sd) or sd<=0:continue
        m=float(z.mean())
        for r in rs:r["_z"]=(r["_la"]-m)/sd;out.append(r)
    return out
def index(v):
    lev=sorted(set(v)); d={x:i for i,x in enumerate(lev)}
    return np.array([d[x] for x in v],int),lev
def design(rows,limb_mod=False):
    raw=[]
    for r in rows:
        raw.append((r["_z"],r["assigned_intensity_percent"]=="50",r["assigned_intensity_percent"]=="75",
                    r["muscle"]=="VL",r["limb"]=="ACL",r["timepoint"]=="2wks",r["timepoint"]=="6wks",
                    r["timepoint"]=="3months",r["timepoint"]=="6months"))
    A=np.array(raw,float); means=A[:,1:9].mean(0)
    X=np.column_stack([A[:,0],A[:,1:9]-means])
    names=["z_log_MUAP_within_bout","Intensity_50_vs_20","Intensity_75_vs_20","Muscle_VL_vs_VM",
           "Limb_ACL_vs_Opp","Time_2wks_vs_Preop","Time_6wks_vs_Preop","Time_3months_vs_Preop","Time_6months_vs_Preop"]
    if limb_mod:X=np.column_stack([X,A[:,0]*A[:,4]]);names.append("zMUAP_x_LimbACL")
    return X,names
def ig(rng,shape,rate):return 1/rng.gamma(shape=shape,scale=1/rate)
def init(y,X,pi,mi,bi,rng):
    p=X.shape[1]; beta=np.linalg.solve(X.T@X+np.eye(p)*1e-6,X.T@y)+rng.normal(0,.01,p)
    res=y-X@beta
    def re_init(idx,n,res,w):
        c=np.bincount(idx,minlength=n).astype(float); s=np.bincount(idx,weights=res,minlength=n)
        return (s/np.maximum(c,1))*w
    npart,nmu,nb=pi.max()+1,mi.max()+1,bi.max()+1
    pr=re_init(pi,npart,res,.25); res=res-pr[pi]
    mr=re_init(mi,nmu,res,.35); res=res-mr[mi]
    br=re_init(bi,nb,res,.35); res=res-br[bi]
    vv=lambda x:max(float(np.var(x,ddof=1)),.05)
    return beta,pr,mr,br,vv(res),vv(pr),vv(mr),vv(br)
def chain(yraw,X,pi,mi,bi,seed):
    rng=np.random.default_rng(seed); n=len(yraw); p=X.shape[1]
    ym=float(yraw.mean()); ys=float(yraw.std(ddof=1)); y=(yraw-ym)/ys
    npart,nmu,nb=pi.max()+1,mi.max()+1,bi.max()+1
    pc=np.bincount(pi,minlength=npart).astype(float);mc=np.bincount(mi,minlength=nmu).astype(float);bc=np.bincount(bi,minlength=nb).astype(float)
    beta,pr,mr,br,se2,sp2,sm2,sb2=init(y,X,pi,mi,bi,rng)
    prior=np.eye(p)/(BETA_PRIOR_SD**2); XtX=X.T@X; keep=[]
    for it in range(N_ITER):
        ay=y-pr[pi]-mr[mi]-br[bi]; P=XtX/se2+prior; rhs=X.T@ay/se2
        L=np.linalg.cholesky(P); pm=np.linalg.solve(L.T,np.linalg.solve(L,rhs)); beta=pm+np.linalg.solve(L.T,rng.normal(size=p))
        for _ in range(N_SWEEPS):
            rr=y-X@beta-mr[mi]-br[bi]; v=1/(pc/se2+1/sp2); m=v*np.bincount(pi,weights=rr,minlength=npart)/se2;pr=rng.normal(m,np.sqrt(v))
            rr=y-X@beta-pr[pi]-br[bi]; v=1/(mc/se2+1/sm2); m=v*np.bincount(mi,weights=rr,minlength=nmu)/se2;mr=rng.normal(m,np.sqrt(v))
            rr=y-X@beta-pr[pi]-mr[mi]; v=1/(bc/se2+1/sb2); m=v*np.bincount(bi,weights=rr,minlength=nb)/se2;br=rng.normal(m,np.sqrt(v))
        rr=y-X@beta-pr[pi]-mr[mi]-br[bi]
        se2=ig(rng,IG_SHAPE+n/2,IG_RATE+float(rr@rr)/2)
        sp2=ig(rng,IG_SHAPE+npart/2,IG_RATE+float(pr@pr)/2)
        sm2=ig(rng,IG_SHAPE+nmu/2,IG_RATE+float(mr@mr)/2)
        sb2=ig(rng,IG_SHAPE+nb/2,IG_RATE+float(br@br)/2)
        if it>=BURN_IN and (it-BURN_IN)%THIN==0:keep.append(beta*ys)
    return np.asarray(keep)
def seed(name,c):return BASE_SEED+sum((i+1)*ord(ch) for i,ch in enumerate(name))+10000*c
def rhat(C):
    m,n=C.shape;h=n//2;S=np.vstack([C[:,:h],C[:,n-h:]]);n2=S.shape[1]
    W=float(np.mean(np.var(S,axis=1,ddof=1)))
    if W<=0:return 1.0
    B=float(n2*np.var(S.mean(1),ddof=1)); return math.sqrt(max((((n2-1)/n2)*W+B/n2)/W,1))
def summ(x):
    return float(x.mean()),float(x.std(ddof=1)),float(np.quantile(x,.025)),float(np.quantile(x,.975)),float((x>0).mean()),float((x<0).mean())
def run(name,rows,outcome,limb=False):
    X,names=design(rows,limb); y=np.array([r["_y"] for r in rows])
    pi,pl=index([r["participant_folder"] for r in rows]);mi,ml=index([r["_mu"] for r in rows]);bi,bl=index([r["stage2e_bout_uid"] for r in rows])
    C=[chain(y,X,pi,mi,bi,seed(name,c)) for c in range(1,N_CHAINS+1)]
    rowsout=[];draw={}
    for j,p in enumerate(names):
        M=np.vstack([c[:,j] for c in C]); x=M.ravel(); a=summ(x)
        rowsout.append(dict(model=name,outcome=outcome,parameter=p,posterior_mean=a[0],posterior_sd=a[1],lower_95_CrI=a[2],upper_95_CrI=a[3],Pr_gt_0=a[4],Pr_lt_0=a[5],split_R_hat=rhat(M),n_rows=len(rows),n_participants=len(pl),n_recording_MUs=len(ml),n_bouts=len(bl)))
        draw[p]=M
    return rowsout,draw
def write(path,rows,fields):
    with open(path,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

src=read_rows(DATA)
sets={"primary":src,
      "sensitivity_5plus":[r for r in src if r["stage2e_sensitivity_5plus_eligible"]=="Yes"],
      "sensitivity_A_or_B":[r for r in src if r["stage2e_sensitivity_A_or_B_eligible"]=="Yes"],
      "sensitivity_A_clean":[r for r in src if r["stage2e_sensitivity_Aclean_eligible"]=="Yes"]}
P={k:prep(v,PRIMARY) for k,v in sets.items()};P["sensitivity_mean_instantaneous_rate"]=prep(src,ALT)
models=[
 ("primary_overall",P["primary"],PRIMARY,False),
 ("secondary_limb_modifier",P["primary"],PRIMARY,True),
 ("sensitivity_5plus_firings",P["sensitivity_5plus"],PRIMARY,False),
 ("sensitivity_A_or_B",P["sensitivity_A_or_B"],PRIMARY,False),
 ("sensitivity_A_clean",P["sensitivity_A_clean"],PRIMARY,False),
 ("sensitivity_mean_instantaneous_rate",P["sensitivity_mean_instantaneous_rate"],ALT,False)]
allrows=[];draws={}
for name,rows,out,limb in models:
    print(f"Fitting {name}: {len(rows)} rows")
    rr,dd=run(name,rows,out,limb);allrows+=rr;draws[name]=dd
fields=["model","outcome","parameter","posterior_mean","posterior_sd","lower_95_CrI","upper_95_CrI","Pr_gt_0","Pr_lt_0","split_R_hat","n_rows","n_participants","n_recording_MUs","n_bouts"]
write(OUT/"Stage2F_v3_model_posterior_summary.csv",allrows,fields)

M=draws["secondary_limb_modifier"]["z_log_MUAP_within_bout"]
I=draws["secondary_limb_modifier"]["zMUAP_x_LimbACL"]
slopes=[]
for lab,D in [("Opp",M),("ACL",M+I),("ACL_minus_Opp",I)]:
    x=D.ravel();a=summ(x);slopes.append(dict(model="secondary_limb_modifier",slope=lab,posterior_mean_Hz_per_1SD_logMUAP=a[0],posterior_sd=a[1],lower_95_CrI=a[2],upper_95_CrI=a[3],Pr_gt_0=a[4],Pr_lt_0=a[5]))
write(OUT/"Stage2F_v3_limb_modifier_slopes.csv",slopes,list(slopes[0]))

key=[]
for name,param in [("primary_overall","z_log_MUAP_within_bout"),("secondary_limb_modifier","z_log_MUAP_within_bout"),("secondary_limb_modifier","zMUAP_x_LimbACL"),
 ("sensitivity_5plus_firings","z_log_MUAP_within_bout"),("sensitivity_A_or_B","z_log_MUAP_within_bout"),("sensitivity_A_clean","z_log_MUAP_within_bout"),("sensitivity_mean_instantaneous_rate","z_log_MUAP_within_bout")]:
    M=draws[name][param]
    for c in range(M.shape[0]):
        for d,v in enumerate(M[c],1):key.append(dict(model=name,parameter=param,chain=c+1,draw=d,value_Hz=float(v)))
write(OUT/"Stage2F_v3_key_posterior_draws.csv",key,["model","parameter","chain","draw","value_Hz"])
print("Done.")
