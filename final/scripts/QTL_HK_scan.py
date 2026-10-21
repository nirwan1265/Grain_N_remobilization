import pandas as pd, numpy as np, json
# R/qtl csvr cross file (B64 x CML14, Asp and Pro BLUPs)
raw=pd.read_csv("data/QTL/Nitrogen_QTL_map.csv",header=None,dtype=str)
hdr=raw.iloc[0].tolist()
mk=pd.DataFrame({"m":hdr[3:],"chr":raw.iloc[1,3:].astype(int).values,"cM":pd.to_numeric(raw.iloc[2,3:]).values,"col":range(3,len(hdr))})
mk["pos"]=mk.m.str.extract(r"_(\d+)$")[0].astype(int)
code={"AA":0,"AB":1,"BB":2}
G=raw.iloc[3:,3:].replace(code).apply(pd.to_numeric,errors="coerce").values  # n x M
ids=raw.iloc[3:,0].values; ph=raw.iloc[3:,1:3].apply(pd.to_numeric).values
n=G.shape[0]; eps=1e-4
def T(r): return np.array([[(1-r)**2,2*r*(1-r),r**2],[r*(1-r),(1-r)**2+r**2,r*(1-r)],[r**2,2*r*(1-r),(1-r)**2]])
def emis(o):
    e=np.full((len(o),3),1.0)
    ok=~np.isnan(o)
    for k in range(3): e[ok,k]=np.where(o[ok]==k,1-eps,eps/2)
    return e
grid=[]; probs=[]
for c in range(1,11):
    m=mk[mk.chr==c].sort_values("cM").reset_index(drop=True)
    cm=m.cM.values
    # grid: markers + pseudomarkers every 1 cM
    pts=np.unique(np.round(np.concatenate([cm,np.arange(cm.min(),cm.max(),1.0)]),6))
    ismk={round(x,6):i for i,x in enumerate(cm)}
    obs=np.full((n,len(pts)),np.nan)
    for j,x in enumerate(pts):
        if round(x,6) in ismk: obs[:,j]=G[:,m.col.values[ismk[round(x,6)]]-3]
    d=np.diff(pts); r=0.5*(1-np.exp(-2*d/100)); Ts=[T(x) for x in r]
    L=len(pts); P=np.zeros((n,L,3))
    for i in range(n):
        E=emis(obs[i]); a=np.zeros((L,3)); b=np.ones((L,3))
        a[0]=np.array([.25,.5,.25])*E[0]; a[0]/=a[0].sum()
        for j in range(1,L): a[j]=(a[j-1]@Ts[j-1])*E[j]; a[j]/=a[j].sum()
        for j in range(L-2,-1,-1): b[j]=Ts[j]@(E[j+1]*b[j+1]); b[j]/=b[j].sum()
        pp=a*b; P[i]=pp/pp.sum(1,keepdims=True)
    for j,x in enumerate(pts):
        near=m.iloc[(np.abs(cm-x)).argmin()]
        grid.append((c,x,round(x,6) in ismk,near.m,near.pos))
    probs.append(P)
P=np.concatenate(probs,axis=1)  # n x positions x 3
grid=pd.DataFrame(grid,columns=["chr","cM","is_marker","nearest_marker","nearest_pos"])
np.save("final/results/qtl/P.npy",P); grid.to_csv("final/results/qtl/grid.csv",index=False)
rng=np.random.default_rng(1); NP=1000; res={}
for t,name in enumerate(["Aspartate","Proline"]):
    y=ph[:,t]; ok=~np.isnan(y); yy=y[ok]; nn=len(yy)
    Y=np.column_stack([yy]+[rng.permutation(yy) for _ in range(NP)]); Y=Y-Y.mean(0)
    rss0=(Y**2).sum(0)
    lod=np.zeros((P.shape[1],NP+1))
    for j in range(P.shape[1]):
        X=np.column_stack([np.ones(nn),P[ok,j,0],P[ok,j,1]])
        Q,_=np.linalg.qr(X); proj=Q.T@Y; rss1=rss0-(proj[1:]**2).sum(0)  # Y centred: first col of Q is constant
        lod[j]=nn/2*np.log10(rss0/rss1)
    thr=np.quantile(lod[:,1:].max(0),0.95)
    grid[name]=lod[:,0]; res[name]={"n":int(nn),"thr95":float(thr)}
grid.to_csv("final/results/qtl/lod.csv",index=False); json.dump(res,open("final/results/qtl/res.json","w")); print(res)
