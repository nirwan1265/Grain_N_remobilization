# LD of the chr6 Q.E lead SNP 6-7041072 (AGPv4 7,395,538; v5 12,596,405, in Zm00001eb261470 aaap40)
# with all MAF>=0.05 HapMap 3.2.1 SNPs in chr6 v4 7.0-7.8 Mb; per-gene max r2 for the AAAP/LHT cluster.
# Output: final/results/goodman_buckler_amino_acid_gwas/LD/chr6_QE_lead_6-7041072_r2.csv
import gzip, numpy as np, sys
VCF="/sessions/rcw-01fbiuktcbzxwfqvxivo36ex/mnt/Sarah_amino_acid/genotype/hmp321_282_agpv4_chr6.vcf.gz"
LO,HI=7_000_000,7_800_000
pos=[];G=[];ra=[]
with gzip.open(VCF,"rt") as f:
    for l in f:
        if l[0]=="#":
            if l.startswith("#CHROM"): n=len(l.split("\t"))-9
            continue
        t=l.split("\t",10); p=int(t[1])
        if p<LO: continue
        if p>HI: break
        if "," in t[4] or len(t[3])!=1 or len(t[4])!=1: continue
        gt=l.rstrip("\n").split("\t")[9:]
        d=np.array([ (np.nan if g[0]=="." else int(g[0])+int(g[2])) for g in gt],float)
        pos.append(p);G.append(d);ra.append(t[3]+t[4])
np.savez_compressed("./chr6_7.0_7.8.npz",pos=np.array(pos),G=np.array(G),ra=np.array(ra))
print(len(pos),n)
import numpy as np, pandas as pd
W="./"
z=np.load(W+"chr6_7.0_7.8.npz"); pos=z["pos"]; G=z["G"]
# chain lift chr6 v4 -> v5 (same-chr + strand, highest score)
CH="/sessions/rcw-01fbiuktcbzxwfqvxivo36ex/mnt/Sarah_amino_acid/chain/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain.txt"
blocks=[]
with open(CH) as f:
    cur=None
    for l in f:
        t=l.split()
        if not t: cur=None; continue
        if t[0]=="chain":
            sc=float(t[1]); tn=t[2]; ts=int(t[5]); qn=t[7]; qst=t[9]; qs=int(t[10])
            ok = tn.replace("chr","")=="6" and qn.replace("chr","")=="6" and qst=="+" and t[4]=="+"
            cur=[sc,ts,qs] if ok else None; continue
        if cur is None: continue
        size=int(t[0]); blocks.append((cur[1],cur[1]+size,cur[2],cur[0]))
        if len(t)==3: cur[1]+=size+int(t[1]); cur[2]+=size+int(t[2])
B=np.array(blocks); B=B[(B[:,1]>6.9e6)&(B[:,0]<7.9e6)]
def lift(p):
    m=B[(B[:,0]<=p)&(B[:,1]>p)]
    if len(m)==0: return np.nan
    b=m[np.argmax(m[:,3])]; return int(b[2]+(p-b[0]))
v5=np.array([lift(p) for p in pos])
# genes v5
g=[]
for l in open("/sessions/rcw-01fbiuktcbzxwfqvxivo36ex/mnt/Grain_N_remobilization/final/results/annotation/Zm-B73-REFERENCE-NAM-5.0_genes_only.gff3"):
    if l.startswith("#"): continue
    f=l.split("\t")
    if len(f)<9 or f[2]!="gene" or f[0] not in("chr6","6"): continue
    s,e=int(f[3]),int(f[4])
    if e<12.2e6 or s>13.1e6: continue
    g.append((f[8].split("ID=")[1].split(";")[0],s,e,f[6]))
g=pd.DataFrame(g,columns=["gid","s","e","strand"])
maf=np.nanmean(G,1)/2; maf=np.minimum(maf,1-maf)
keep=maf>=0.05
lead=np.where(pos==7395538)[0][0]; print("lead maf",maf[lead],"v5",v5[lead])
def r2(i,J):
    x=G[i]; out=[]
    for j in J:
        y=G[j]; ok=~np.isnan(x)&~np.isnan(y)
        if ok.sum()<50 or np.std(x[ok])==0 or np.std(y[ok])==0: out.append(np.nan); continue
        out.append(np.corrcoef(x[ok],y[ok])[0,1]**2)
    return np.array(out)
J=np.where(keep)[0]; R=r2(lead,J)
d=pd.DataFrame({"pos_v4":pos[J],"pos_v5":v5[J],"maf":maf[J],"r2":R})
d.to_csv(W+"lead_r2.csv",index=False)
print("SNPs maf>=.05",len(d),"r2>=0.8",(d.r2>=0.8).sum()-1,"r2>=0.4",(d.r2>=0.4).sum()-1)
print(d[d.r2>=0.4].sort_values("pos_v5").to_string(index=False))
# per gene max r2
for _,r in g.iterrows():
    m=d[(d.pos_v5>=r.s-2000)&(d.pos_v5<=r.e+2000)]
    print(r.gid,r.s,r.e,r.strand,len(m),"max r2 %.2f"%(m.r2.max() if len(m) else np.nan))
his=np.where(pos==7518003)[0]; print("His SNP r2 with lead", r2(lead,his) if len(his) else "absent")
