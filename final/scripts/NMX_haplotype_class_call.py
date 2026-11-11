#!/usr/bin/env python3
"""Steps 3-5 of the NMX-hsftf9 haplotype check (sanity analysis, not for the PlantPhys paper).

  kinship  build a LOCO kinship (chromosomes 1-8 and 10, 20,000 random SNPs, MAF >= 0.05)
           for the 3,107 soil-N landraces and cache its eigendecomposition
  classify (step 3-4) diagnostic SNPs from the assemblies (gbs_alleles/*.tsv + the interval
           type per genome from snp_alleles_by_type.tsv), then a long-version score per landrace
  models   (step 5) mixed models with kinship, each site variable ~ predictor, for the predictors
           class score, lead SNP (151,742,636) and 152 Mb SNP (151,983,481); plus joint models
           soil N ~ predictor + other site variables and elevation ~ predictor + other site variables

usage: python3 NMX_haplotype_class_call.py <kinship|classify|models> <repo> <bed_prefix> <class_call_dir> <cache_dir>
"""
import sys, os, glob, csv, numpy as np, pandas as pd
mode, repo, bed, CC, CACHE = sys.argv[1:6]
os.makedirs(CACHE, exist_ok=True); OUT = os.path.join(CC, "results"); os.makedirs(OUT, exist_ok=True)
SN = os.path.join(repo, "final/results/soil_nitrogen_gwas")
fam = [l.split()[1] for l in open(bed + ".fam")]; n = len(fam)
bim = pd.read_csv(bed + ".bim", sep=r"\s+", header=None, names=["chr", "id", "cm", "pos", "a1", "a2"])
raw = np.memmap(bed + ".bed", dtype=np.uint8, mode="r", offset=3).reshape(len(bim), (n + 3) // 4)
lut = np.array([2, np.nan, 1, 0])                                    # dosage of a1
env = pd.read_csv(os.path.join(SN, "control_phenotypes_maize_landrace.csv")).drop_duplicates("Taxa").set_index("Taxa")
keepset = set(open(os.path.join(SN, "LD/landrace_N_keep.txt")).read().split())
keep = np.array([i for i, t in enumerate(fam) if t in keepset and t in env.index]); taxa = [fam[i] for i in keep]
E = env.loc[taxa]
def G(idx):
    idx = np.atleast_1d(idx); b = np.asarray(raw[idx])
    g = np.stack([(b >> s) & 3 for s in (0, 2, 4, 6)], axis=-1).reshape(len(idx), -1)[:, :n]
    return lut[g][:, keep]
def snp(pos): return G(bim.index[(bim.chr == 9) & (bim.pos == pos)][0])[0]

if mode == "kinship":
    rng = np.random.default_rng(1); cand = bim.index[bim.chr != 9].values
    M = G(np.sort(rng.choice(cand, 40000, replace=False))); f = np.nanmean(M, 1) / 2
    M = M[(f >= 0.05) & (f <= 0.95)][:20000]; f = np.nanmean(M, 1) / 2
    Z = M - 2 * f[:, None]; Z[np.isnan(Z)] = 0
    K = Z.T @ Z / (2 * (f * (1 - f)).sum())
    lam, U = np.linalg.eigh(K); np.savez(os.path.join(CACHE, "kin.npz"), lam=lam, U=U, taxa=np.array(taxa))
    print("kinship", K.shape, "from", M.shape[0], "SNPs; eigenvalues", lam.min().round(3), lam.max().round(1))

elif mode == "classify":
    ty = {r["genome"]: r["type"] for r in csv.DictReader(open(os.path.join(CC, "..", "snp_alleles_by_type.tsv")), delimiter="\t")}
    reg = pd.read_csv(os.path.join(CC, "gbs_snps_region.tsv"), sep="\t").set_index("pos")
    A = {}
    for f in glob.glob(os.path.join(CC, "gbs_alleles", "*.tsv")):
        d = pd.read_csv(f, sep="\t"); g = d.genome.iloc[0]
        if not g.startswith("Zm-") or ty.get(g, "NA") in ("NA", ""): continue
        line = g.split("-")[1]
        if line in A or line == "B73_AB10": continue                   # one assembly per line
        d = d[(d.n_good == 1) & d.allele.isin(list("ACGT"))]
        A[line] = (("long" if ty[g] == "long" else "short"), dict(zip(d.pos, d.allele)))
    nl = sum(1 for v in A.values() if v[0] == "long"); ns = len(A) - nl
    rows = []
    for pos, r in reg.iterrows():
        al = {"long": [], "short": []}
        for cls, gt in A.values():
            if pos in gt and gt[pos] in (r.a1, r.a2): al[cls].append(gt[pos])
        if len(al["long"]) < 0.8 * nl or len(al["short"]) < 0.8 * ns: continue
        for a, b in ((r.a1, r.a2), (r.a2, r.a1)):
            fl = al["long"].count(a) / len(al["long"]); fs = al["short"].count(b) / len(al["short"])
            if fl >= 0.95 and fs >= 0.95:
                rows.append(dict(pos=pos, snp=r.snp, long_allele=a, short_allele=b, n_long=len(al["long"]), n_short=len(al["short"]), frac_long=round(fl, 3), frac_short=round(fs, 3)))
    D = pd.DataFrame(rows); D.to_csv(os.path.join(OUT, "diagnostic_snps.csv"), index=False)
    print(f"assembly lines: {nl} long, {ns} short; diagnostic SNPs: {len(D)} of {len(reg)}")
    if len(D) == 0: sys.exit()
    print("diagnostic SNP span:", D.pos.min(), "-", D.pos.max())
    X = G(np.array([bim.index[(bim.chr == 9) & (bim.pos == p)][0] for p in D.pos]))
    a1 = np.array([bim.a1[bim.index[(bim.chr == 9) & (bim.pos == p)][0]] for p in D.pos])
    Lc = np.where((a1 == D.long_allele.values)[:, None], X, 2 - X)   # copies of the long allele
    score = np.nanmean(Lc, 0); called = (~np.isnan(Lc)).sum(0)
    lead, s152 = snp(151742636), snp(151983481)
    cls = np.select([score >= 1.8, score <= 0.2, (score >= 0.8) & (score <= 1.2)], ["long/long", "short/short", "het"], "ambiguous")
    C = pd.DataFrame({"taxa": taxa, "long_score": score.round(3), "n_called": called, "class": cls,
                      "lead_G": lead, "snp152_C": s152}, index=taxa)
    C.to_csv(os.path.join(OUT, "landrace_haplotype_class.csv"), index=False)
    print(C["class"].value_counts().to_string())
    print("r2 score vs lead %.3f, vs 152 Mb %.3f" % (np.corrcoef(score, np.nan_to_num(lead, nan=np.nanmean(lead)))[0, 1] ** 2,
                                                     np.corrcoef(score, np.nan_to_num(s152, nan=np.nanmean(s152)))[0, 1] ** 2))
    mix = C[((C.lead_G == 0) & (C.snp152_C == 2)) | ((C.lead_G == 2) & (C.snp152_C == 0))]
    print(f"\n{len(mix)} landraces with mixed homozygous combinations; their class:"); print(mix["class"].value_counts().to_string())
    mix.to_csv(os.path.join(OUT, "mixed_landraces_class.csv"), index=False)

elif mode == "models":
    from scipy import stats, optimize
    k = np.load(os.path.join(CACHE, "kin.npz")); lam, U = k["lam"], k["U"]; assert list(k["taxa"]) == taxa
    preds = {"lead SNP": snp(151742636), "152 Mb SNP": snp(151983481)}
    fc = os.path.join(OUT, "landrace_haplotype_class.csv")
    if os.path.exists(fc): preds = {"haplotype class score": pd.read_csv(fc).set_index("taxa").loc[taxa, "long_score"].values, **preds}
    V = ["N_0_5", "N_5_15", "elevation", "aridity", "pH", "CEC", "latitude", "longitude"]
    def z(v): v = np.asarray(v, float); return (v - np.nanmean(v)) / np.nanstd(v)
    def lmm(y, Xc, x):
        """y ~ Xc + x + g, g ~ N(0, s2 K); ML delta on the null model, then GLS Wald test of x."""
        m = ~np.isnan(y) & ~np.isnan(x) & ~np.isnan(Xc).any(1)
        assert m.all(), "missing values"
        yt, Ct, xt = U.T @ y, U.T @ Xc, U.T @ x
        def nll(ld):
            d = lam + np.exp(ld); w = 1 / d
            b = np.linalg.solve((Ct * w[:, None]).T @ Ct, (Ct * w[:, None]).T @ yt); r = yt - Ct @ b
            s2 = (r * r * w).sum() / len(y); return 0.5 * (len(y) * np.log(s2) + np.log(d).sum())
        ld = optimize.minimize_scalar(nll, bounds=(-8, 8), method="bounded").x
        w = 1 / (lam + np.exp(ld)); X = np.column_stack([Ct, xt])
        XtWX = (X * w[:, None]).T @ X; b = np.linalg.solve(XtWX, (X * w[:, None]).T @ yt); r = yt - X @ b
        s2 = (r * r * w).sum() / (len(y) - X.shape[1]); se = np.sqrt(np.diag(np.linalg.inv(XtWX)) * s2)
        t = b[-1] / se[-1]; p = 2 * stats.t.sf(abs(t), len(y) - X.shape[1])
        return b[-1], t, p, np.exp(ld)
    rows = []; one = np.ones((len(taxa), 1))
    for pn, x in preds.items():
        x = np.where(np.isnan(x), np.nanmean(x), x)
        for v in V:
            b, t, p, dl = lmm(z(E[v]), one, z(x)); rows.append(dict(predictor=pn, response=v, covariates="none", beta_sd=round(b, 3), t=round(t, 2), P=p, delta=round(dl, 3)))
        for v, others in [("N_0_5", ["elevation", "aridity", "pH", "CEC", "latitude", "longitude"]),
                          ("elevation", ["N_0_5", "aridity", "pH", "CEC", "latitude", "longitude"]),
                          ("aridity", ["N_0_5", "elevation", "pH", "CEC", "latitude", "longitude"])]:
            Xc = np.column_stack([one] + [z(E[o]) for o in others])
            b, t, p, dl = lmm(z(E[v]), Xc, z(x)); rows.append(dict(predictor=pn, response=v, covariates="other site variables", beta_sd=round(b, 3), t=round(t, 2), P=p, delta=round(dl, 3)))
    R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "haplotype_env_mixed_models.csv"), index=False)
    R["P"] = R.P.map(lambda v: f"{v:.2e}")
    pd.set_option("display.width", 200); print(R.to_string(index=False))

if mode == "regional":
    """Regional scan of the GBS SNPs in chr9:150.7-153.0 Mb for soil N (and aridity, pH):
    (a) plain mixed model, (b) conditioned on the 152 Mb SNP (the only clean marker of the
    interval version in landraces), (c) within 152 Mb SNP homozygotes (CC = long, TT = short),
    each with its own kinship. delta estimated once per model on the null (P3D)."""
    from scipy import stats, optimize
    kz = np.load(os.path.join(CACHE, "kin.npz")); lamA, UA = kz["lam"], kz["U"]
    Kfull = (UA * lamA) @ UA.T
    reg = bim[(bim.chr == 9) & (bim.pos >= 150_700_000) & (bim.pos <= 153_000_000)]
    Xr = G(reg.index.values); f = np.nanmean(Xr, 1) / 2
    ok = (np.minimum(f, 1 - f) >= 0.01); Xr, reg = Xr[ok], reg[ok]
    Xr = np.where(np.isnan(Xr), np.nanmean(Xr, 1, keepdims=True), Xr)
    c152 = snp(151983481)
    def z(v): v = np.asarray(v, float); return (v - v.mean()) / v.std()
    def scan(y, cov, sel):
        K = Kfull[np.ix_(sel, sel)]; lam, U = np.linalg.eigh(K); lam = np.clip(lam, 0, None)
        y = y[sel]; C = np.column_stack([np.ones(sel.sum())] + [c[sel] for c in cov])
        yt, Ct = U.T @ y, U.T @ C
        def nll(ld):
            d = lam + np.exp(ld); w = 1 / d; b = np.linalg.solve((Ct * w[:, None]).T @ Ct, (Ct * w[:, None]).T @ yt); r = yt - Ct @ b
            return 0.5 * (len(y) * np.log((r * r * w).sum() / len(y)) + np.log(d).sum())
        w = 1 / (lam + np.exp(optimize.minimize_scalar(nll, bounds=(-8, 8), method="bounded").x))
        Xt = U.T @ Xr[:, sel].T; P = np.full(Xr.shape[0], np.nan)
        CtW = Ct * w[:, None]; A = CtW.T @ Ct
        for j in range(Xr.shape[0]):
            x = Xt[:, j]
            if Xr[j, sel].std() == 0: continue
            X = np.column_stack([Ct, x]); XtWX = (X * w[:, None]).T @ X
            try: b = np.linalg.solve(XtWX, (X * w[:, None]).T @ yt)
            except np.linalg.LinAlgError: continue
            r = yt - X @ b; s2 = (r * r * w).sum() / (len(y) - X.shape[1]); se = np.sqrt(np.linalg.inv(XtWX)[-1, -1] * s2)
            P[j] = 2 * stats.t.sf(abs(b[-1] / se), len(y) - X.shape[1])
        return P
    allsel = np.ones(len(taxa), bool)
    out = pd.DataFrame({"pos": reg.pos.values, "snp": reg.id.values})
    for v in ["N_0_5", "aridity", "pH"]:
        y = z(E[v])
        out[f"{v}_plain"] = scan(y, [], allsel)
        out[f"{v}_cond152"] = scan(y, [c152], allsel)
    y = z(E["N_0_5"])
    out["N_0_5_within_CC_long"] = scan(y, [], c152 == 2)
    out["N_0_5_within_TT_short"] = scan(y, [], c152 == 0)
    out.to_csv(os.path.join(OUT, "regional_scan_conditional.csv"), index=False)
    print("SNPs", len(out), "| landraces CC", int((c152 == 2).sum()), "TT", int((c152 == 0).sum()))
    for c in out.columns[2:]:
        i = out[c].idxmin(); print(f"{c:28s} min P {out[c].min():.2e} at {out.pos[i]:,}  | SNPs P<1e-4: {(out[c] < 1e-4).sum()}")
    for p in (151742636, 151983481, 151737489):
        r = out[out.pos == p]
        if len(r): print(p, {c: f"{r[c].values[0]:.1e}" for c in out.columns[2:]})
