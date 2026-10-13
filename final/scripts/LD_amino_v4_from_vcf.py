#!/usr/bin/env python3
"""
LD intervals for the Goodman-Buckler amino acid GWAS hits, computed from the
HapMap 3.2.1 282-panel VCFs (hmp321_282_agpv4_chr*.vcf.gz; POS = B73 AGPv4, which
is the coordinate system of the GWAS 'Pos' column).

For each lead SNP: r^2 (squared correlation of allele dosages, missing genotypes
dropped pairwise) with every biallelic SNP within +/-WINDOW bp that has MAF >= MINMAF
(INFO MAF field).

LD interval. With ~12 SNPs/kb, single distant SNPs reach r^2 >= 0.4 by chance or
through population structure (0.6% of SNPs 100-250 kb from a lead do), so the outermost
partner is not a usable boundary. Instead the interval is grown outward from the lead
in 5-kb steps, on each side separately, and a step is added while >= 10% of its SNPs
(>= 5 SNPs per step) are in LD (r^2 >= R2) with the lead; growth stops at the first
step that fails. Steps with < 5 SNPs are passed over without stopping. 10% is ~16x the
distant background rate, so a step passes by chance with P < 1e-5. The outermost-
partner interval is also written, for reference only. Coordinates stay AGPv4; they are
lifted to v5 afterwards (LD_amino_liftover_genes.py).

The VCFs are BGZF without an index; random access uses the BGZF block structure
(block offsets cached per chromosome, binary search on the first record of a block).
The run is resumable and stops after --budget seconds, so it can be called
repeatedly: progress is kept in <out>/done.txt and results appended to
<out>/ld_intervals_v4.tsv.

usage: LD_amino_v4_from_vcf.py <vcf_dir> <leads.tsv chr pos id> <out_dir>
                               [--r2 0.4] [--window 250000] [--minmaf 0.05] [--budget 150]
"""
import sys, os, struct, zlib, time, argparse
import numpy as np
ap = argparse.ArgumentParser()
ap.add_argument("vcfdir"); ap.add_argument("leads"); ap.add_argument("out")
ap.add_argument("--r2", type=float, default=0.4); ap.add_argument("--window", type=int, default=250000)
ap.add_argument("--minmaf", type=float, default=0.05); ap.add_argument("--budget", type=float, default=150)
a = ap.parse_args(); T0 = time.time(); os.makedirs(a.out, exist_ok=True)
MAGIC = b"\x1f\x8b\x08\x04"

def block_offsets(fn, chrom):
    cache = os.path.join(a.out, f"blocks_chr{chrom}.npy")
    if os.path.exists(cache): return np.load(cache)
    offs = []; off = 0; f = open(fn, "rb"); size = os.path.getsize(fn)
    while off < size:
        f.seek(off); hd = f.read(18)
        if len(hd) < 18 or hd[:4] != MAGIC: break
        offs.append(off); off += struct.unpack("<H", hd[16:18])[0] + 1
    offs = np.array(offs, dtype=np.int64); np.save(cache, offs); return offs

def read_block(f, off):
    f.seek(off); hd = f.read(18); bs = struct.unpack("<H", hd[16:18])[0] + 1
    f.seek(off); return zlib.decompress(f.read(bs)[18:-8], -15)

def first_pos(f, off):
    d = read_block(f, off); i = d.find(b"\n")
    while i >= 0:
        j = d.find(b"\n", i + 1)
        if j < 0: return None
        line = d[i + 1:j]
        if line and not line.startswith(b"#"):
            p = line.split(b"\t", 2)
            try: return int(p[1])
            except Exception: return None
        i = j
    return None

def seek_block(f, offs, target):
    lo, hi = 0, len(offs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2; p = first_pos(f, offs[mid])
        if p is None or p < target: lo = mid
        else: hi = mid
    return max(lo - 1, 0)

GTMAP = {b"0/0": 0, b"0|0": 0, b"0/1": 1, b"0|1": 1, b"1|0": 1, b"1/0": 1, b"1/1": 2, b"1|1": 2}
def stream(f, offs, bi, start, end):
    """yield (pos, dosage array) for SNPs with start <= pos <= end"""
    buf = b""; first = bi > 0
    while bi < len(offs):
        buf += read_block(f, offs[bi]); bi += 1
        if first:                       # drop the partial record at a mid-file block start
            i = buf.find(b"\n")
            if i < 0: continue
            buf = buf[i + 1:]; first = False
        *lines, buf = buf.split(b"\n")
        for l in lines:
            if not l or l[:1] == b"#": continue
            p = l.split(b"\t", 8)
            pos = int(p[1])
            if pos < start: continue
            if pos > end: return
            if len(p[3]) != 1 or len(p[4]) != 1: continue
            k = p[7].find(b";MAF=")
            if k >= 0:
                e = p[7].find(b";", k + 5); maf = float(p[7][k + 5:e if e > 0 else None])
                if maf < a.minmaf: continue
            g = p[8].split(b"\t")[1:]
            yield pos, np.array([GTMAP.get(x[:3], -1) for x in g], dtype=np.int8)

leads = [l.split() for l in open(a.leads) if l.strip()]
done = set(open(os.path.join(a.out, "done.txt")).read().split()) if os.path.exists(os.path.join(a.out, "done.txt")) else set()
fo = open(os.path.join(a.out, "ld_intervals_v4.tsv"), "a")
if os.path.getsize(os.path.join(a.out, "ld_intervals_v4.tsv")) == 0:
    fo.write("lead\tchr\tpos_v4\tfound\tn_window_snps\tn_partners\tmax_r2\touter_start\touter_end\tinterval_start\tinterval_end\n")
fd = open(os.path.join(a.out, "done.txt"), "a")
W = a.window; n_new = 0; STEP = 5000; FRAC = 0.10
by_chr = {}
for c, p, i in leads: by_chr.setdefault(c, []).append((int(p), i))
for c in sorted(by_chr, key=int):
    todo = sorted(x for x in by_chr[c] if x[1] not in done)
    if not todo: continue
    fn = os.path.join(a.vcfdir, f"hmp321_282_agpv4_chr{c}.vcf.gz"); f = open(fn, "rb")
    offs = block_offsets(fn, c)
    k = 0
    while k < len(todo):
        if time.time() - T0 > a.budget:
            print(f"budget reached: {n_new} new leads; {len(done)} done of {len(leads)}"); sys.exit(0)
        # group consecutive leads whose windows overlap -> one streaming pass
        grp = [todo[k]]
        while k + len(grp) < len(todo) and todo[k + len(grp)][0] - grp[-1][0] <= 2 * W and len(grp) < 40:
            grp.append(todo[k + len(grp)])
        s, e = grp[0][0] - W, grp[-1][0] + W
        P, G = [], []
        for pos, g in stream(f, offs, seek_block(f, offs, s), s, e):
            P.append(pos); G.append(g)
        P = np.array(P); G = np.array(G, dtype=np.float32) if G else np.zeros((0, 1), np.float32)
        G[G < 0] = np.nan
        for lp, lid in grp:
            sel = np.where((P >= lp - W) & (P <= lp + W))[0]
            j = np.where(P[sel] == lp)[0]
            if len(j) == 0:
                fo.write(f"{lid}\t{c}\t{lp}\tFALSE\t{len(sel)}\tNA\tNA\tNA\tNA\tNA\tNA\n")
            else:
                X = G[sel]; x = X[j[0]]; ok = ~np.isnan(X) & ~np.isnan(x)
                xm = np.where(ok, x, 0); Xm = np.where(ok, X, 0); m = ok.sum(1)
                with np.errstate(invalid="ignore", divide="ignore"):
                    mx, mX = xm.sum(1) / m, Xm.sum(1) / m
                    cov = (xm * Xm).sum(1) / m - mx * mX
                    r2 = cov ** 2 / (((xm * xm).sum(1) / m - mx ** 2) * ((Xm * Xm).sum(1) / m - mX ** 2))
                r2[j[0]] = np.nan
                wp = P[sel]; part = np.where(r2 >= a.r2)[0]
                ps = np.r_[wp[part], lp]; bs_, be_ = ps.min(), ps.max()
                inld = r2 >= a.r2
                def grow(sign):
                    edge, k = 0, 1
                    while k * STEP <= W:
                        d = (wp - lp) * sign
                        sb = (d > (k - 1) * STEP) & (d <= k * STEP)
                        nb_ = sb.sum()
                        if nb_ >= 5:
                            if inld[sb].mean() >= FRAC: edge = k * STEP
                            else: break
                        k += 1
                    return edge
                ia, ib = max(1, lp - grow(-1)), lp + grow(1)
                mr = np.nanmax(r2) if np.isfinite(r2).any() else float("nan")
                fo.write(f"{lid}\t{c}\t{lp}\tTRUE\t{len(sel)}\t{len(part)}\t{mr:.4f}\t{bs_}\t{be_}\t{ia}\t{ib}\n")
            fd.write(lid + "\n"); done.add(lid); n_new += 1
        fo.flush(); fd.flush()
        k += len(grp)
print(f"finished: {n_new} new leads; {len(done)} done of {len(leads)}")
