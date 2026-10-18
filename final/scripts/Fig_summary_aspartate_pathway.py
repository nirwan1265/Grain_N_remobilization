# Summary figure (final/main_figs/Fig_summary_aspartate_pathway.{svg,png,pdf}): the aspartate pathway with each population's candidate genes on their steps.
# Writes an editable SVG (text kept as text) and PNG/PDF renders.
import cairosvg, sys
OUT = sys.argv[1]
W, H = 1800, 1180
FONT = "Arial, Helvetica, Inter, sans-serif"
el = []
def t(x, y, s, size=17, w="normal", anchor="start", fill="#111", style="normal"):
    el.append(f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" font-weight="{w}" '
              f'font-style="{style}" text-anchor="{anchor}" fill="{fill}">{s}</text>')
def rect(x, y, w, h, fill="white", stroke="#111", sw=1.5, r=10, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    el.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')
def arrow(x1, y1, x2, y2, dash=None, sw=2.2, col="#111", head=True):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = ' marker-end="url(#ah)"' if head else ""
    el.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{col}" stroke-width="{sw}"{d}{m}/>')
def node(x, y, s, w=None):
    w = w or (len(s) * 10.5 + 34)
    rect(x - w / 2, y - 22, w, 44, fill="white", sw=2, r=22)
    t(x, y + 6, s, 18, "bold", "middle")
def badge(x, y, s):
    w = len(s) * 11 + 14
    rect(x, y - 15, w, 22, fill="#111", stroke="#111", sw=1, r=5)
    t(x + w / 2, y + 1, s, 14, "bold", "middle", "white")
    return x + w + 5
def enzyme(x, y, name, genes, badges, anchor="start"):
    # name on one line; gene IDs + population badges on the next
    t(x, y, name, 16, "bold", anchor, "#111")
    gx = x
    t(gx, y + 22, genes, 14, "normal", anchor, "#555")
    bx = x + len(genes) * 7.6 + 8 if anchor == "start" else x + 8
    for b in badges:
        bx = badge(bx, y + 19, b)

el.append('<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
          '<path d="M0,0 L10,5 L0,10 z" fill="#111"/></marker></defs>')
el.append(f'<rect width="{W}" height="{H}" fill="white"/>')

# ---------------- left: populations ----------------
t(30, 48, "Populations and phenotypes", 21, "bold")
pops = [("GB", "Goodman–Buckler panel", "280 inbred lines", "Grain free amino acid GWAS (71 traits)", "Free pool dominated by Asn and Pro"),
        ("SN", "Maize landraces", "3,107 georeferenced accessions", "Environmental GWAS of soil N", "12 SNPs pass specificity tests"),
        ("IC", "Indian Chief, C0 vs C14", "Reciprocal recurrent selection", "F<tspan baseline-shift='sub' font-size='12'>ST</tspan> scan, top-5% windows", "Free AA and Pro up, Asp and Glu down"),
        ("J", "Jarvis, C0 vs C14", "Reciprocal recurrent selection", "F<tspan baseline-shift='sub' font-size='12'>ST</tspan> scan, top-5% windows", "Asp fraction up (nominal only)"),
        ("QTL", "B64 × CML14 F<tspan baseline-shift='sub' font-size='12'>2:3</tspan>", "Biparental cross", "QTL mapping of free Asp and Pro", "Shared Asp/Pro QTL on chr 3")]
y0 = 70
for code, a, b, c, d in pops:
    rect(30, y0, 370, 172, fill="#f4f4f4", stroke="#999", sw=1.2)
    badge(46, y0 + 34, code)
    t(46 + len(code) * 11 + 26, y0 + 35, a, 18, "bold")
    t(46, y0 + 70, b, 15, fill="#444")
    t(46, y0 + 100, c, 15)
    t(46, y0 + 130, d, 15, "bold")
    y0 += 188

# ---------------- centre: pathway ----------------
cx0 = 440
t(cx0, 48, "Candidate genes on aspartate and glutamate metabolism", 21, "bold")
# homoserine band
rect(452, 252, 330, 518, fill="#e9e9e9", stroke="#777", sw=1.5, r=14, dash="7 5")
el.append(f'<text x="476" y="640" transform="rotate(-90 476 640)" font-family="{FONT}" font-size="16" font-weight="bold" text-anchor="middle" fill="#333">L-homoserine biosynthesis</text>')
t(468, 782 - 2, "", 1)
# nodes
OAA, ASP, ASPP, ASA, HSER = (617, 150), (617, 300), (617, 430), (617, 560), (617, 690)
TOG, GLU, GLN = (1110, 150), (1110, 300), (1310, 300)
ASN = (1000, 470)
LYS = (1310, 560)
THR, ILE, MET, SAM = (617, 870), (617, 990), (900, 870), (1180, 870)
t(860, 110, "TCA cycle", 16, "normal", "middle", "#666", "italic")
arrow(705, 150, 1010, 150, dash="2 6", sw=1.6, col="#888", head=False)
for (x, y), s in [(OAA, "Oxaloacetate"), (ASP, "Asp"), (ASPP, "Aspartyl-4-P"), (ASA, "Asp-semialdehyde"),
                  (HSER, "Homoserine"), (TOG, "2-Oxoglutarate"), (GLU, "Glu"), (GLN, "Gln"), (ASN, "Asn"),
                  (LYS, "Lys"), (THR, "Thr"), (ILE, "Ile"), (MET, "Met"), (SAM, "SAM")]:
    node(x, y, s)
# arrows + enzymes
arrow(OAA[0], OAA[1] + 22, ASP[0], ASP[1] - 24)                         # OAA -> Asp
t(600, 222, "Asp aminotransferase", 16, "bold", "end"); t(560, 244, "152450", 14, "normal", "end", "#555"); badge(568, 241, "J")
arrow(TOG[0], TOG[1] + 22, GLU[0], GLU[1] - 24)                         # 2-OG -> Glu
enzyme(1128, 205, "Glu dehydrogenase", "060520, 162300", ["IC"])
arrow(GLU[0] + 30, GLU[1], GLN[0] - 34, GLN[1])                         # Glu -> Gln
t(1135, 345, "GS / GOGAT", 16, "bold"); t(1135, 367, "162900", 14, "normal", "start", "#555"); badge(1191, 364, "IC"); t(1135, 391, "329710, 360480", 14, "normal", "start", "#555"); badge(1135+14*7.6+8, 388, "GB")
arrow(ASP[0], ASP[1] + 22, ASPP[0], ASPP[1] - 24, sw=3)                  # Asp -> Asp-P
enzyme(468, 355, "Asp kinase", "094670", ["J"])
t(468, 399, "220040", 14, "normal", "start", "#555"); badge(468 + 6 * 7.6 + 8, 396, "SN")
arrow(ASPP[0], ASPP[1] + 22, ASA[0], ASA[1] - 24, sw=3)
enzyme(640, 487, "Asp-semialdehyde DH", "214360", ["IC"])
arrow(ASA[0], ASA[1] + 22, HSER[0], HSER[1] - 24, sw=3)
enzyme(640, 617, "Homoserine DH", "405410", ["GB"])
# Asp -> Asn (+ Gln amide donor), Asn -> OAA via omega-amidase
arrow(ASP[0] + 40, ASP[1] + 14, ASN[0] - 30, ASN[1] - 24)
enzyme(800, 330, "Asn synthetase", "036060, 036070", ["GB"])
el.append(f'<polyline points="{GLN[0]},{GLN[1]+22} {GLN[0]},{ASN[1]} {ASN[0]+42},{ASN[1]}" fill="none" stroke="#111" stroke-width="1.6" stroke-dasharray="5 5" marker-end="url(#ah)"/>')
t(1230, 462, "amide N", 14, "normal", "middle", "#555", "italic")
t(1048, 506, "AVT6C transporter", 15, "bold"); t(1048, 528, "368680", 14, "normal", "start", "#555"); badge(1102, 525, "GB")
el.append(f'<polyline points="{ASN[0]+40},{ASN[1]-22} {ASN[0]+40},190 {OAA[0]+92},{OAA[1]+14}" fill="none" stroke="#111" stroke-width="1.8" stroke-dasharray="8 5" marker-end="url(#ah)"/>')
t(1025, 214, "Asn transamination,", 16, "bold", "end"); t(1025, 236, "ω-amidase", 16, "bold", "end"); t(980, 258, "005710", 14, "normal", "end", "#555"); badge(986, 255, "SN")
# ASA -> Lys
arrow(ASA[0] + 90, ASA[1], LYS[0] - 34, LYS[1], sw=2.2)
enzyme(800, 600, "DHDPS, DHDPR, DAP decarboxylase", "389290, 176170, 406820", ["GB"])
t(1180, 582, "several steps", 14, "normal", "middle", "#666", "italic")
# Homoserine -> Thr -> Ile ; Homoserine -> Met -> SAM
arrow(HSER[0], HSER[1] + 22, THR[0], THR[1] - 24)
enzyme(468, 805 + 8, "Thr synthase", "294790", ["IC"])
arrow(THR[0], THR[1] + 22, ILE[0], ILE[1] - 24)
arrow(HSER[0] + 50, HSER[1] + 20, MET[0] - 24, MET[1] - 24)
enzyme(780, 760, "Met synthase", "052750", ["IC", "SN"])
arrow(MET[0] + 30, MET[1], SAM[0] - 36, SAM[1])
enzyme(960, 915, "SAM salvage", "pathway hit by", ["GB", "SN", "IC"])

# ---------------- right: other candidates ----------------
rx = 1400
t(rx, 48, "Other candidate genes", 21, "bold")
t(rx, 72, "(outside the aspartate pathway)", 15, "normal", fill="#555")
groups = [("Protein turnover", [("Peptidase A1 (190800), 8 traits", ["GB"]), ("SAG39 cysteine protease (214780)", ["IC"]),
                               ("ATG18a-like autophagy (352300)", ["IC"]), ("NPL4-like (146560)", ["J"])]),
          ("ER membrane proteins", [("EMC subunit 4 (394360)", ["GB"]), ("Rhomboid-like 19 (394380), 9 traits", ["GB"])]),
          ("Transport", [("NPF4.6, UmamiTs, AVT3A", ["GB"]), ("LHT-like, ProT-like", ["J"]), ("AMT2-like (214890)", ["IC"])])]
gy = 100
for title, items in groups:
    hgt = 52 + 34 * len(items)
    rect(rx, gy, 380, hgt, fill="#f4f4f4", stroke="#999", sw=1.2)
    t(rx + 16, gy + 32, title, 17, "bold")
    iy = gy + 66
    for s, bs in items:
        bx = rx + 16
        for b in bs: bx = badge(bx, iy, b)
        t(bx + 4, iy + 1, s, 14.5)
        iy += 34
    gy += hgt + 22

# ---------------- footnote ----------------
t(30, 1100, "Genes are positional candidates: within 25 kb of a significant GWAS SNP, or in a top-5% F<tspan baseline-shift='sub' font-size='11'>ST</tspan> window. Numbers are Zm00001eb gene IDs (B73 RefGen_v5).", 15, fill="#333")
t(30, 1128, "Shaded box: L-homoserine biosynthesis, the only CornCyc pathway with a candidate gene in all four analyses (permutation P = 0.0012; not significant after correction for 468 pathways).", 15, fill="#333")
t(30, 1156, "Thick arrows, steps of L-homoserine biosynthesis. Dashed arrows, amide-N donation from Gln to Asn synthesis, and Asn transamination with ω-amidase returning carbon to oxaloacetate.", 15, fill="#333")

svg = f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">' + "\n".join(el) + "</svg>"
open(OUT + ".svg", "w").write(svg)
cairosvg.svg2png(bytestring=svg.encode(), write_to=OUT + ".png", output_width=W * 2)
cairosvg.svg2pdf(bytestring=svg.encode(), write_to=OUT + ".pdf")
print("ok")
