# Summary figure, triangle layout (final/main_figs/Fig_summary_triangle.{svg,png,pdf}).
# Three analyses sit at the corners (GB grain FAA GWAS, landrace soil-N GWAS, IC/J F_ST scans).
# The L-homoserine pathway sits in the centre. A small triangle glyph next to each gene ID
# shows which analysis found it (filled corner). Other candidates are listed outside each corner.
# Writes an editable SVG (all text kept as text) plus PNG and PDF renders.
# Usage: python3 make_triangle.py <out_prefix>
import math, sys
import cairosvg

OUT = sys.argv[1]
W, H = 1800, 1800
FONT = "Arial, Helvetica, Inter, sans-serif"
GREY = "#8a8a8a"
el = []

def tw(s, size):  # approximate Arial advance widths
    narrow = set(" ,.;:'|iIlj()[]-ft")
    wide = set("MWmw")
    w = 0
    for ch in s:
        w += 0.28 if ch in narrow else (0.83 if ch in wide else (0.67 if ch.isupper() else 0.556))
    return w * size

def t(x, y, s, size=16, w="normal", anchor="start", fill="#111", style="normal", rot=None):
    r = f' transform="rotate({rot} {x} {y})"' if rot is not None else ""
    el.append(f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" font-weight="{w}" '
              f'font-style="{style}" text-anchor="{anchor}" fill="{fill}"{r}>{s}</text>')

def rect(x, y, w, h, fill="white", stroke="#111", sw=1.5, r=10, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    el.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')

def line(x1, y1, x2, y2, col="#111", sw=2.2, dash=None, head=False, marker="ah"):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = f' marker-end="url(#{marker})"' if head else ""
    el.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{col}" stroke-width="{sw}"{d}{m}/>')

def node(x, y, s, w=None, faded=False):
    w = w or (len(s) * 10.5 + 34)
    col = GREY if faded else "#111"
    rect(x - w / 2, y - 22, w, 44, fill="white", stroke=col, sw=2, r=22)
    t(x, y + 6, s, 18, "bold", "middle", col)
    return (x - w / 2, x + w / 2)

# ---- the population glyph: small triangle, filled corner(s) = analysis that found the gene ----
GS, GH = 21, 18  # glyph width / height

def glyph(x, y, corners, scale=1.0):
    """x,y = bottom-left of glyph. corners subset of {'GB','SN','ICJ'}. returns width"""
    s, h = GS * scale, GH * scale
    V = {"GB": (x + s / 2, y - h), "SN": (x, y), "ICJ": (x + s, y)}
    pts = " ".join(f"{px:.1f},{py:.1f}" for px, py in V.values())
    el.append(f'<polygon points="{pts}" fill="white" stroke="#111" stroke-width="{1.3*scale:.2f}" stroke-linejoin="round"/>')
    for c in corners:
        vx, vy = V[c]
        o = [V[k] for k in V if k != c]
        p1 = (vx + (o[0][0] - vx) * 0.55, vy + (o[0][1] - vy) * 0.55)
        p2 = (vx + (o[1][0] - vx) * 0.55, vy + (o[1][1] - vy) * 0.55)
        el.append(f'<polygon points="{vx:.1f},{vy:.1f} {p1[0]:.1f},{p1[1]:.1f} {p2[0]:.1f},{p2[1]:.1f}" fill="#111" stroke="#111" stroke-width="{1.3*scale:.2f}" stroke-linejoin="round"/>')
    return s

def item_width(items, size=14):
    wsum = 0
    for txt, corners, tag in items:
        wsum += tw(txt, size) + 6 + GS + (tw(tag, 12) + 4 if tag else 0) + 14
    return wsum - 14

def gline(x, y, items, anchor="start", size=14, col="#444"):
    """items = [(gene id text, corners, tag)] ; tag = 'IC' or 'J' printed small after glyph"""
    if anchor == "end":
        x = x - item_width(items, size)
    elif anchor == "middle":
        x = x - item_width(items, size) / 2
    for txt, corners, tag in items:
        t(x, y, txt, size, "normal", "start", col)
        x += tw(txt, size) + 6
        x += glyph(x, y + 1, corners)
        if tag:
            t(x + 3, y, tag, 12, "bold", "start", "#111")
            x += tw(tag, 12) + 4
        x += 14

def enzyme(x, y, name, items, anchor="start"):
    t(x, y, name, 16, "bold", anchor)
    gline(x, y + 22, items, anchor)

# ---------------------------------------------------------------- canvas
el.append('<defs>'
          '<marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#111"/></marker>'
          '<marker id="ahg" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="#8a8a8a"/></marker>'
          '</defs>')
el.append(f'<rect width="{W}" height="{H}" fill="white"/>')

# ---------------------------------------------------------------- big triangle
A, B, C = (900, 250), (130, 1470), (1670, 1470)
el.append(f'<polygon points="{A[0]},{A[1]} {B[0]},{B[1]} {C[0]},{C[1]}" fill="#fafafa" stroke="#b5b5b5" stroke-width="3" stroke-linejoin="round"/>')

def edge_label(P, Q, frac, text, offset, size=15):
    x = P[0] + (Q[0] - P[0]) * frac
    y = P[1] + (Q[1] - P[1]) * frac
    ang = math.degrees(math.atan2(Q[1] - P[1], Q[0] - P[0]))
    if ang > 90 or ang < -90:
        ang += 180
    # outward normal
    nx, ny = (Q[1] - P[1]), -(Q[0] - P[0])
    n = math.hypot(nx, ny); nx, ny = nx / n, ny / n
    cxm, cym = (A[0] + B[0] + C[0]) / 3, (A[1] + B[1] + C[1]) / 3
    if (x + nx - cxm) ** 2 + (y + ny - cym) ** 2 < (x - cxm) ** 2 + (y - cym) ** 2:
        nx, ny = -nx, -ny
    x += nx * offset; y += ny * offset
    t(x, y, text, size, "normal", "middle", "#333", "italic", rot=ang)

edge_label(A, B, 0.42, "Shared gene, SAM synthase (354640)", 16)
edge_label(A, C, 0.42, "Shared gene, class I Gln amidotransferase (217220)", 16)
edge_label(B, C, 0.50, "Shared gene, Met synthase (052750)", 30)

# ---------------------------------------------------------------- corner icons
R = 92
def corner_circle(P):
    el.append(f'<circle cx="{P[0]}" cy="{P[1]}" r="{R}" fill="white" stroke="#111" stroke-width="2.6"/>')

# GB, a maize kernel with free amino acids
corner_circle(A)
cx, cy = A
el.append(f'<path d="M{cx-40},{cy-22} Q{cx-42},{cy-62} {cx},{cy-62} Q{cx+42},{cy-62} {cx+40},{cy-22} '
          f'L{cx+13},{cy+54} Q{cx},{cy+64} {cx-13},{cy+54} Z" fill="white" stroke="#111" stroke-width="2.4" stroke-linejoin="round"/>')
el.append(f'<path d="M{cx-11},{cy+46} Q{cx-21},{cy+2} {cx},{cy-14} Q{cx+21},{cy+2} {cx+11},{cy+46} Q{cx},{cy+53} {cx-11},{cy+46} Z" fill="#d9d9d9" stroke="#111" stroke-width="1.6"/>')
for dx, dy in [(-24, -38), (-6, -46), (14, -40), (28, -26), (-28, -18), (22, -8), (-18, 8)]:
    el.append(f'<circle cx="{cx+dx}" cy="{cy+dy}" r="4" fill="#111"/>')

# SN, soil profile with a seedling
corner_circle(B)
cx, cy = B
el.append(f'<line x1="{cx-62}" y1="{cy+4}" x2="{cx+62}" y2="{cy+4}" stroke="#111" stroke-width="2.6"/>')
el.append(f'<path d="M{cx-60},{cy+4} L{cx+60},{cy+4} L{cx+48},{cy+52} Q{cx},{cy+72} {cx-48},{cy+52} Z" fill="#e6e6e6" stroke="none"/>')
el.append(f'<line x1="{cx-56}" y1="{cy+26}" x2="{cx+56}" y2="{cy+26}" stroke="#777" stroke-width="1.4" stroke-dasharray="4 4"/>')
el.append(f'<line x1="{cx}" y1="{cy+4}" x2="{cx}" y2="{cy-52}" stroke="#111" stroke-width="2.6"/>')
el.append(f'<path d="M{cx},{cy-22} Q{cx-34},{cy-30} {cx-42},{cy-56} Q{cx-12},{cy-50} {cx},{cy-22} Z" fill="white" stroke="#111" stroke-width="2"/>')
el.append(f'<path d="M{cx},{cy-34} Q{cx+32},{cy-42} {cx+40},{cy-68} Q{cx+10},{cy-62} {cx},{cy-34} Z" fill="white" stroke="#111" stroke-width="2"/>')
for (x2, y2) in [(-18, 40), (-5, 48), (8, 44), (20, 34)]:
    el.append(f'<path d="M{cx},{cy+4} Q{cx+x2*0.4},{cy+20} {cx+x2},{cy+y2}" fill="none" stroke="#111" stroke-width="1.6"/>')
for dx, dy in [(-42, 16), (36, 16), (-34, 40), (40, 40)]:
    el.append(f'<circle cx="{cx+dx}" cy="{cy+dy}" r="4" fill="#111"/>')

# IC/J, two ears (Cycle 0 to Cycle 14)
corner_circle(C)
cx, cy = C
for ex, lab in [(cx - 27, "IC"), (cx + 27, "J")]:
    ey = cy - 12
    el.append(f'<path d="M{ex},{ey+46} Q{ex-26},{ey+30} {ex-22},{ey-6} Q{ex-14},{ey+24} {ex},{ey+40} Z" fill="#d9d9d9" stroke="#111" stroke-width="1.6"/>')
    el.append(f'<ellipse cx="{ex}" cy="{ey}" rx="16" ry="46" fill="white" stroke="#111" stroke-width="2.2"/>')
    for k in range(-4, 5):
        yy = ey + k * 9
        half = 16 * math.sqrt(max(0, 1 - (k * 9 / 46) ** 2))
        el.append(f'<line x1="{ex-half:.1f}" y1="{yy}" x2="{ex+half:.1f}" y2="{yy}" stroke="#111" stroke-width="0.9"/>')
    el.append(f'<line x1="{ex}" y1="{ey-46}" x2="{ex}" y2="{ey+46}" stroke="#111" stroke-width="0.9"/>')
    t(ex, cy + 62, lab, 14, "bold", "middle")

# ---------------------------------------------------------------- corner text
t(900, 58, "Goodman–Buckler panel", 22, "bold", "middle")
t(900, 86, "Grain free amino acid GWAS, 280 inbreds, 71 traits", 15, "normal", "middle", "#333")
t(900, 110, "Free pool dominated by Asn and Pro", 15, "bold", "middle")

t(30, 1600, "Maize landraces", 22, "bold")
t(30, 1628, "Soil N GWAS, 3,107 georeferenced accessions", 15, "normal", "start", "#333")
t(30, 1652, "12 SNPs pass both specificity tests", 15, "bold")

t(1770, 1600, "Indian Chief and Jarvis", 22, "bold", "end")
t(1770 - tw("FST scan, Cycle 0 vs Cycle 14", 15), 1628, "F<tspan baseline-shift='sub' font-size='11'>ST</tspan> scan, Cycle 0 vs Cycle 14", 15, "normal", "start", "#333")
t(1770, 1652, "IC free AA and Pro up, Asp and Glu down", 15, "bold", "end")
t(1770, 1676, "J Asp fraction up (nominal only)", 15, "bold", "end")

# ---------------------------------------------------------------- centre pathway
X = 900
OAA, ASP, ASPP, ASA, HSER = (X, 640), (X, 750), (X, 865), (X, 980), (X, 1095)
ASN, GLN, LYS = (1150, 750), (1150, 875), (1300, 980)
THR, ILE, MET, SAM, SMM = (760, 1215), (760, 1335), (1040, 1215), (1235, 1215), (1040, 1335)

# shaded L-homoserine box
rect(668, 788, 345, 334, fill="#ececec", stroke="#777", sw=1.5, r=14, dash="7 5")
t(690, 955, "L-homoserine biosynthesis", 16, "bold", "middle", "#333", rot=-90)

# arrows (drawn first so nodes sit on top)
line(X, 662, X, 726, head=True)                       # OAA -> Asp
line(X, 772, X, 841, sw=3.4, head=True)               # Asp -> Asp-4-P
line(X, 887, X, 956, sw=3.4, head=True)               # -> ASA
line(X, 1002, X, 1071, sw=3.4, head=True)             # -> Hser
line(942, 750, 1110, 750, head=True)                  # Asp -> Asn
line(ASN[0], 853, ASN[0], 775, col=GREY, dash="6 5", head=True, marker="ahg")   # Gln amide N -> Asn
line(1122, 730, 978, 650, dash="6 5", head=True)      # Asn -> OAA via omega-amidase
line(1002, 980, 1262, 980, head=True)                 # ASA -> Lys
line(880, 1117, 778, 1191, head=True)                 # Hser -> Thr
line(920, 1117, 1020, 1191, head=True)                # Hser -> Met
line(THR[0], 1237, THR[0], 1311, head=True)           # Thr -> Ile
line(1077, 1215, 1196, 1215, head=True)               # Met -> SAM
line(MET[0], 1237, MET[0], 1311, head=True)           # Met -> SMM

node(*OAA, "Oxaloacetate", faded=True)
node(*ASP, "Asp")
node(*ASPP, "Aspartyl-4-P")
node(*ASA, "Asp-semialdehyde", w=200)
node(*HSER, "Homoserine")
node(*ASN, "Asn")
node(*GLN, "Gln", faded=True)
node(*LYS, "Lys")
node(*THR, "Thr"); node(*ILE, "Ile"); node(*MET, "Met"); node(*SAM, "SAM")
node(*SMM, "S-methyl-Met", w=150)

# enzyme labels
enzyme(884, 683, "Asp aminotransferase", [("152450", ["ICJ"], "J")], "end")
enzyme(884, 812, "Asp kinase", [("094670", ["ICJ"], "J"), ("220040", ["SN"], None)], "end")
enzyme(884, 915, "Asp-semialdehyde DH", [("214360", ["ICJ"], "IC")], "end")
enzyme(884, 1030, "Homoserine DH", [("405410", ["GB"], None)], "end")
t(1026, 740, "Asn synthetase domain", 14, "bold", "middle")
gline(1026, 775, [("036060, 036070", ["GB"], None)], "middle")
t(1002, 608, "ω-amidase", 16, "bold")
gline(1002, 630, [("005710", ["SN"], None)])
t(1166, 820, "amide N", 14, "normal", "start", GREY, "italic")
t(1132, 970, "several steps", 13, "normal", "middle", "#666", "italic")
t(1020, 1006, "DHDPS, DHDPR, DAPDC", 16, "bold")
gline(1020, 1028, [("389290, 176170, 406820", ["GB"], None)])
enzyme(772, 1142, "Thr synthase", [("294790", ["ICJ"], "IC")], "end")
enzyme(1030, 1142, "Met synthase", [("052750", ["SN", "ICJ"], "IC")], "start")
t(1136, 1199, "SAM synthase", 15, "bold", "middle")
gline(1136, 1243, [("354640", ["GB", "SN"], None)], "middle")
enzyme(1056, 1282, "Met S-methyltransferase", [("356170", ["SN"], None)], "start")

t(900, 505, "Aspartate-family", 18, "bold", "middle", "#555")
t(900, 529, "amino acid metabolism", 18, "bold", "middle", "#555")

# ---------------------------------------------------------------- outside lists
def listbox(x, y, w, title, groups, h=None):
    rows = sum(1 + len(g[1]) for g in groups)
    h = h or 52 + rows * 25 + len(groups) * 6
    rect(x, y, w, h, fill="white", stroke="#999", sw=1.3)
    t(x + 16, y + 30, title, 17, "bold")
    yy = y + 58
    for gname, items in groups:
        t(x + 16, yy, gname, 14, "bold", "start", "#555", "italic")
        yy += 24
        for txt, corners, tag in items:
            gx = x + 18
            glyph(gx, yy + 1, corners)
            if tag:
                t(gx + GS + 3, yy, tag, 12, "bold")
                gx += tw(tag, 12) + 4
            t(gx + GS + 10, yy, txt, 14, "normal", "start", "#222")
            yy += 24
        yy += 6
    return h

listbox(1095, 150, 400, "Other grain FAA GWAS candidates", [
    ("Protein turnover and ER membrane", [("Peptidase A1 (190800), 8 traits", ["GB"], None),
                                          ("EMC4 (394360), rhomboid-like 19 (394380), 9 traits", ["GB"], None)]),
    ("N assimilation", [("GOGAT (329710, 360480), GS (059870)", ["GB"], None)]),
    ("Transport", [("NPF4.6 (277160)", ["GB"], None),
                   ("UmamiTs (379220, 207630, 352400)", ["GB"], None),
                   ("AVT6C (368680), AVT3A-like (429100)", ["GB"], None)]),
])

listbox(1395, 640, 385, "Other IC and J candidates", [
    ("Protein turnover", [("SAG39 cysteine protease (214780)", ["ICJ"], "IC"),
                          ("ATG18a-like (352300), ATG18g (287690)", ["ICJ"], "IC"),
                          ("NPL4-like 1 (146560)", ["ICJ"], "J")]),
    ("N assimilation", [("GDH (060520, 162300), GS (162900)", ["ICJ"], "IC")]),
    ("Transport", [("AMT2-like (214890), NPF3.1 (287950)", ["ICJ"], "IC"),
                   ("LHT8-like (132750), ProT-like (159580)", ["ICJ"], "J")]),
])

listbox(20, 640, 385, "Other soil N candidates", [
    ("Within or near the 12 SNPs", [("Syndetin (022740)", ["SN"], None),
                                    ("RLT1 homeobox-DDT (359990)", ["SN"], None),
                                    ("TTL1-like zinc finger (130210)", ["SN"], None),
                                    ("Protein kinase (226000)", ["SN"], None),
                                    ("HB transcription factor (265710)", ["SN"], None),
                                    ("Transcription repressor (001920)", ["SN"], None)]),
])

# ---------------------------------------------------------------- legend (top left) and QTL tag
rect(40, 150, 400, 250, fill="white", stroke="#999", sw=1.3)
t(56, 180, "How to read the gene marks", 17, "bold")
gx, gy = 100, 300
glyph(gx, gy, [], scale=4.2)
t(gx + GS * 4.2 / 2, gy - GH * 4.2 - 10, "GB", 14, "bold", "middle")
t(gx - 8, gy + 20, "SN", 14, "bold", "middle")
t(gx + GS * 4.2 + 10, gy + 20, "IC / J", 14, "bold", "middle")
t(240, 230, "A filled corner shows which", 14, "normal", "start", "#222")
t(240, 250, "analysis found the gene.", 14, "normal", "start", "#222")
x0 = 250
x0 += glyph(x0, 290, ["GB"]) + 8; t(x0, 289, "GB only", 13, "normal", "start", "#222")
x0 = 250
x0 += glyph(x0, 318, ["SN", "ICJ"]) + 8; t(x0, 317, "SN and IC", 13, "normal", "start", "#222")
t(56, 385, "IC or J is printed after marks at the IC / J corner.", 13, "normal", "start", "#555", "italic")

rect(470, 300, 300, 62, fill="white", stroke="#999", sw=1.3, dash="5 4")
t(486, 326, "B64 × CML14 biparental cross", 14, "bold")
t(486, 348, "Shared Asp and Pro QTL on chr 3", 14, "normal", "start", "#222")

# ---------------------------------------------------------------- footnotes
fy = 1716
t(30, fy, "Genes are positional candidates, within 25 kb of a significant GWAS SNP or in a top-5% F<tspan baseline-shift='sub' font-size='10'>ST</tspan> window. Numbers are Zm00001eb gene IDs (B73 RefGen_v5).", 15, fill="#333")
t(30, fy + 26, "Shaded box, the CornCyc L-homoserine biosynthesis pathway, the only pathway with a candidate in all four gene sets and the top-ranked pathway in the recurrence test.", 15, fill="#333")
t(30, fy + 52, "Thick arrows, L-homoserine biosynthesis. Dashed arrows, amide N from Gln to Asn, and the ω-amidase route returning Asn carbon to oxaloacetate. Grey, context not found in any analysis.", 15, fill="#333")

svg = f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">' + "\n".join(el) + "</svg>"
open(OUT + ".svg", "w").write(svg)
cairosvg.svg2png(bytestring=svg.encode(), write_to=OUT + ".png", output_width=W * 2)
cairosvg.svg2pdf(bytestring=svg.encode(), write_to=OUT + ".pdf")
print("ok")
