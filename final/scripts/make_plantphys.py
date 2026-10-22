"""Build a Plant Physiology (OUP template) version of main.tex without touching it.

Usage: python3 make_plantphys.py <main.tex> <outdir>
Writes <outdir>/PlantPhys_manuscript.tex and <outdir>/PlantPhys_supplement.tex.
Figure paths are rewritten relative to final/plant_physiology/.
"""
import re, sys, os

src, out = sys.argv[1], sys.argv[2]
s = open(src, encoding="utf-8").read()


def balanced(text, start):
    """text[start] == '{'; return index just past the matching brace."""
    depth = 0
    for i in range(start, len(text)):
        c = text[i]
        if c == "\\":
            continue
        if c == "{" and text[i - 1] != "\\":
            depth += 1
        elif c == "}" and text[i - 1] != "\\":
            depth -= 1
            if depth == 0:
                return i + 1
    raise ValueError("unbalanced")


def arg_after(text, cmd, start=0):
    i = text.index(cmd, start) + len(cmd)
    while text[i] != "{":
        i += 1
    j = balanced(text, i)
    return text[i + 1 : j - 1], j


def strip_comments(t):
    return re.sub(r"(?<!\\)%.*", "", t)


# ------------------------------------------------------------------ pieces
abstract = s[s.index("\\begin{abstract}") + 16 : s.index("\\end{abstract}")].strip()
keywords, _ = arg_after(s, "\\keywords", re.search(r"^\\keywords\{", s, re.M).start())
title_orig, _ = arg_after(s, "\\title", re.search(r"^\\title\{", s, re.M).start())

body_start = s.index("\\section{Introduction}")
data_av = s.index("\\section{Data availability}")
ack = s.index("\\section{Acknowledgements}")
bibl = s.index("\\bibliography{")
supp = s.index("\\section{Supplement}")

body = s[body_start:data_av]
data_av_txt = s[data_av:ack].replace("\\section{Data availability}", "").strip()
ack_txt = s[ack:bibl].replace("\\section{Acknowledgements}", "").strip()
supp_txt = s[supp:s.index("\\end{document}")]

# ------------------------------------------------------- figure inventories
fig_re = re.compile(r"\\begin\{figure\*?\}.*?\\end\{figure\*?\}", re.S)


def parse_fig(block):
    cap, _ = arg_after(block, "\\caption")
    lab = re.search(r"\\label\{([^}]*)\}", block).group(1)
    g = re.search(r"\\includegraphics(\[[^\]]*\])?\{([^}]*)\}", block)
    return dict(caption=cap.strip(), label=lab, opt=g.group(1) or "", path=g.group(2))


main_figs = [parse_fig(b) for b in fig_re.findall(body)]
supp_figs = [parse_fig(b) for b in fig_re.findall(supp_txt)]
main_num = {f["label"]: str(i + 1) for i, f in enumerate(main_figs)}
supp_num = {f["label"]: f"S{i + 1}" for i, f in enumerate(supp_figs)}
tab_labels = re.findall(r"\\begin\{table\*?\}.*?\\label\{([^}]*)\}", body, re.S)
tab_num = {l: str(i + 1) for i, l in enumerate(tab_labels)}

# supplemental table legends from the S# Table paragraphs
supp_tabs = []
for m in re.finditer(r"\\label\{S(\d+)_Table\}\s*\n(.*?)(?=\n\s*\n|\n% ---|\Z)", supp_txt, re.S):
    supp_tabs.append((int(m.group(1)), m.group(2).strip()))
supp_tabs.sort()


# ---------------------------------------------------------- text converters
US = [(r"\banalysed\b", "analyzed"), (r"\breanalysed\b", "reanalyzed"), (r"\bcentre\b", "center"),
      (r"\bgrey\b", "gray"), (r"\bGrey\b", "Gray"), (r"\bcoloured\b", "colored"), (r"\bcolour\b", "color"),
      (r"\blabelled\b", "labeled"), (r"\bmodelled\b", "modeled")]


def convert(t, in_supp=False):
    # supplemental figure refs, with panel letters kept
    t = re.sub(r"Supplementary Figs?\.?~?\\ref\{([^}]*)\}",
               lambda m: "Supplemental Figure~" + supp_num[m.group(1)], t)
    t = re.sub(r"\\ref\{((?:suppfig|fig:supp)[^}]*)\}", lambda m: supp_num[m.group(1)], t)
    # supplemental tables
    t = re.sub(r"Supplementary Tables(~| )", r"Supplemental Tables\1", t)
    t = re.sub(r"Supplementary Table(~| )", r"Supplemental Table\1", t)
    t = t.replace("Supplementary Table~S", "Supplemental Table~S")
    # main figures: Fig.~\ref -> Figure~\ref (numbers resolved in the supplement only)
    if in_supp:
        t = re.sub(r"Fig\.~\\ref\{([^}]*)\}", lambda m: "Figure~" + main_num[m.group(1)], t)
        t = re.sub(r"Table~\\ref\{([^}]*)\}", lambda m: "Table~" + tab_num[m.group(1)], t)
    else:
        t = re.sub(r"Fig\.~\\ref\{", r"Figure~\\ref{", t)
    t = t.replace("\\citealp{", "\\citealp{")
    for a, b in US:
        t = re.sub(a, b, t)
    return t


# ------------------------------------------------------------------ body
body = re.sub(r"\\lettrine\[[^\]]*\]\{\\color\{color2\}N\}\{itrogen\}", "Nitrogen", body)
body = body.replace("\\section{Methods}", "\\section{Materials and methods}")
body = re.sub(r"(\\includegraphics(?:\[[^\]]*\])?\{)final/", r"\1../", body)
body = convert(body)

# figure legends get an alt-text placeholder (Plant Physiology requires alt text under each legend)
def add_alt(m):
    blk = m.group(0)
    return blk.replace("\\end{figure", "% Alt text: TODO (required by Plant Physiology for each figure)\n\\end{figure", 1)
body = fig_re.sub(add_alt, body)

data_av_txt = convert(data_av_txt).replace("supplementary tables", "Supplemental Data")

# ------------------------------------------------- acknowledgments / funding
funding = (
    "Part of this work and N.T. were supported by the U.S. Department of Energy, Office of Science, Biological and "
    "Environmental Research program, Early Career Award DE-SC0021889. This work was also supported in part by the "
    "Science and Technologies for Phosphorus Sustainability (STEPS) Center, a National Science Foundation "
    "Science and Technology Center (CBET-2019435). S.F. was supported by the Marcus S. Zuber Assistantship and "
    "the Oak Ridge Institute for Science and Education CERCA Project. S.F.-G. was supported by USDA-ARS."
)
ack_only = (
    "We thank members of the Rub\\'en Rell\\'an-\\'Alvarez and Jim Holland labs for their contribution to fieldwork "
    "and sampling, Drs. Daniel Kick and Ha Ngoc Duong for their guidance in conducting the GWAS, and members of the "
    "Flint-Garcia lab for pollination and harvest of the QTL mapping population. Mention of trade names or commercial "
    "products in this publication is solely for the purpose of providing specific information and does not imply "
    "recommendation or endorsement by the U.S. Department of Agriculture. USDA is an equal opportunity provider and employer."
)

# ----------------------------------------------------- supplemental data list
def short_title(cap):
    c = cap.lstrip()
    if c.startswith("\\textbf{"):
        t, _ = arg_after(c, "\\textbf")
        t = t.strip()
    else:
        t = re.split(r"(?<=\.)\s", cap.strip(), 1)[0]
    t = t.strip()
    return t if t.endswith(".") else t + "."

supp_list = []
for f in supp_figs:
    supp_list.append(f"\\textbf{{Supplemental Figure {supp_num[f['label']]}.}} {convert(short_title(f['caption']), True)}")
for n, txt in supp_tabs:
    supp_list.append(f"\\textbf{{Supplemental Table S{n}.}} {convert(short_title(txt), True)}")
supp_list_tex = "\n\n".join(supp_list)

# --------------------------------------------------------------- abstract
abstract_pp = (
    "Grain nitrogen in maize depends on nitrogen uptake, assimilation, and allocation, but it is unclear whether the same "
    "metabolic functions underlie natural variation in grain amino acid composition, genetic variation associated with soil "
    "nitrogen, and response to selection. We tested this in three populations, each with a distinct nitrogen-related phenotype. "
    "In the Goodman--Buckler Association Panel, grain free amino acids were dominated by asparagine and proline and were highly "
    "heritable, and genome-wide association mapping of 71 traits implicated genes for amino acid metabolism and transport. "
    "The loci shared across the most traits lay near a peptidase A1 gene and endoplasmic reticulum membrane protein genes. "
    "In 3,107 georeferenced landraces from across the Americas, association mapping of soil nitrogen at collection sites "
    "identified variants near an aspartate kinase isozyme and an omega-amidase that passed tests for confounding with "
    "geography and climate. Indian Chief and Jarvis responded to 14 cycles of reciprocal recurrent selection for yield with "
    "different changes in nitrogen allocation. The two populations diverged genetically, and Indian Chief accumulated more "
    "free amino acids, mostly proline, while protein-bound amino acids were unchanged. Their candidate regions included "
    "protein-turnover genes in Indian Chief and amino acid importers in Jarvis, consistent with these population-specific "
    "phenotypes. Each of the four analyses highlighted a different enzyme of the aspartate-to-homoserine pathway, the only "
    "pathway that remained significant after correction across all pathways. Together, these results point to "
    "aspartate-family metabolism as a shared source of candidate genes for how nitrogen is partitioned among grain amino "
    "acid pools in maize."
)
n_abs = len(re.sub(r"[\\{}]", " ", abstract_pp).split())

TITLE = "The aspartate-to-homoserine pathway recurs across maize grain amino acids, soil nitrogen, and yield selection"
SHORT = "Aspartate pathway and maize grain nitrogen"
assert len(TITLE) <= 110, len(TITLE)
assert len(SHORT) <= 50

# --------------------------------------------------------------- manuscript
ms = r"""%% Plant Physiology submission version of main.tex (built by make_plantphys.py).
%% Edit main.tex as the source of truth, or edit this file directly if you prefer;
%% the two are not kept in sync automatically.
%% Compile from final/plant_physiology/:  pdflatex, bibtex, pdflatex, pdflatex
\documentclass[namedate,unnumsec,webpdf,traditional,large]{oup-authoring-template}
\onecolumn
\usepackage{tabularx}
\usepackage{booktabs}
\usepackage[switch]{lineno}
\newcolumntype{b}{X}
\newcolumntype{s}{>{\hsize=.5\hsize}X}
\newcommand{\societylogo}{}

\begin{document}

\journaltitle{Plant Physiology}
\DOI{DOI added during production}
\copyrightyear{2026}
\pubyear{2026}
\access{Advance Access Publication Date: Day Month Year}
\appnotes{Research Article}
\firstpage{1}

%% Title: Plant Physiology allows at most 110 characters (this one is %(tlen)d).
%% Title in main.tex: %(torig)s
\title[%(short)s]{%(title)s}

\author[1,$\dagger$]{Nirwan Tandukar}
\author[2,3,$\dagger$]{Sarah Fitzsimmons}
\author[4,1]{Sherry Flint-Garcia}
\author[2,5]{Ruthie Angelovici}
\author[1,$\ast$]{Rub\'en Rell\'an-\'Alvarez}

\address[1]{\orgdiv{Department of Molecular and Structural Biochemistry, N.C. Plant Sciences Initiative, Genetics and Genomics Program}, \orgname{North Carolina State University}, \orgaddress{\state{Raleigh, NC}, \country{USA}}}
\address[2]{\orgdiv{Division of Biological Sciences}, \orgname{University of Missouri}, \orgaddress{\state{Columbia, MO 65211}, \country{USA}}}
\address[3]{\orgdiv{Department of Agronomy and Plant Genetics}, \orgname{University of Minnesota}, \orgaddress{\state{St. Paul, MN}, \country{USA}}}
\address[4]{\orgname{USDA-ARS Plant Genetics Research Unit}, \orgaddress{\state{Columbia, MO 65211}, \country{USA}}}
\address[5]{\orgdiv{Department of Plant Biology}, \orgname{Michigan State University}, \orgaddress{\state{East Lansing, MI}, \country{USA}}}

\corresp[$\dagger$]{These authors contributed equally to this work. $^{\ast}$Corresponding authors. Sherry Flint-Garcia (\href{mailto:sherry.flint-garcia@usda.gov}{sherry.flint-garcia@usda.gov}), Ruthie Angelovici (\href{mailto:angelovicir@missouri.edu}{angelovicir@missouri.edu}), and Rub\'en Rell\'an-\'Alvarez (\href{mailto:rrellan@ncsu.edu}{rrellan@ncsu.edu})}

%% Plant Physiology: abstract at most 250 words, no citations or abbreviations (this one is %(nabs)d words).
\abstract{%(abstract)s}

\keywords{%(keywords)s}

\maketitle
\linenumbers

%(body)s
\section{Accession numbers}
Sequence data from this article can be found in the NCBI Sequence Read Archive under BioProject PRJNAXXXXXX. Gene identifiers refer to the B73 RefGen\_v5 annotation (Zm00001eb).
%% TODO: fill in the BioProject accession

\section{Supplemental data}
The following materials are available in the online version of this article.

%(supp_list)s

\section{Acknowledgments}
%(ack)s

\section{Author contributions}
%% TODO: Plant Physiology requires one sentence per author (designed research, performed research, contributed tools, analyzed data, wrote the paper).
TODO.

\section{Funding}
%(funding)s

\section{Conflict of interest statement}
None declared.

\section{Data availability}
%(data)s

\bibliographystyle{oup-abbrvnat}
\bibliography{../../sn-bibliography}

\end{document}
""" % dict(tlen=len(TITLE), torig=title_orig.replace("\n", " "), short=SHORT, title=TITLE, nabs=n_abs,
           abstract=abstract_pp, keywords=keywords, body=body.strip() + "\n", supp_list=supp_list_tex,
           ack=ack_only, funding=funding, data=data_av_txt)

# --------------------------------------------------------------- supplement
sf = []
for f in supp_figs:
    path = re.sub(r"^final/", "../", f["path"])
    opt = f["opt"] if f["opt"] else "[width=\\textwidth]"
    sf.append("\\begin{figure}[p]\n\\centering\n\\includegraphics%s{%s}\n\\caption*{\\textbf{Supplemental Figure %s.} %s}\n\\end{figure}\n\\clearpage"
              % (opt, path, supp_num[f["label"]], convert(f["caption"], True)))
st = []
for n, txt in supp_tabs:
    st.append("\\noindent\\textbf{Supplemental Table S%d.} %s\n" % (n, convert(txt, True)))

sup = r"""%% Plant Physiology Supplemental Data (built by make_plantphys.py)
\documentclass[11pt]{article}
\usepackage[margin=2cm]{geometry}
\usepackage{graphicx}
\usepackage{caption}
\usepackage[authoryear,round]{natbib}
\usepackage{amsmath,amssymb}
\usepackage[hidelinks]{hyperref}
\captionsetup{labelformat=empty,font=small}
\begin{document}
\begin{center}
{\Large\textbf{Supplemental Data}}\\[6pt]
%(title)s\\[4pt]
Tandukar, Fitzsimmons, Flint-Garcia, Angelovici, and Rell\'an-\'Alvarez
\end{center}

\section*{Supplemental Figures}
%(figs)s

\section*{Supplemental Tables}
Supplemental Tables S1--S%(ntab)d are provided as separate CSV files.

%(tabs)s

\bibliographystyle{oup-abbrvnat}
\bibliography{../../sn-bibliography}
\end{document}
""" % dict(title=TITLE, figs="\n\n".join(sf), ntab=len(supp_tabs), tabs="\n".join(st))

os.makedirs(out, exist_ok=True)
open(os.path.join(out, "PlantPhys_manuscript.tex"), "w", encoding="utf-8").write(ms)
open(os.path.join(out, "PlantPhys_supplement.tex"), "w", encoding="utf-8").write(sup)
print("main figs", main_num)
print("supp figs", supp_num)
print("supp tables", len(supp_tabs), "abstract words", n_abs, "title chars", len(TITLE))
leftover = re.findall(r"Supplementary [A-Za-z.]+", ms + sup)
print("leftover 'Supplementary':", leftover[:10])
