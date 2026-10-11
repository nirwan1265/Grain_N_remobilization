suppressPackageStartupMessages(library(data.table))
D   <- "/home/claude/env/landrace"
GFF <- "/home/claude/rerun/genes_only.gff3"
GO  <- "/home/claude/rerun/shiny/amino_gwas_app/data/supplementary/SuppTable1_GWAS_annotation_soilN_amino_with_GO.csv"
OUT <- "/home/claude/env/out"; dir.create(OUT, showWarnings=FALSE)
M <- 495755; WINDOW <- 25000
MODELS <- c("MLM","BLINK","FarmCPU")
ENV <- c("latitude","abs_latitude","longitude","elevation","pH","CEC","aridity")
GEOG <- c("latitude","abs_latitude","longitude","elevation")

rdx <- function(tr,md){d<-fread(file.path(D,sprintf("%s__%s.csv.gz",tr,md)),showProgress=FALSE)
  d[,.(trait=tr, model=md, SNP, Chr, Pos, P.value, MAF, FDR, fdr_sig, bonferroni)]}
traits <- unique(sub("__.*","",list.files(D,pattern="\\.csv\\.gz$")))
all <- rbindlist(lapply(traits,function(t) rbindlist(lapply(MODELS,function(m) rdx(t,m)))))

## ---- soil N: per-model and cross-model counts at FDR -------------------------
sn <- all[trait=="N_0_5"]
w <- dcast(sn[fdr_sig==TRUE], SNP+Chr+Pos ~ model, value.var="P.value", fun.aggregate=min, fill=NA)
for (m in MODELS) if (!m %in% names(w)) w[, (m):=NA_real_]
setnames(w, MODELS, paste0("P_",MODELS))
pc <- paste0("P_",MODELS)
w[, n_models := rowSums(!is.na(.SD)), .SDcols=pc]
w[, min_P    := do.call(pmin, c(.SD, na.rm=TRUE)), .SDcols=pc]
maf <- unique(sn[!is.na(MAF), .(SNP, MAF)], by="SNP"); w <- maf[w, on="SNP"]

cat("================ SOIL N (N_0_5) at FDR 5% ================\n")
cat(sprintf("per model: %s\n", paste(sprintf("%s=%d", MODELS,
    sapply(MODELS, function(m) sn[model==m & fdr_sig==TRUE, .N])), collapse="  ")))
cat(sprintf("union %d | >=2 models %d | all 3 models %d\n",
    nrow(w), w[n_models>=2,.N], w[n_models==3,.N]))
cat("\ncross-model breakdown:\n")
print(w[, .N, by=.(n_models)][order(n_models)])

## ---- env significance per SNP ------------------------------------------------
envsig <- all[trait %in% ENV & fdr_sig==TRUE, .(SNP, trait)]
envtab <- dcast(unique(envsig), SNP ~ trait, fun.aggregate=length, value.var="trait")
for (e in ENV) if (!e %in% names(envtab)) envtab[, (e):=0L]
w <- merge(w, envtab, by="SNP", all.x=TRUE)
for (e in ENV) w[is.na(get(e)), (e):=0L]
w[, n_env   := rowSums(.SD), .SDcols=ENV]
w[, n_geog  := rowSums(.SD), .SDcols=GEOG]
w[, env_hits := apply(.SD, 1, function(r) paste(ENV[r>0], collapse=", ")), .SDcols=ENV]
## N_5_15 replication
n515 <- unique(all[trait=="N_5_15" & fdr_sig==TRUE]$SNP)
w[, replicates_N_5_15 := SNP %in% n515]

cat(sprintf("\nof %d soil-N SNPs: %d hit no environmental trait (%.1f%%), %d hit >=1 geography trait\n",
    nrow(w), w[n_env==0,.N], 100*w[n_env==0,.N]/nrow(w), w[n_geog>0,.N]))
cat(sprintf("replicating at N_5_15: %d (%.1f%%)\n", w[replicates_N_5_15==TRUE,.N],
    100*w[replicates_N_5_15==TRUE,.N]/nrow(w)))

## ---- genes within 25 kb (positions are B73 v5: SNP name carries the v5 coord) -
gff <- fread(cmd=sprintf("awk -F'\\t' '$3==\"gene\"' %s", shQuote(GFF)), header=FALSE, sep="\t")
genes <- gff[, .(Chr=sub("^[Cc]hr","",V1), Start=V4, End=V5,
                 GeneID=sub("^ID=([^;]+).*","\\1",V9))]
genes <- genes[Chr %in% as.character(1:10)][, Chr := as.integer(Chr)]
setkey(genes, Chr, Start, End)
q <- w[, .(Chr, s=pmax(1,Pos-WINDOW), e=Pos+WINDOW, idx=.I)]
ov <- foverlaps(q, genes, by.x=c("Chr","s","e"), by.y=c("Chr","Start","End"), nomatch=0L)
ann <- cbind(w[ov$idx], ov[, .(GeneID, GeneStart=Start, GeneEnd=End)])
ann[, Relation := fifelse(Pos>=GeneStart & Pos<=GeneEnd, "within",
                   fifelse(Pos<GeneStart, "upstream", "downstream"))]
ann[, Distance_bp := fifelse(Relation=="within", 0,
      fifelse(Relation=="upstream", GeneStart-Pos, Pos-GeneEnd))]
## genome-wide PANTHER export -- the amino-acid annotation table covers only the
## GB candidate genes and leaves most soil-N genes unannotated
go <- rbindlist(lapply(c("GO_all_maize_1.txt","GO_all_maize_2.txt"),
   function(f) fread(file.path("/home/claude/rerun/final/results/GO", f),
                     sep="\t", header=FALSE, fill=TRUE)), fill=TRUE)
setnames(go, c("long","GeneID","Desc","PantherFamily","ProtClass","GO_MF","GO_BP","GO_CC"))
go[, GeneName := tstrsplit(Desc, ";", fixed=TRUE)[[1]]]
go <- unique(go[, .(GeneID, GeneName, PantherFamily, ProtClass, GO_BP, GO_MF)])[, .SD[1], by=GeneID]
ann <- merge(ann, go, by="GeneID", all.x=TRUE, sort=FALSE)
setorder(ann, min_P)
fwrite(ann, file.path(OUT,"soilN_FDR_all_SNPs_with_env.csv"))
fwrite(w[order(min_P)], file.path(OUT,"soilN_FDR_SNP_env_matrix.csv"))

cat(sprintf("\nSNPs with >=1 gene within 25 kb: %d of %d | candidate genes %d\n",
    uniqueN(ann$SNP), nrow(w), uniqueN(ann$GeneID)))

## ---- the specific set: soil-N, no environmental trait ------------------------
spec <- w[n_env==0]
cat(sprintf("\n=== SOIL-N-SPECIFIC SET: %d SNPs ===\n", nrow(spec)))
cat(sprintf("  of these >=2 models: %d | replicating at N_5_15: %d | both: %d\n",
    spec[n_models>=2,.N], spec[replicates_N_5_15==TRUE,.N],
    spec[n_models>=2 & replicates_N_5_15==TRUE,.N]))
cat("\ntop 15 by P, with genes:\n")
top <- merge(spec[order(min_P)][1:15], unique(ann[,.(SNP, GeneID, GeneName, PantherFamily, Relation, Distance_bp)]),
             by="SNP", all.x=TRUE, allow.cartesian=TRUE)
print(top[order(min_P), .(SNP, n_models, minP=signif(min_P,3), MAF=round(MAF,3),
          rep=replicates_N_5_15, GeneID, Gene=substr(GeneName,1,46))], nrows=40)

cat("\n=== the two focal loci from the current manuscript ===\n")
for (s in c("S1_16128632","S5_22338979")) {
  r <- w[SNP==s]
  if (!nrow(r)) { cat(sprintf("  %-14s NOT FDR-significant in this run\n", s)); next }
  g <- unique(ann[SNP==s, .(GeneID, Relation, Distance_bp, GeneName, PantherFamily)])
  cat(sprintf("  %-14s n_models=%d  minP=%.2e  MAF=%.3f  env_hits=[%s]  N_5_15=%s\n",
      s, r$n_models, r$min_P, r$MAF, r$env_hits, r$replicates_N_5_15))
  for (i in seq_len(nrow(g))) cat(sprintf("      %s  %s %dbp  %s\n", g$GeneID[i], g$Relation[i],
      g$Distance_bp[i], substr(paste(g$GeneName[i],"|",g$PantherFamily[i]),1,78)))
}
