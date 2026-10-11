suppressPackageStartupMessages(library(data.table))
D <- "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/results/soil_nitrogen_gwas/env_landrace"; M <- 495755; W <- 25000
MODELS <- c("MLM","BLINK","FarmCPU")
ENV <- c("latitude","abs_latitude","longitude","elevation","pH","CEC","aridity")
rd <- function(t,m) fread(file.path(D,sprintf("%s__%s.csv.gz",t,m)),showProgress=FALSE)[,
        .(trait=t, model=m, SNP, Chr, Pos, P=P.value, MAF, fdr_sig)]
traits <- unique(sub("__.*","",list.files(D,pattern="\\.csv\\.gz$")))
all <- rbindlist(lapply(traits,function(t) rbindlist(lapply(MODELS,function(m) rd(t,m)))))
s <- all[trait=="N_5_15" & fdr_sig==TRUE]
cat("N_5_15 per model:", paste(sprintf("%s=%d",MODELS,sapply(MODELS,function(m) s[model==m,.N])),collapse="  "), "\n")
w <- s[, .(Chr=Chr[1], Pos=Pos[1], n_models=uniqueN(model), models=paste(sort(unique(model)),collapse="+"),
           min_P=min(P), MAF=max(MAF,na.rm=TRUE)), by=SNP]
w[!is.finite(MAF), MAF:=NA]
n05 <- unique(all[trait=="N_0_5" & fdr_sig==TRUE]$SNP)
w[, in_N_0_5 := SNP %in% n05]
es <- unique(all[trait %in% ENV & fdr_sig==TRUE, .(SNP,trait)])
w <- merge(w, es[, .(env_hits=paste(sort(trait),collapse=", ")), by=SNP], by="SNP", all.x=TRUE)
w[is.na(env_hits), env_hits:=""]
cat(sprintf("union %d | >=2 models %d | shared with N_0_5 %d | no env trait %d\n",
    nrow(w), w[n_models>=2,.N], w[in_N_0_5==TRUE,.N], w[env_hits=="",.N]))

gff <- fread(cmd="awk -F'\\t' '$3==\"gene\"' /Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/results/annotation/Zm-B73-REFERENCE-NAM-5.0_genes_only.gff3", header=FALSE, sep="\t")
G <- gff[, .(Chr=suppressWarnings(as.integer(sub("^chr","",V1))), Start=V4, End=V5, GeneID=sub("^ID=([^;]+).*","\\1",V9))][!is.na(Chr)]
setkey(G, Chr, Start, End)
q <- w[, .(Chr, s=pmax(1,Pos-W), e=Pos+W, idx=.I)]
ov <- foverlaps(q, G, by.x=c("Chr","s","e"), by.y=c("Chr","Start","End"), nomatch=0L)
a <- cbind(w[ov$idx], ov[, .(GeneID, GeneStart=Start, GeneEnd=End)])
a[, Distance_bp := fifelse(Pos>=GeneStart & Pos<=GeneEnd, 0, fifelse(Pos<GeneStart, GeneStart-Pos, Pos-GeneEnd))]
go <- rbindlist(lapply(c("GO_all_maize_1.txt","GO_all_maize_2.txt"), function(f)
  fread(file.path("/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/results/GO",f), sep="\t", header=FALSE, fill=TRUE)), fill=TRUE)
setnames(go, c("long","GeneID","Desc","PantherFamily","ProtClass","GO_MF","GO_BP","GO_CC"))
go[, GeneName := tstrsplit(Desc,";",fixed=TRUE)[[1]]]
go <- unique(go[, .(GeneID, GeneName, PantherFamily)])[, .SD[1], by=GeneID]
a <- merge(a, go, by="GeneID", all.x=TRUE)
setorder(a, min_P, Distance_bp)
dir.create("out", showWarnings=FALSE)
fwrite(a[, .(SNP, Chr, Pos, n_models, models, min_P, MAF, in_N_0_5, env_hits,
             GeneID, Distance_bp, GeneName, PantherFamily)], "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/results/soil_nitrogen_gwas/N_5_15_FDR_SNP_genes.csv")
writeLines(sort(unique(a$GeneID)), "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/results/soil_nitrogen_gwas/N_5_15_FDR_geneIDs.txt")
core <- a[in_N_0_5==TRUE]
writeLines(sort(unique(core$GeneID)), "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/results/soil_nitrogen_gwas/N_both_depths_core_geneIDs.txt")
cat(sprintf("genes within 25 kb: %d | core (both depths) genes: %d\n", uniqueN(a$GeneID), uniqueN(core$GeneID)))

## nitrogen / amino-acid keyword scan, same as the amino-acid section
PAT <- "nitrate|nitrite|ammoni|NRT|NPF|glutamine synth|glutamate synth|glutamate dehydrog|asparagine synth|asparaginase|aspartate|aspartokinase|homoserine|threonine|methionine|cystathionine|lysine|amino acid|urea|DUR3|AMT|permease|aminotransferase|transaminase|GOGAT|nitrilase|amidase"
hit <- a[grepl(PAT, paste(GeneName, PantherFamily), ignore.case=TRUE)]
cat("\n=== N / amino-acid annotated genes at 5-15 cm FDR 5% ===\n")
print(unique(hit[, .(SNP, mod=n_models, P=signif(min_P,2), both=in_N_0_5, env=env_hits,
       GeneID, kb=round(Distance_bp/1000,1), gene=substr(GeneName,1,48))])[order(P)], nrows=60)
