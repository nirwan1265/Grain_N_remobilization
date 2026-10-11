#!/usr/bin/env Rscript
################################################################################
## 01_GWAS_models.R -- MLM, BLINK and FarmCPU for one trait, on the cluster.
##
##   Rscript 01_GWAS_models.R <trait> [models|collect]
##
## Server version of final/scripts/GWAS_models.R.  The laptop version reads
## per-chromosome HapMap files; those do not exist here, so this one reads the
## numeric matrix written by 00_vcf_to_matrix.py and hands GAPIT GD/GM directly.
## Markers are in VCF order, so each chromosome is a contiguous byte range and
## only that chromosome is ever in memory.
##
## MLM and BLINK run per chromosome (as the published scans did); FarmCPU runs
## once genome-wide through a filebacked big.matrix, which is the only way to
## give it cofactors from every chromosome.
################################################################################

gv <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }

DIR_GENO   <- gv("GN_GENO",   "/share/maize/ntanduk/GrainN/GWAS_env/geno")
PHENO_FILE <- gv("GN_PHENO",  "/share/maize/ntanduk/GrainN/GWAS_env/control_phenotypes_maize_landrace.csv")
PCA_FILE   <- gv("GN_PCA",    "/share/maize/ntanduk/GrainN/GWAS_env/PCA_maize_romerro.csv")
OUTROOT    <- gv("GN_OUTDIR", "/share/maize/ntanduk/GrainN/GWAS_env/out")
NPC        <- as.integer(gv("GN_NPC", "3"))
CHRS       <- as.integer(strsplit(gv("GN_CHRS", paste(1:10, collapse=",")), ",")[[1]])
MATCH      <- gv("GN_MATCH_TAXA", "N_0_5")   # "" disables panel matching

TRAIT_COL <- c(latitude="latitude", lat="latitude",
               abs_latitude="abs_latitude", abslat="abs_latitude",
               longitude="longitude", lon="longitude",
               elevation="elevation", elev="elevation",
               pH="pH", ph="pH", CEC="CEC", cec="CEC",
               aridity="aridity", aridity_index="aridity",
               N_0_5="N_0_5", `N_0-5`="N_0_5", soilN="N_0_5",
               N_5_15="N_5_15", `N_5-15`="N_5_15")
VALID   <- c("MLM", "MLMM", "BLINK", "FarmCPU")
DEFAULT <- c("MLM", "BLINK", "FarmCPU")

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  cat("usage: Rscript 01_GWAS_models.R <trait> [models|collect]\n",
      "traits: ", paste(unique(names(TRAIT_COL)), collapse=", "), "\n",
      "models: ", paste(VALID, collapse=", "), "  default ",
      paste(DEFAULT, collapse=","), "\n", sep="")
  quit(status = 1)
}
TRAIT <- args[1]
if (!TRAIT %in% names(TRAIT_COL)) stop("unknown trait '", TRAIT, "'")
COL <- TRAIT_COL[[TRAIT]]; TAG <- gsub("[^A-Za-z0-9]", "_", TRAIT)
COLLECT_ONLY <- length(args) >= 2 && tolower(args[2]) == "collect"
MODELS <- if (COLLECT_ONLY) c("MLM","MLMM","BLINK") else
          if (length(args) >= 2) trimws(strsplit(args[2], ",")[[1]]) else DEFAULT
bad <- setdiff(MODELS, VALID); if (length(bad)) stop("unknown model(s): ", paste(bad, collapse=", "))

OUTDIR <- file.path(OUTROOT, TAG); dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
sink(file.path(OUTDIR, sprintf("GWAS_models_%s.log", TAG)), split = TRUE)
cat("================================================================\n")
cat(" trait : ", TRAIT, " (column '", COL, "')\n", sep="")
cat(" models: ", paste(MODELS, collapse=", "),
    if (COLLECT_ONLY) "  (collect only)" else "", "\n", sep="")
cat(" geno  : ", DIR_GENO, "\n out   : ", OUTDIR, "\n started: ", format(Sys.time()), "\n", sep="")
cat("================================================================\n\n")

suppressPackageStartupMessages({library(data.table); library(bigmemory); library(biganalytics)})
GAPIT_SOURCE <- Sys.getenv("GN_GAPIT_SOURCE")
if (requireNamespace("GAPIT", quietly = TRUE)) {
  suppressPackageStartupMessages(library(GAPIT))
  cat("GAPIT: installed package", as.character(utils::packageVersion("GAPIT")), "\n")
} else if (nzchar(GAPIT_SOURCE)) {
  for (p in c("multtest","gplots","genetics","EMMREML","scatterplot3d","ape","lme4","MASS","compiler","snowfall"))
    suppressWarnings(suppressPackageStartupMessages(require(p, character.only=TRUE, quietly=TRUE)))
  cat("GAPIT: sourcing", GAPIT_SOURCE, "\n"); source(GAPIT_SOURCE)
  if (!exists("GAPIT")) stop("sourcing ", GAPIT_SOURCE, " did not define GAPIT()")
} else stop("GAPIT not available; install it or set GN_GAPIT_SOURCE")

## ------------------------------------------------------------------ inputs --
need <- file.path(DIR_GENO, c("dims.txt","taxa.txt","GM.csv","GD.int8.bin"))
miss <- need[!file.exists(need)]
if (length(miss)) stop("missing from ", DIR_GENO, ":\n  ",
                       paste(basename(miss), collapse="\n  "),
                       "\nRun 00_vcf_to_matrix.py first.")
dims <- scan(file.path(DIR_GENO,"dims.txt"), quiet=TRUE)
n <- as.integer(dims[1]); m <- as.integer(dims[2])
taxa <- readLines(file.path(DIR_GENO,"taxa.txt"))
GM <- fread(file.path(DIR_GENO,"GM.csv"))
stopifnot(nrow(GM) == m, length(taxa) == n)
cat(sprintf("genotypes: %d taxa x %s markers\n", n, format(m, big.mark=",")))
EXPECTED_SNPS <- as.integer(gv("GN_EXPECTED_SNPS", as.character(m)))

Y <- fread(PHENO_FILE); setnames(Y, 1, "Taxa")
if (!COL %in% names(Y)) stop("column '", COL, "' not in ", PHENO_FILE)
Y <- Y[, c("Taxa", COL), with=FALSE]; setnames(Y, 2, TAG); Y <- Y[!is.na(get(TAG))]
cat(sprintf("phenotyped taxa: %d\n", nrow(Y)))
if (nzchar(MATCH) && COL != MATCH && MATCH %in% names(fread(PHENO_FILE, nrows=0))) {
  keep <- fread(PHENO_FILE, select=c("Taxa", MATCH))[!is.na(get(MATCH)), Taxa]
  b <- nrow(Y); Y <- Y[Taxa %in% keep]
  cat(sprintf("matched to the %s panel: %d -> %d taxa\n", MATCH, b, nrow(Y)))
}
P_all <- fread(PCA_FILE); setnames(P_all, 1, "Taxa")
evs <- paste0("EV", seq_len(NPC))
if (!all(evs %in% names(P_all))) stop("PCA file needs ", paste(evs, collapse=", "))
P_all <- P_all[, c("Taxa", evs), with=FALSE]
P <- P_all[Taxa %in% Y$Taxa]
if (anyNA(P)) stop("PCA file has NA for some phenotyped taxa")
cat(sprintf("covariates: %d PCs for %d taxa\n\n", NPC, nrow(P)))

## --------------------------------------------- read one chromosome of GD ----
read_chr <- function(chr) {
  idx <- which(GM$Chromosome == chr)
  if (!length(idx)) return(NULL)
  if (any(diff(idx) != 1L)) stop("markers for chromosome ", chr, " are not contiguous")
  con <- file(file.path(DIR_GENO,"GD.int8.bin"), "rb")
  on.exit(close(con))
  seek(con, where = as.double(idx[1]-1) * n, origin = "start")
  v <- readBin(con, integer(), n = n*length(idx), size = 1L, signed = TRUE)
  g <- matrix(v, nrow = n)
  colnames(g) <- GM$SNP[idx]
  list(GD = data.frame(taxa = taxa, g, check.names = FALSE, stringsAsFactors = FALSE),
       GM = as.data.frame(GM[idx]))
}

## ------------------------------------------------------- MLM / MLMM / BLINK --
gapit_models <- intersect(MODELS, c("MLM","MLMM","BLINK"))
if (length(gapit_models) && !COLLECT_ONLY) {
  owd <- getwd(); t0 <- Sys.time()
  for (chr in CHRS) {
    wdc <- file.path(OUTDIR, "gapit", paste0("chr", chr))
    dir.create(wdc, recursive=TRUE, showWarnings=FALSE)
    cat(sprintf("  chr%-2d loading ... ", chr)); utils::flush.console()
    gg <- read_chr(chr)
    if (is.null(gg)) { cat("no markers, skipped\n"); next }
    cat(sprintf("%s markers; running %s\n", format(nrow(gg$GM), big.mark=","),
                paste(gapit_models, collapse="+"))); utils::flush.console()
    setwd(wdc)
    ok <- try(GAPIT(Y = as.data.frame(Y), GD = gg$GD, GM = gg$GM,
                    CV = as.data.frame(P), model = gapit_models,
                    SNP.MAF = 0, file.output = TRUE), silent = FALSE)
    setwd(owd)
    if (inherits(ok,"try-error"))
      warning("GAPIT failed on chromosome ", chr, " -- it will be missing from the scan")
    rm(gg, ok); invisible(gc())
  }
  cat(sprintf("GAPIT finished in %.1f min\n\n",
              as.numeric(difftime(Sys.time(), t0, units="mins"))))
}

## ------------------------------------------------------------------ FarmCPU --
if ("FarmCPU" %in% MODELS && !COLLECT_ONLY) {
  .bm <- function(x, ...) { mm <- x[,]; if (!is.matrix(mm)) mm <- matrix(mm, nrow=nrow(x)); mm }
  suppressWarnings(suppressMessages(try(methods::setMethod("as.matrix","big.matrix",.bm), silent=TRUE)))
  assign("as.matrix.big.matrix", .bm, envir = globalenv())
  .chk <- try(as.matrix(as.big.matrix(matrix(1:4,2,2), type="char")), silent=TRUE)
  if (inherits(.chk,"try-error") || !identical(dim(.chk), c(2L,2L)))
    stop("the as.matrix patch did not take -- FarmCPU would fail in FarmCPU.Burger")
  cat("as.matrix patch for big.matrix: verified\n")

  desc <- file.path(DIR_GENO,"GD.desc"); rebuild <- TRUE
  if (file.exists(desc)) {
    GD <- attach.big.matrix(desc)
    if (ncol(GD) == m + 1L) { rebuild <- FALSE; cat("attaching existing big.matrix\n") }
  }
  lock <- file.path(DIR_GENO, ".rebuild.lock")
  if (rebuild) {
    if (file.exists(lock)) stop("another job is building the big.matrix (", lock, ")")
    writeLines(c(format(Sys.time()), Sys.getenv("LSB_JOBID")), lock)
    for (f in c("GD.desc","GD.bin")) { fp <- file.path(DIR_GENO,f); if (file.exists(fp)) file.remove(fp) }
    cat("building the big.matrix (one-off)...\n")
    options(bigmemory.typecast.warning = FALSE)
    GD <- filebacked.big.matrix(nrow=n, ncol=m+1L, type="char",
            backingpath=DIR_GENO, backingfile="GD.bin", descriptorfile="GD.desc")
    GD[,1] <- 0L                                    # placeholder taxa column
    con <- file(file.path(DIR_GENO,"GD.int8.bin"), "rb"); ch <- 20000L
    for (s in seq(1L, m, by=ch)) {
      kk <- min(ch, m-s+1L)
      v <- readBin(con, integer(), n=n*kk, size=1L, signed=TRUE)
      GD[, (s+1L):(s+kk)] <- matrix(v, nrow=n)
      cat(sprintf("\r  %d/%d", s+kk-1L, m)); utils::flush.console()
    }
    close(con); cat("\n"); unlink(lock)
  }
  stopifnot(ncol(GD) == m + 1L, nrow(GD) == n)

  ord <- data.table(Taxa = taxa)
  Yf <- Y[ord, on="Taxa"]; Pf <- P_all[ord, on="Taxa"]
  if (anyNA(Pf[,-1])) stop("PCA file lacks PCs for ", sum(!complete.cases(Pf[,-1])), " genotyped taxa")
  cat(sprintf("FarmCPU: %d of %d taxa phenotyped\n", sum(!is.na(Yf[[TAG]])), n))
  FC <- tryCatch(get("FarmCPU", envir=asNamespace("GAPIT")),
        error=function(e) tryCatch(get("FarmCPU"), error=function(e2) stop("FarmCPU not found")))
  fwd <- file.path(OUTDIR,"farmcpu"); dir.create(fwd, recursive=TRUE, showWarnings=FALSE)
  owd2 <- getwd(); setwd(fwd); t0 <- Sys.time()
  myFarm <- FC(Y=as.data.frame(Yf), GD=GD, GM=as.data.frame(GM),
               CV=as.matrix(Pf[,-1,with=FALSE]), method.bin="optimum",
               bin.size=c(5e5,5e6,5e7), bin.selection=seq(10,100,10),
               maxLoop=10, QTN.threshold=0.01, MAF.calculate=FALSE, file.output=FALSE)
  setwd(owd2)
  cat(sprintf("FarmCPU finished in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
  saveRDS(myFarm, file.path(fwd, sprintf("FarmCPU_%s_GENOMEWIDE.rds", TAG)))
  res <- as.data.table(myFarm$GWAS)
  setnames(res, c("Chromosome","Position","p.value"), c("Chr","Pos","P.value"), skip_absent=TRUE)
  if (!"MAF" %in% names(res)) res[, MAF := NA_real_]
  if (!"nobs" %in% names(res)) res[, nobs := sum(!is.na(Yf[[TAG]]))]
  if ("effect" %in% names(res)) setnames(res, "effect", "Effect")
  if (!"Effect" %in% names(res)) res[, Effect := NA_real_]
  res[, `H&B.P.Value` := p.adjust(P.value, method="BH")]
  fwrite(res[, .(SNP,Chr,Pos,P.value,MAF,nobs,`H&B.P.Value`,Effect)],
         file.path(OUTDIR, sprintf("raw_GWAS_FarmCPU_3PC_%s.txt", TAG)), sep="\t")
  sq <- myFarm$seqQTN
  if (length(sq)) {
    qt <- cbind(data.table(rank=seq_along(sq)), GM[sq]); qt$P.value <- res$P.value[sq]
    fwrite(qt, file.path(OUTDIR, sprintf("FarmCPU_pseudoQTN_%s.csv", TAG)))
    cat("pseudo-QTNs:", length(sq), "across", uniqueN(qt$Chromosome), "chromosomes\n")
  }
}

## --------------------------------------------------------------- collection --
## GAPIT can write the same scan twice per chromosome under "(NYC)" / "(Kansas)"
## suffixes -- identical p-values, differing only in whether the effect column
## is "Effect" or "effect".  Take exactly one file per chromosome.
if (length(gapit_models)) {
  wd <- file.path(OUTDIR, "gapit")
  if (dir.exists(wd)) {
    cat("\ncollecting GAPIT output from ", wd, "\n", sep="")
    for (mod in gapit_models) {
      hits <- list.files(wd, recursive=TRUE, full.names=TRUE,
                         pattern=sprintf("GWAS_Results\\.%s\\..*\\.csv$", mod))
      if (!length(hits)) { if (!COLLECT_ONLY) warning("no result files for ", mod); next }
      bydir <- split(hits, dirname(hits))
      chosen <- vapply(bydir, function(fs) {
        if (length(fs)==1L) return(fs)
        sc <- vapply(fs, function(f)
          as.numeric("Effect" %in% names(fread(f, nrows=0))), numeric(1))
        fs[which.max(sc)] }, character(1))
      cat(sprintf("  %-8s %d chromosome dir(s)%s\n", mod, length(chosen),
          if (length(hits)>length(chosen)) sprintf(", %d duplicate(s) ignored", length(hits)-length(chosen)) else ""))
      r <- rbindlist(lapply(chosen, function(f) {
             x <- fread(f); setnames(x, names(x), sub("^effect$","Effect", names(x)))
             setnames(x, names(x), sub("^Chromosome$","Chr", names(x)))
             setnames(x, names(x), sub("^Position *$","Pos", names(x))); x }), fill=TRUE)
      if (sum(duplicated(r$SNP))) { warning(mod, ": duplicated SNPs across folders"); r <- unique(r, by="SNP") }
      for (nd in c("MAF","nobs","Effect")) if (!nd %in% names(r)) r[[nd]] <- NA
      if (!"H&B.P.Value" %in% names(r)) r[["H&B.P.Value"]] <- p.adjust(r$P.value, method="BH")
      setorder(r, Chr, Pos)
      outf <- file.path(OUTDIR, sprintf("raw_GWAS_%s_3PC_%s.txt", mod, TAG))
      fwrite(r[, c("SNP","Chr","Pos","P.value","MAF","nobs","H&B.P.Value","Effect"), with=FALSE],
             outf, sep="\t")
      msg <- if (nrow(r)==EXPECTED_SNPS) "OK" else
             sprintf("*** %d markers, expected %d -- %d MISSING ***",
                     nrow(r), EXPECTED_SNPS, EXPECTED_SNPS-nrow(r))
      cat(sprintf("           -> %s  %s\n", basename(outf), msg))
      if (nrow(r)!=EXPECTED_SNPS) print(r[, .N, by=Chr][order(Chr)])
      PB <- 0.05/EXPECTED_SNPS; PS <- 1/EXPECTED_SNPS
      cat(sprintf("              Bonferroni %.3e: %d | suggestive %.3e: %d\n",
                  PB, sum(r$P.value<=PB, na.rm=TRUE), PS, sum(r$P.value<=PS, na.rm=TRUE)))
    }
  }
}
cat("\nfinished: ", format(Sys.time()), "\n", sep="")
sink()
