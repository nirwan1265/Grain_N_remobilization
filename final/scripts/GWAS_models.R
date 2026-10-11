#!/usr/bin/env Rscript
################################################################################
## GWAS_models.R -- run MLM, MLMM, BLINK and FarmCPU on one trait, genome-wide.
##
##   Rscript GWAS_models.R <trait> [model1,model2,...]
##
##   trait  : latitude | abs_latitude | longitude | elevation |
##            pH | CEC | aridity | N_0_5 | N_5-15
##   models : optional, comma separated.  Default MLM,BLINK,FarmCPU -- the three
##            the consensus actually uses, so the control is procedure-matched
##            to the soil-N analysis.  MLMM is valid but opt-in: it is excluded
##            from the consensus because its soil-N scan is missing all of
##            chromosome 2, and it is the slowest of the four.
##
##   Rscript GWAS_models.R latitude                  # MLM, BLINK, FarmCPU
##   Rscript GWAS_models.R N_5-15 FarmCPU            # just one
##   Rscript GWAS_models.R latitude MLM,MLMM,BLINK,FarmCPU   # all four
##
## Built from the existing pipeline:
##   * covariates, PCA file and hapmap layout are those of
##     Grain_Nitrogen/src/2_run_GWAS_modesl.R
##   * the FarmCPU block is the one from
##     ~/Desktop/farmcpu_genomewide/02_farmcpu_genomewide.R, including all four
##     GAPIT workarounds, and it REUSES the big.matrix already built in
##     ~/Desktop/farmcpu_gw -- nothing is rebuilt.
##
## Two deliberate departures from 2_run_GWAS_modesl.R, both bugs there:
##   1. that script does `out_dir <- paste0(getwd(), ...)` then `setwd(out_dir)`
##      inside the chromosome loop, so the output directories nest
##      (chr1/FarmCPU/N/chr2/FarmCPU/N/chr3/...).  Here every path is absolute
##      and the working directory is restored on exit.
##   2. it runs the loop per chromosome.  MLM is unaffected, but MLMM, BLINK and
##      FarmCPU are multi-locus: run one chromosome at a time they can only draw
##      cofactors from that chromosome.  GAPIT's own file.G interface reads all
##      ten chromosome files and analyses them jointly, so that is used here.
##
## Outputs, one per model, in the layout the downstream scripts already read:
##   <OUTDIR>/raw_GWAS_<MODEL>_3PC_<trait>.txt
##     SNP  Chr  Pos  P.value  MAF  nobs  H&B.P.Value  Effect
################################################################################

## ---------------------------------------------------------------- settings --
gv <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }

PHENO_FILE <- gv("GN_PHENO", path.expand(
  "~/Documents/Research/data/Phenotypes/control_phenotypes_maize_landrace.csv"))
PCA_FILE   <- gv("GN_PCA", path.expand(
  "~/Documents/Github/Grain_Nitrogen/results/PCA/PCA_maize_romerro.csv"))
HAPMAP_DIR <- gv("GN_HAPMAP", path.expand(
  "~/Documents/Research/data/maize_genotype/RomeroNavarro/hapmap/"))
BIGMAT_DIR <- gv("GN_BIGMAT", path.expand("~/Desktop/farmcpu_gw"))
OUTDIR     <- gv("GN_OUTDIR", path.expand("~/Desktop/control_gwas"))

NPC           <- as.integer(gv("GN_NPC", "3"))     # 3 PCs, as in the published scans
EXPECTED_SNPS <- as.integer(gv("GN_EXPECTED_SNPS", "495755"))
HAPMAP_PREFIX <- gv("GN_HAPMAP_PREFIX", "chr")     # files are chr1.txt .. chr10.txt
HAPMAP_EXT    <- gv("GN_HAPMAP_EXT", "txt")
CHR_FROM      <- as.integer(gv("GN_CHR_FROM", "1"))
CHR_TO        <- as.integer(gv("GN_CHR_TO", "10"))

## trait name -> column in PHENO_FILE.  Aliases so the names in the write-up work.
TRAIT_COL <- c(
  latitude = "latitude", lat = "latitude",
  abs_latitude = "abs_latitude", abslat = "abs_latitude",
  longitude = "longitude", lon = "longitude",
  elevation = "elevation", elev = "elevation",
  pH = "pH", ph = "pH",
  CEC = "CEC", cec = "CEC",
  aridity = "aridity", aridity_index = "aridity",
  N_0_5 = "N_0_5", `N_0-5` = "N_0_5", soilN = "N_0_5",
  N_5_15 = "N_5_15", `N_5-15` = "N_5_15")

VALID_MODELS   <- c("MLM", "MLMM", "BLINK", "FarmCPU")
## the three used by soilN_crossmodel_consensus.R
DEFAULT_MODELS <- c("MLM", "BLINK", "FarmCPU")

## ------------------------------------------------------------------- args ---
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  cat("usage: Rscript GWAS_models.R <trait> [models]\n\n",
      "traits : ", paste(unique(names(TRAIT_COL)), collapse = ", "), "\n",
      "models : ", paste(VALID_MODELS, collapse = ", "), "  (comma separated)\n",
      "         default ", paste(DEFAULT_MODELS, collapse = ","),
      " -- the three the consensus uses.\n",
      "         MLMM is opt-in: excluded from the consensus, and the slowest.\n",
      sep = "")
  quit(status = 1)
}
TRAIT <- args[1]
if (!TRAIT %in% names(TRAIT_COL)) stop("unknown trait '", TRAIT, "'. Known: ",
                                       paste(unique(names(TRAIT_COL)), collapse = ", "))
COL    <- TRAIT_COL[[TRAIT]]
## "collect" re-does only the last step -- gathering GAPIT's per-chromosome CSVs
## into the raw_GWAS_* tables.  Use it after a collection bug; it never re-runs
## a scan, so hours of GAPIT are not thrown away.
COLLECT_ONLY <- length(args) >= 2 && tolower(args[2]) == "collect"
MODELS <- if (COLLECT_ONLY) c("MLM", "MLMM", "BLINK") else
          if (length(args) >= 2) trimws(strsplit(args[2], ",")[[1]]) else DEFAULT_MODELS
bad <- setdiff(MODELS, VALID_MODELS)
if (length(bad)) stop("unknown model(s): ", paste(bad, collapse = ", "),
                      "\nvalid: ", paste(VALID_MODELS, collapse = ", "))

TAG <- gsub("[^A-Za-z0-9]", "_", TRAIT)
OUTDIR <- file.path(OUTDIR, TAG)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUTDIR, sprintf("GWAS_models_%s.log", TAG))
sink(LOG, split = TRUE)

cat("================================================================\n")
cat(" trait  : ", TRAIT, "  (column '", COL, "')\n", sep = "")
cat(" models : ", paste(MODELS, collapse = ", "),
    if (COLLECT_ONLY) "   (collect only -- no scan is run)" else
    if (length(args) < 2) "   (default)" else "", "\n", sep = "")
if ("MLMM" %in% MODELS)
  cat("          note: MLMM is not used by the consensus -- running it is optional.\n")
cat(" out    : ", OUTDIR, "\n", sep = "")
cat(" started: ", format(Sys.time()), "\n", sep = "")
cat(" pid    : ", Sys.getpid(), "\n", sep = "")
cat("================================================================\n\n")
## This script is single-threaded.  The only multi-core work is whatever the
## BLAS underneath R does for matrix algebra -- print it so the usage is known.
si <- sessionInfo()
cat("BLAS  : ", if (!is.null(si$BLAS)) basename(si$BLAS) else "unknown", "\n", sep = "")
cat("LAPACK: ", if (!is.null(si$LAPACK)) basename(si$LAPACK) else "unknown", "\n", sep = "")
cat("To use more of the machine, run several TRAITS at once rather than\n",
    "looking for a thread count here -- the traits are independent.\n",
    "Memory, not cores, is the limit: each run holds one chromosome's hapmap.\n\n", sep = "")

suppressPackageStartupMessages({
  library(data.table); library(bigmemory); library(biganalytics)
})
## GAPIT may be an installed package, or just its function file sourced into the
## global environment -- which is how it is usually used on a cluster where
## remotes/devtools are unavailable.  Both work here.
##   export GN_GAPIT_SOURCE=~/gapit_functions.txt
GAPIT_SOURCE <- Sys.getenv("GN_GAPIT_SOURCE")
if (requireNamespace("GAPIT", quietly = TRUE)) {
  suppressPackageStartupMessages(library(GAPIT))
  cat("GAPIT: installed package, version ",
      as.character(utils::packageVersion("GAPIT")), "\n", sep = "")
} else if (nzchar(GAPIT_SOURCE)) {
  ## sourcing the function file does not attach GAPIT's dependencies, so do it
  for (p in c("multtest", "gplots", "genetics", "EMMREML", "scatterplot3d",
              "ape", "lme4", "MASS", "compiler", "snowfall")) {
    if (!suppressWarnings(suppressPackageStartupMessages(
          require(p, character.only = TRUE, quietly = TRUE))))
      warning("GAPIT dependency not available: ", p)
  }
  cat("GAPIT: sourcing functions from ", GAPIT_SOURCE, "\n", sep = "")
  source(GAPIT_SOURCE)
  if (!exists("GAPIT")) stop("sourcing ", GAPIT_SOURCE, " did not define GAPIT()")
} else {
  stop("GAPIT is not available. Either\n",
       "  * install it   -- see install_requirements.R, which does not need\n",
       "                    remotes or devtools if you give it a tarball:\n",
       "                      GN_GAPIT_TARBALL=~/GAPIT.tar.gz Rscript install_requirements.R\n",
       "  * or source it -- put gapit_functions.txt on the server and set\n",
       "                      export GN_GAPIT_SOURCE=~/gapit_functions.txt")
}

## -------------------------------------------------- phenotype + covariates --
for (f in c(PHENO_FILE, PCA_FILE)) if (!file.exists(f)) stop("not found: ", f)

Y <- fread(PHENO_FILE); setnames(Y, 1, "Taxa")
if (!COL %in% names(Y)) stop("column '", COL, "' not in ", PHENO_FILE,
                             "\navailable: ", paste(names(Y), collapse = ", "))
Y <- Y[, c("Taxa", COL), with = FALSE]
setnames(Y, 2, TAG)                       # GAPIT names its output after this
Y <- Y[!is.na(get(TAG))]
cat(sprintf("phenotyped taxa: %d\n", nrow(Y)))

## Match the analysis set to the trait being controlled for.
##
## Soil N is scored for 3,107 accessions; latitude and longitude for all 3,539.
## Left alone, a control scan runs on a larger panel than soil N did, so it has
## more power -- and "this SNP is also a latitude hit" then partly means "the
## latitude scan simply had more accessions".  Restricting every trait to the
## taxa that have soil N makes the panel, the kinship and the PCs identical
## across scans, so the only thing that differs is the phenotype.
##   GN_MATCH_TAXA=N_0_5   (recommended)   GN_MATCH_TAXA=""  to disable
MATCH <- gv("GN_MATCH_TAXA", "N_0_5")
if (nzchar(MATCH) && COL != MATCH) {
  hdr <- names(fread(PHENO_FILE, nrows = 0))
  if (!MATCH %in% hdr) {
    warning("GN_MATCH_TAXA='", MATCH, "' is not a column in ", PHENO_FILE,
            " -- running on all phenotyped taxa")
  } else {
    ref <- fread(PHENO_FILE, select = c("Taxa", MATCH))
    keep <- ref[!is.na(get(MATCH)), Taxa]
    before <- nrow(Y)
    Y <- Y[Taxa %in% keep]
    cat(sprintf("matched to the %s panel: %d -> %d taxa\n", MATCH, before, nrow(Y)))
  }
} else if (!nzchar(MATCH)) {
  cat("GN_MATCH_TAXA is empty -- running on every taxon with this phenotype\n")
}
if (nrow(Y) < 100) stop("fewer than 100 phenotyped taxa -- check the trait column")

P_all <- fread(PCA_FILE); setnames(P_all, 1, "Taxa")
evs <- paste0("EV", seq_len(NPC))
if (!all(evs %in% names(P_all))) stop("the PCA file needs ", paste(evs, collapse = ", "))
P_all <- P_all[, c("Taxa", evs), with = FALSE]
## GAPIT gets the covariates for the phenotyped taxa; FarmCPU needs them for
## every genotyped taxon, phenotyped or not, so it uses P_all.
P <- P_all[Taxa %in% Y$Taxa]
if (anyNA(P)) stop("the PCA file has NA for some phenotyped taxa")
cat(sprintf("covariates: %d PCs for %d phenotyped taxa (%d in the PCA file)\n\n",
            NPC, nrow(P), nrow(P_all)))

## How far is this trait from the trait it is a control for?  Skipped when the
## trait IS soil N -- and the reference column is renamed before the merge, so a
## trait whose own name matches the reference cannot collide with it.
if (COL != "N_0_5") {
  ph_hdr <- names(fread(PHENO_FILE, nrows = 0))
  if ("N_0_5" %in% ph_hdr) {
    ref <- fread(PHENO_FILE, select = c("Taxa", "N_0_5"))
    setnames(ref, "N_0_5", ".refN")
    m  <- merge(Y, ref, by = "Taxa")
    ok <- stats::complete.cases(m[[TAG]], m$.refN)
    if (sum(ok) > 50)
      cat(sprintf("correlation with soil N 0-5 cm: r = %+.3f (n = %d)\n\n",
                  cor(m[[TAG]][ok], m$.refN[ok]), sum(ok)))
  }
} else {
  cat("this run IS the soil N 0-5 cm scan -- use it to check the pipeline\n",
      "reproduces the published counts (MLM 1, BLINK 38, FarmCPU 33 at\n",
      "Bonferroni 1.01e-07) before trusting the control traits.\n\n", sep = "")
}

## ------------------------------------------------------------- GAPIT models --
## MLM, MLMM and BLINK go through GAPIT with a hapmap genotype.
##
## GAPIT3 dropped the multi-file genotype interface (file.G / file.from / file.to),
## so there is no way to hand it all ten chromosomes in one call without holding
## the whole 495,755-marker hapmap in memory at once.  The script checks which
## interface the installed GAPIT actually has and uses it if present; otherwise it
## falls back to the per-chromosome loop of 2_run_GWAS_modesl.R, which is the path
## that produced the published scans.
##
## Per chromosome is exact for MLM (single-locus).  For MLMM and BLINK it is a
## compromise -- they can only draw cofactors from the chromosome in hand -- but
## it is what the soil-N scans did, and the point of this run is to be
## procedure-matched to them.  FarmCPU below is genome-wide either way.
gapit_models <- intersect(MODELS, c("MLM", "MLMM", "BLINK"))
if (length(gapit_models) && !COLLECT_ONLY) {
  if (!dir.exists(HAPMAP_DIR)) stop("hapmap directory not found: ", HAPMAP_DIR)
  hm <- file.path(HAPMAP_DIR, sprintf("%s%d.%s", HAPMAP_PREFIX, CHR_FROM:CHR_TO, HAPMAP_EXT))
  names(hm) <- CHR_FROM:CHR_TO
  miss <- hm[!file.exists(hm)]
  if (length(miss)) stop("missing hapmap files:\n  ", paste(miss, collapse = "\n  "))

  gf <- names(formals(GAPIT))
  has_multifile <- all(c("file.G", "file.from", "file.to") %in% gf)
  cat("installed GAPIT accepts a multi-file genotype interface: ",
      if (has_multifile) "yes" else "no (GAPIT3) -- using the per-chromosome loop",
      "\n", sep = "")
  cat("GAPIT models: ", paste(gapit_models, collapse = ", "), "\n", sep = "")

  read_hapmap <- function(f) {
    ## as in 2_run_GWAS_modesl.R: the hapmap header row is read as data
    read.delim(f, header = FALSE, sep = "\t", quote = "", comment.char = "",
               check.names = FALSE, stringsAsFactors = FALSE)
  }

  owd <- getwd()
  t0 <- Sys.time()
  if (has_multifile) {
    wd <- file.path(OUTDIR, "gapit"); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
    setwd(wd)
    invisible(GAPIT(Y = as.data.frame(Y), CV = as.data.frame(P),
                    file.G = HAPMAP_PREFIX, file.Ext.G = HAPMAP_EXT,
                    file.path = HAPMAP_DIR, file.from = CHR_FROM, file.to = CHR_TO,
                    file.fragment = 512, model = gapit_models,
                    SNP.MAF = 0, file.output = TRUE))
    setwd(owd)
  } else {
    for (chr in CHR_FROM:CHR_TO) {
      wdc <- file.path(OUTDIR, "gapit", paste0("chr", chr))
      dir.create(wdc, recursive = TRUE, showWarnings = FALSE)
      cat(sprintf("  chr%-2d reading hapmap ... ", chr)); utils::flush.console()
      G <- read_hapmap(hm[[as.character(chr)]])
      cat(sprintf("%d markers; running %s\n", nrow(G) - 1L,
                  paste(gapit_models, collapse = "+"))); utils::flush.console()
      setwd(wdc)
      ok <- try(GAPIT(Y = as.data.frame(Y), G = G, CV = as.data.frame(P),
                      model = gapit_models, SNP.MAF = 0, file.output = TRUE),
                silent = FALSE)
      setwd(owd)
      if (inherits(ok, "try-error"))
        warning("GAPIT failed on chromosome ", chr,
                " -- that chromosome will be missing from the concatenated scan")
      rm(G, ok); invisible(gc())
    }
  }
  cat(sprintf("GAPIT finished in %.1f min\n\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## ---------------------------------------------------------------- FarmCPU ---
## Called directly, not through GAPIT(), on the big.matrix in BIGMAT_DIR.
## Everything below is the block from 02_farmcpu_genomewide.R; the comments
## there explain each of the four workarounds and are kept.
if ("FarmCPU" %in% MODELS && !COLLECT_ONLY) {
  need <- file.path(BIGMAT_DIR, c("dims.txt", "taxa.txt", "GM.csv"))
  miss <- need[!file.exists(need)]
  if (length(miss))
    stop("missing from ", BIGMAT_DIR, ":\n  ", paste(basename(miss), collapse = "\n  "),
         "\nRun 01_vcf_to_bigmatrix.py first (see ~/Desktop/farmcpu_genomewide/).")

  ## patch: as.matrix() has no method for big.matrix, and FarmCPU.Burger /
  ## FarmCPU.Remove call it at five sites.
  .bm_as_matrix <- function(x, ...) {
    mm <- x[, ]; if (!is.matrix(mm)) mm <- matrix(mm, nrow = nrow(x)); mm
  }
  suppressWarnings(suppressMessages(
    try(methods::setMethod("as.matrix", "big.matrix", .bm_as_matrix), silent = TRUE)))
  assign("as.matrix.big.matrix", .bm_as_matrix, envir = globalenv())
  .tb <- as.big.matrix(matrix(c(1L, 2L, 3L, 4L), 2, 2), type = "char")
  .chk <- try(as.matrix(.tb), silent = TRUE)
  if (inherits(.chk, "try-error") || !identical(dim(.chk), c(2L, 2L)))
    stop("the as.matrix patch did not take -- FarmCPU would fail in FarmCPU.Burger.")
  rm(.tb, .chk); cat("as.matrix() patch for big.matrix: verified\n")

  dims <- scan(file.path(BIGMAT_DIR, "dims.txt"), quiet = TRUE)
  n <- as.integer(dims[1]); m <- as.integer(dims[2])
  taxa <- readLines(file.path(BIGMAT_DIR, "taxa.txt"))
  GM   <- fread(file.path(BIGMAT_DIR, "GM.csv"))
  stopifnot(nrow(GM) == m, length(taxa) == n)
  cat(sprintf("genotypes: %d taxa x %d SNPs\n", n, m))

  ## the placeholder first column: FarmCPU's orientation test fires when
  ## ncol(GD) == nm, which is exactly the taxa-free layout its header asks for.
  desc <- file.path(BIGMAT_DIR, "GD.desc")
  rebuild <- TRUE
  if (file.exists(desc)) {
    GD <- attach.big.matrix(desc)
    if (ncol(GD) == m + 1L) { rebuild <- FALSE; cat("attaching the existing big.matrix\n") }
    else cat("existing big.matrix has", ncol(GD), "columns, need", m + 1L, "- rebuilding\n")
  }
  ## Several traits can share this big.matrix READ-only quite happily -- it is
  ## memory-mapped, so concurrent runs cost one copy of it between them.  Two
  ## processes REBUILDING it at the same time would write over each other and
  ## silently corrupt the genotypes, so that part is serialised with a lock.
  lock <- file.path(BIGMAT_DIR, ".rebuild.lock")
  if (rebuild) {
    if (file.exists(lock))
      stop("another run is rebuilding the big.matrix in ", BIGMAT_DIR, "\n",
           "  (lock file: ", lock, ")\n",
           "  Wait for it to finish, then start this one. If no run is active,\n",
           "  the previous one died -- delete the lock and rerun.")
    writeLines(c(format(Sys.time()), Sys.info()[["nodename"]], as.character(Sys.getpid())), lock)
    on.exit(unlink(lock), add = TRUE)
  }
  if (rebuild) {
    if (!file.exists(file.path(BIGMAT_DIR, "GD.int8.bin")))
      stop("GD.int8.bin is missing and GD.desc is unusable -- rerun 01_vcf_to_bigmatrix.py")
    for (f in c("GD.desc", "GD.bin")) {
      fp <- file.path(BIGMAT_DIR, f); if (file.exists(fp)) file.remove(fp)
    }
    cat("filling the big.matrix (one-off, a few minutes)...\n")
    options(bigmemory.typecast.warning = FALSE)
    GD <- filebacked.big.matrix(nrow = n, ncol = m + 1L, type = "char",
                                backingpath = BIGMAT_DIR, backingfile = "GD.bin",
                                descriptorfile = "GD.desc")
    GD[, 1] <- 0L
    con <- file(file.path(BIGMAT_DIR, "GD.int8.bin"), "rb"); chunk <- 20000L
    for (s in seq(1L, m, by = chunk)) {
      k <- min(chunk, m - s + 1L)
      v <- readBin(con, integer(), n = n * k, size = 1L, signed = TRUE)
      GD[, (s + 1L):(s + k)] <- matrix(v, nrow = n)
      cat(sprintf("\r  %d/%d", s + k - 1L, m)); utils::flush.console()
    }
    close(con); cat("\n"); unlink(lock)
  }
  stopifnot(ncol(GD) == m + 1L, nrow(GD) == n)

  ## phenotype and covariates in the genotype taxa order
  ord <- data.table(Taxa = taxa)
  Yf <- Y[ord, on = "Taxa"]
  Pf <- P_all[ord, on = "Taxa"]
  stopifnot(identical(Yf$Taxa, taxa), identical(Pf$Taxa, taxa))
  if (anyNA(Pf[, -1]))
    stop("the PCA file has no PCs for ", sum(!complete.cases(Pf[, -1])),
         " of the ", n, " genotyped taxa. FarmCPU needs covariates for every\n",
         "  row of the genotype matrix. Either extend the PCA file or drop those\n",
         "  taxa from the big.matrix.")
  cat(sprintf("FarmCPU: %d of %d taxa phenotyped\n", sum(!is.na(Yf[[TAG]])), n))

  ## FarmCPU() directly, not GAPIT(): GAPIT.Genotype() assumes column 1 of GD is
  ## taxa names and would both eat marker 1 and materialise the whole matrix.
  FarmCPU_fun <- tryCatch(get("FarmCPU", envir = asNamespace("GAPIT")),
                   error = function(e) tryCatch(get("FarmCPU"),
                   error = function(e2) stop("cannot find the FarmCPU function in GAPIT")))

  fwd <- file.path(OUTDIR, "farmcpu"); dir.create(fwd, recursive = TRUE, showWarnings = FALSE)
  owd2 <- getwd(); setwd(fwd)
  CVmat <- as.matrix(Pf[, -1, with = FALSE])
  stopifnot(nrow(CVmat) == n, ncol(CVmat) == NPC, !anyNA(CVmat))

  t0 <- Sys.time()
  myFarm <- FarmCPU_fun(
    Y = as.data.frame(Yf), GD = GD, GM = as.data.frame(GM), CV = CVmat,
    method.bin = "optimum", bin.size = c(5e5, 5e6, 5e7),
    bin.selection = seq(10, 100, 10), maxLoop = 10,
    QTN.threshold = 0.01, MAF.calculate = FALSE, file.output = FALSE)
  setwd(owd2)
  cat(sprintf("FarmCPU finished in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  saveRDS(myFarm, file.path(fwd, sprintf("FarmCPU_%s_GENOMEWIDE.rds", TAG)))
  res <- as.data.table(myFarm$GWAS)
  setnames(res, c("Chromosome", "Position", "p.value"), c("Chr", "Pos", "P.value"),
           skip_absent = TRUE)
  if (!"MAF" %in% names(res)) res[, MAF := NA_real_]
  if (!"nobs" %in% names(res)) res[, nobs := sum(!is.na(Yf[[TAG]]))]
  res[, `H&B.P.Value` := p.adjust(P.value, method = "BH")]
  if (!"effect" %in% names(res) && !"Effect" %in% names(res)) res[, Effect := NA_real_]
  if ("effect" %in% names(res)) setnames(res, "effect", "Effect")
  fwrite(res[, .(SNP, Chr, Pos, P.value, MAF, nobs, `H&B.P.Value`, Effect)],
         file.path(OUTDIR, sprintf("raw_GWAS_FarmCPU_3PC_%s.txt", TAG)), sep = "\t")

  sq <- myFarm$seqQTN
  if (length(sq)) {
    qt <- cbind(data.table(rank = seq_along(sq)), GM[sq])
    qt$P.value <- res$P.value[sq]
    fwrite(qt, file.path(OUTDIR, sprintf("FarmCPU_pseudoQTN_%s.csv", TAG)))
    cat("pseudo-QTNs retained: ", length(sq), " across ",
        length(unique(qt$Chromosome)), " chromosomes\n", sep = "")
  } else cat("FarmCPU retained no pseudo-QTNs -- it reduced to a GLM.\n")
}

## --------------------------------------- collect GAPIT output into raw files --
## GAPIT writes GAPIT.Association.GWAS_Results.<model>.<trait>.csv into whatever
## directory it ran in, so there is one per chromosome.
##
## It can write the SAME scan more than once per chromosome under different
## suffixes -- "<trait>(NYC)" and "<trait>(Kansas)" -- which are GAPIT's internal
## labels, not separate analyses: identical p-values, differing only in whether
## the effect column is spelled "Effect" or "effect".  Taking both would double
## every chromosome and leave the effect column half empty, so exactly one file
## is chosen per chromosome and the choice is reported.
if (length(gapit_models)) {
  wd <- file.path(OUTDIR, "gapit")
  if (!dir.exists(wd)) {
    warning("no GAPIT output directory at ", wd)
  } else {
  cat("\ncollecting GAPIT output from ", wd, "\n", sep = "")
  for (mod in gapit_models) {
    hits <- list.files(wd, recursive = TRUE, full.names = TRUE,
                       pattern = sprintf("GWAS_Results\\.%s\\..*\\.csv$", mod))
    if (!length(hits)) {
      if (!COLLECT_ONLY)
        warning("no GWAS result file found for ", mod, " under ", wd)
      next
    }
    ## one file per chromosome directory
    bydir <- split(hits, dirname(hits))
    chosen <- vapply(bydir, function(fs) {
      if (length(fs) == 1L) return(fs)
      ## prefer the copy with the canonical "Effect" column
      sc <- vapply(fs, function(f) {
        h <- tryCatch(names(fread(f, nrows = 0)), error = function(e) character(0))
        as.numeric("Effect" %in% h)
      }, numeric(1))
      fs[which.max(sc)]
    }, character(1))
    ndup <- length(hits) - length(chosen)
    cat(sprintf("  %-8s %d chromosome dir(s)%s\n", mod, length(chosen),
                if (ndup > 0) sprintf(", %d duplicate copy/copies ignored", ndup) else ""))

    parts <- lapply(chosen, function(f) {
      x <- fread(f)
      setnames(x, names(x), sub("^effect$", "Effect", names(x)))
      setnames(x, names(x), sub("^Chromosome$", "Chr", names(x)))
      setnames(x, names(x), sub("^Position *$", "Pos", names(x)))
      x
    })
    r <- rbindlist(parts, fill = TRUE)
    ## a SNP appearing twice now would mean two chromosome folders overlapped,
    ## which is a real problem rather than a naming quirk -- say so.
    dup <- sum(duplicated(r$SNP))
    if (dup > 0) {
      warning(mod, ": ", dup, " duplicated SNPs across chromosome folders -- ",
              "keeping the first of each, but check ", wd)
      r <- unique(r, by = "SNP")
    }
    for (need in c("MAF", "nobs", "Effect")) if (!need %in% names(r)) r[[need]] <- NA
    if (!"H&B.P.Value" %in% names(r)) r[["H&B.P.Value"]] <- p.adjust(r$P.value, method = "BH")
    setorder(r, Chr, Pos)
    out <- file.path(OUTDIR, sprintf("raw_GWAS_%s_3PC_%s.txt", mod, TAG))
    fwrite(r[, c("SNP","Chr","Pos","P.value","MAF","nobs","H&B.P.Value","Effect"), with = FALSE],
           out, sep = "\t")
    msg <- if (nrow(r) == EXPECTED_SNPS) "OK" else
      sprintf("*** %d markers, expected %d -- %d MISSING ***",
              nrow(r), EXPECTED_SNPS, EXPECTED_SNPS - nrow(r))
    cat(sprintf("           -> %s  %s\n", basename(out), msg))
    cat(sprintf("              effect column populated for %d of %d rows\n",
                sum(!is.na(r$Effect)), nrow(r)))
    if (nrow(r) != EXPECTED_SNPS) {
      cat("    markers present per chromosome:\n"); print(r[, .N, by = Chr][order(Chr)])
      cat("    (this is the check that would have caught MLMM dropping all of\n",
          "     chromosome 2 in the published soil-N run.)\n", sep = "")
    }
    ## quick read on the scan, so the pipeline check needs no extra step
    PB <- 0.05 / EXPECTED_SNPS; PS <- 1 / EXPECTED_SNPS
    cat(sprintf("              Bonferroni %.3e: %d SNPs | suggestive %.3e: %d SNPs\n",
                PB, sum(r$P.value <= PB, na.rm = TRUE),
                PS, sum(r$P.value <= PS, na.rm = TRUE)))
  }
  }
}

cat("\nfinished: ", format(Sys.time()), "\n", sep = "")
sink()
cat("next:  GN_DIR_GAPIT=", OUTDIR,
    " GN_GAPIT_TMPL='raw_GWAS_%s_3PC_", TAG, ".txt' \\\n",
    "       Rscript soilN_crossmodel_consensus.R\n", sep = "")
