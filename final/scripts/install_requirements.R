#!/usr/bin/env Rscript
################################################################################
## Everything GWAS_models.R needs, installed into a personal library.
##
##   Rscript install_requirements.R          # install what is missing
##   Rscript install_requirements.R check    # report only, install nothing
##
## GAPIT needs R >= 4.3.0.  On a cluster, load an R module that new before
## running this (module avail R).
################################################################################

CHECK_ONLY <- identical(commandArgs(trailingOnly = TRUE)[1], "check")

cat(R.version.string, "\n")
if (getRversion() < "4.3.0")
  cat("*** GAPIT's DESCRIPTION requires R >= 4.3.0 -- load a newer R module ***\n")

## a writable personal library, so nothing needs root
lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(lib)) lib <- file.path(path.expand("~"), "R", paste0(R.version$platform,
                                   "-library", substr(getRversion(), 1, 3)))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
cat("library: ", lib, "\n\n", sep = "")

CRAN_PKGS <- c(
  ## GAPIT's own Imports (jiabowang/GAPIT DESCRIPTION)
  "ape", "bigmemory", "EMMREML", "genetics", "gplots", "htmltools",
  "magrittr", "lme4", "MASS", "plotly", "RcppArmadillo", "scatterplot3d",
  "snowfall",
  ## used by GWAS_models.R itself
  "data.table", "biganalytics",
  ## only needed for the install_github route; harmless if it fails
  "remotes")
BIOC_PKGS <- c("multtest")        # GAPIT imports it; it is NOT on CRAN

have <- function(p) requireNamespace(p, quietly = TRUE)

miss_cran <- CRAN_PKGS[!vapply(CRAN_PKGS, have, logical(1))]
miss_bioc <- BIOC_PKGS[!vapply(BIOC_PKGS, have, logical(1))]
has_gapit <- have("GAPIT")

cat("CRAN missing : ", if (length(miss_cran)) paste(miss_cran, collapse = ", ") else "none", "\n")
cat("Bioc missing : ", if (length(miss_bioc)) paste(miss_bioc, collapse = ", ") else "none", "\n")
cat("GAPIT        : ", if (has_gapit) "installed" else "MISSING", "\n\n")

if (CHECK_ONLY) {
  ok <- !length(miss_cran) && !length(miss_bioc) && has_gapit
  cat(if (ok) "all requirements satisfied\n" else "something is missing -- rerun without 'check'\n")
  quit(status = if (ok) 0 else 1)
}

if (length(miss_cran)) {
  cat("installing from CRAN...\n")
  install.packages(miss_cran, lib = lib, repos = "https://cloud.r-project.org", Ncpus = 4)
}
if (length(miss_bioc)) {
  if (!have("BiocManager"))
    install.packages("BiocManager", lib = lib, repos = "https://cloud.r-project.org")
  cat("installing from Bioconductor...\n")
  BiocManager::install(miss_bioc, lib = lib, ask = FALSE, update = FALSE)
}
if (!has_gapit) {
  ## Three routes, tried in order.  Any one of them is enough.
  tarball <- Sys.getenv("GN_GAPIT_TARBALL")      # route 2, set to a local file
  done <- FALSE

  if (nzchar(tarball) && file.exists(tarball)) {
    cat("installing GAPIT from the local tarball: ", tarball, "\n", sep = "")
    install.packages(tarball, lib = lib, repos = NULL, type = "source")
    done <- have("GAPIT")
  }

  if (!done && have("remotes")) {
    cat("installing GAPIT from GitHub via remotes...\n")
    try(remotes::install_github("jiabowang/GAPIT", lib = lib,
                                dependencies = TRUE, upgrade = "never"), silent = FALSE)
    done <- have("GAPIT")
  }

  if (!done) {
    cat("\n*** GAPIT could not be installed automatically ***\n\n",
        "Pick whichever of these the server allows:\n\n",
        "1. remotes is just a CRAN package -- install it and rerun this script:\n",
        '     install.packages("remotes")\n\n',
        "2. No remotes, no devtools, no GitHub access from the server.\n",
        "   Download the tarball anywhere with a browser, copy it over, then:\n",
        "     curl -L -o GAPIT.tar.gz https://github.com/jiabowang/GAPIT/archive/refs/heads/master.tar.gz\n",
        "     scp GAPIT.tar.gz you@server:~/\n",
        "   and on the server either\n",
        "     tar xzf GAPIT.tar.gz && R CMD INSTALL -l ", lib, " GAPIT-master\n",
        "   or rerun this script with\n",
        "     GN_GAPIT_TARBALL=~/GAPIT.tar.gz Rscript install_requirements.R\n\n",
        "3. Do not install GAPIT at all -- source its functions instead.\n",
        "   Fetch gapit_functions.txt once, put it on the server, and run\n",
        "     export GN_GAPIT_SOURCE=~/gapit_functions.txt\n",
        "   GWAS_models.R picks that up and uses it. The CRAN and Bioconductor\n",
        "   packages above are still required for that route.\n", sep = "")
  }
}

## ---- verify ----------------------------------------------------------------
cat("\n---- verification ----\n")
for (p in c(CRAN_PKGS, BIOC_PKGS, "GAPIT"))
  cat(sprintf("  %-15s %s\n", p, if (have(p)) "ok" else "STILL MISSING"))

gsrc <- Sys.getenv("GN_GAPIT_SOURCE")
if (have("GAPIT") || nzchar(gsrc)) {
  if (have("GAPIT")) suppressPackageStartupMessages(library(GAPIT)) else {
    cat("sourcing GAPIT functions from ", gsrc, "\n", sep = ""); source(gsrc)
  }
  gf <- names(formals(GAPIT))
  cat("\nGAPIT() accepts a multi-file genotype interface (file.G/file.from/file.to): ",
      if (all(c("file.G", "file.from", "file.to") %in% gf)) "yes" else
        "no -- GWAS_models.R will use the per-chromosome loop", "\n", sep = "")
  cat("number of GAPIT() arguments: ", length(gf), "\n", sep = "")
}
