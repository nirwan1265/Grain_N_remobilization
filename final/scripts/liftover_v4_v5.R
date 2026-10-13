################################################################################
## B73 AGPv4 -> Zm-B73-REFERENCE-NAM-5.0 point lift-over from a UCSC chain file.
##
## The Goodman-Buckler amino acid GWAS used the HapMap 3.2.1 282-panel genotypes
## (hmp321_282_agpv4), so its 'Pos' column is AGPv4 (the SNP names carry the
## older AGPv3 coordinate). Genes are B73 RefGen_v5, so positions are lifted
## before any gene is assigned.
##
## Rule: same-chromosome '+' chains only; a position inside an aligned block of
## several chains takes the highest-scoring one; a position in an alignment gap
## moves to the nearest aligned base of the nearest block (flag "gap").
## source() this file, then lift_v4_v5(chr, pos) -> data.table(pos_v5, lift).
################################################################################
read_chain_blocks <- function(chain_file) {
  L <- readLines(chain_file)
  out <- vector("list", 0); k <- 0L
  i <- 1L; n <- length(L)
  while (i <= n) {
    if (startsWith(L[i], "chain")) {
      f <- strsplit(L[i], " ", fixed = TRUE)[[1]]
      score <- as.numeric(f[2]); tchr <- f[3]; tstart <- as.numeric(f[6])
      qchr <- f[8]; qstrand <- f[10]; qstart <- as.numeric(f[11])
      keep <- tchr == qchr && qstrand == "+" && grepl("^[0-9]+$", tchr)
      j <- i + 1L
      while (j <= n && nzchar(L[j]) && !startsWith(L[j], "chain")) j <- j + 1L
      if (keep) {
        m <- do.call(rbind, lapply(strsplit(L[(i + 1L):(j - 1L)], "[ \t]+"), function(x) as.numeric(c(x, 0, 0)[1:3])))
        size <- m[, 1]; dt <- m[, 2]; dq <- m[, 3]
        ts <- tstart + c(0, cumsum(size + dt)[-length(size)])
        qs <- qstart + c(0, cumsum(size + dq)[-length(size)])
        k <- k + 1L
        out[[k]] <- data.table::data.table(chr = as.integer(tchr), tS = ts, tE = ts + size, qS = qs, score = score)
      }
      i <- j
    } else i <- i + 1L
  }
  B <- data.table::rbindlist(out); data.table::setorder(B, chr, tS); B
}

lift_v4_v5 <- function(chr, pos, blocks) {
  res <- data.table::data.table(i = seq_along(pos), chr = as.integer(chr), x = pos - 1)  # 0-based
  res[, `:=`(pos_v5 = NA_real_, lift = NA_character_)]
  for (c in unique(res$chr)) {
    B <- blocks[chr == c]; if (!nrow(B)) { res[chr == c, lift := "no_chain"]; next }
    idx <- res[chr == c, i]; x <- res$x[idx]
    ## candidate blocks containing x: blocks are sorted by tS; check those with tS <= x and tE > x
    j <- findInterval(x, B$tS)
    best <- rep(NA_real_, length(x)); bscore <- rep(-Inf, length(x))
    for (off in 0:30) {                         # overlapping chains: look back a few blocks
      jj <- j - off; ok <- jj >= 1
      hit <- ok & B$tS[pmax(jj, 1)] <= x & B$tE[pmax(jj, 1)] > x & B$score[pmax(jj, 1)] > bscore
      best[hit] <- B$qS[jj[hit]] + (x[hit] - B$tS[jj[hit]]); bscore[hit] <- B$score[jj[hit]]
    }
    gap <- is.na(best)
    if (any(gap)) {
      jg <- pmax(j[gap], 1); jn <- pmin(jg + 1, nrow(B)); xg <- x[gap]
      dl <- abs(xg - (B$tE[jg] - 1)); dr <- abs(B$tS[jn] - xg)
      useR <- dr < dl
      best[gap] <- ifelse(useR, B$qS[jn], B$qS[jg] + (B$tE[jg] - B$tS[jg] - 1))
    }
    res[idx, `:=`(pos_v5 = best + 1, lift = ifelse(gap, "gap", "ok"))]
  }
  res[, .(pos_v5, lift)]
}
