## Contaminated read simulation

NUCLEOTIDES <- c('A', 'T', 'C', 'G')
#adapters <- list("solexa_adapter_1"="AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG")
adapter <- "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"
#adapter <- "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

set.seed(0)

addErrors <-
# given a sequence, add errors based on an Illumina-encoded quality line
function(seq, quality) {
  # Note: I have checked this against
  # http://seqanswers.com/forums/showthread.php?t=1523 and it appears
  # correct - VB
  stopifnot(nchar(seq) == nchar(quality))

  probs.wrong <- sapply(charToRaw(quality), function(x) 1/(10^((as.integer(x) - 64)/10)))
  paste(mapply(function(base, p.wrong) {
    if (rbinom(1, 1, p.wrong))
      return(sample(setdiff(NUCLEOTIDES, base), 1))
    else
      return(base)
  }, unlist(strsplit(seq, '')), probs.wrong), collapse="")
}

generateRandomSeq <-
# generate a random sequence, possibly with contamination of an adapter, it's length uniform
function(length, adapter=NULL, quality=NULL, is.contam=FALSE, min.contam=3) {
  if (!is.contam)
    paste(sample(NUCLEOTIDES, length, replace=TRUE), collapse="")
  else {
    ll <- sample(min.contam:nchar(adapter), 1)
    contam <- substr(adapter, 1, ll)
    pre.error.seq <- paste(paste(sample(NUCLEOTIDES, length-nchar(contam), replace=TRUE), collapse=""), contam, sep="")
    return(list(contam.n=ll, seq=addErrors(pre.error.seq, quality)))
  }
}

contaminateFASTQEntry <- function(con, outfile, rate, adapters, min.contam=3) {
  blocks.processed <- 0
  reads <- readLines(con)
  outlist <- vector('character', length(reads))
  while (blocks.processed*4 < length(reads)) {
    quality <- reads[4*blocks.processed+4]
    seq <- reads[4*blocks.processed+2]
    header <- reads[4*blocks.processed+1]

    if (runif(1) <= rate) {
      # contaminate
      tmp <- generateRandomSeq(nchar(seq), adapter=adapter, quality=quality, is.contam=TRUE, min.contam=min.contam)
      seq <- tmp$seq
      ll <- tmp$contam.n
      header <- sprintf("%s-contaminated-%d", header, ll)
    } else {
      header <- sprintf("%s-uncontaminated", header)
      seq <- generateRandomSeq(nchar(seq), is.contam=FALSE, min.contam=min.contam)
    }

    # output results to vector
    outlist[4*blocks.processed+1] <- header
    outlist[4*blocks.processed+2] <- seq
    outlist[4*blocks.processed+3] <- sprintf("+%s", substr(header, 2, nchar(header)))
    outlist[4*blocks.processed+4] <- quality

    blocks.processed <- blocks.processed + 1
    if (blocks.processed %% 100 == 0)
      message(sprintf("%d blocks processed.", blocks.processed))
  }
  writeLines(outlist, con=file(outfile))
}

ff <- file("test-sample.fastq", open="r")
contaminateFASTQEntry(ff, outfile="test.fastq", 0.8, adapter)
