## Contaminated read simulation

NUCLEOTIDES <- c('A', 'T', 'C', 'G')
#adapters <- list("solexa_adapter_1"="AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG")
adapter <- "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"

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
    ll <- runif(1, min.contam, nchar(adapter))
    contam <- substr(adapter, 1, ll)
    pre.error.seq <- paste(paste(sample(NUCLEOTIDES, length-ll+1, replace=TRUE), collapse=""), contam, sep="")
    return(addErrors(pre.error.seq, quality))
  }
}

contaminateFASTQEntry <- function(con, outfile, rate, adapters, min.contam=3) {
  outlist <- list()
  while (length(block <- readLines(con, n=4, warn=FALSE)) > 0) {
    print(block)
    if (runif(1) <= rate) {
      # contaminate
      header <- paste(block[1], "contaminated", sep="-")
      quality <- block[4]
      seq <- generateRandomSeq(nchar(block[2]), adapter=adapter, quality=quality, is.contam=TRUE, min.contam=min.contam)
      outlist <- c(outlist, header, seq, header, quality)
    } else {
      header <- block[1]
      quality <- block[4]
      seq <- generateRandomSeq(nchar(block[2]), is.contam=FALSE, min.contam=min.contam)
      outlist <- c(outlist, header, seq, header, quality)
    }
  }
  browser()
  writeLines(unlist(outlist), con=file(outfile))
  TRUE
}

ff <- file("small.fastq", open="r")
contaminateFASTQEntry(ff, outfile="test.fastq", 0.8, adapter)
