## analysis.R -- read in sequences to analyze trimming results

library(Biostrings)
library(ShortRead)

# note: we use ShortRead because Biostrings can't:
#
# (1) read variable length FASTQ
# (2) drops quality info

### Settings
# Currently, we're using simulated reads with a contamination rate of
# 40%
sim.contam.rate <- 0.4

# Here are the trimmer we are testing; the names refer to the output
# directories of these trimmers
trimmers <- c("scythe", "btrim", "cutadapt")

# Simulatd read directories
sim.read.dir <- "simulated-reads"

### FASTQ file reading
# This is designed for handling many contamination rates, even though
# we're currently using only one.
sim.reads <- lapply(sim.contam.rate, function(contam) {
  # list replicate FASTQ files
  reps <- list.files(file.path(sim.read.dir, contam), full.names=TRUE)

  # extract rep name
  rep.names <- sapply(basename(reps), function(n) {
    gsub("(\\d+)\\.fastq", "\\1", n)
  })
  
  # read in FASTQ files
  out <- lapply(reps, function(f) {
    readFastq(f)
  })
  names(out) <- rep.names
  out
})
names(sim.reads) <- sim.contam.rate

### Read in trimmed FASTQ files
trimmed.results.dir <- "results"
trimmed.reads <- lapply(trimmers, function(trimmer) {

  # iterate over contam rates for each trimmer, gather files into list
  found.contam.rates <- list.files(file.path(trimmer, trimmed.results.dir), full.name=TRUE)
  
  out <- lapply(found.contam.rates, function(path) {
    ff <- list.files(path, pattern="trimmed-.*", full.names=TRUE)
    trimmed.md <- local({
      tmp <- gsub("trimmed-(\\d+)-(0?\\.?\\d+)\\.fastq", "\\2;;;;\\1",
                  basename(ff))
      reads <- lapply(ff, readFastq)
      names(reads) <- tmp
      reads
    })
  })

  names(out) <- basename(found.contam.rates)
  out
})
names(trimmed.reads) <- trimmers

### Main Analysis

extractContam <-
# extract the contamination from the header.
function(x) {
  if (regexpr("uncontaminated", x) != -1)
    return(0)
  return(as.numeric(gsub(".*-contaminated-([0-9]+)", "\\1", x)))
}

seqWidthDiff <- function(x, y) {
  width(x) - width(y)
}

compareFastqTrim <- function(original, trimmed) {
  original.d <- DataFrame(id=id(original), read=sread(original))
  trimmed.d <- DataFrame(id=id(trimmed), read=sread(trimmed))

  # full outer join of reads
  x <- merge(original.d, trimmed.d, by='id', all=TRUE)
  # change types
  x <- with(x, DataFrame(id=as.character(id),
                         original.read=DNAStringSet(read.x),
                         trimmed.read=DNAStringSet(read.y)))

  # Now, get actual contamination rate per sequence
  x$contam <- sapply(x$id, extractContam)

  # get width diff
  x$n.trimmed <- with(x, width(original.read) - width(trimmed.read))
  x[, -c(2, 3)]
}



