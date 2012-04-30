## analysis.R -- read in sequences to analyze trimming results
## note: assumes POSIX directory delimiter ("/")

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

compareFastqTrim <- function(original, trimmed) {
  original.d <- data.frame(id=as.character(id(original)), read=as.character(sread(original)), stringsAsFactors=FALSE)
  trimmed.d <- data.frame(id=as.character(id(trimmed)), read=as.character(sread(trimmed)), stringsAsFactors=FALSE)

  # full outer join of reads
  x <- merge(original.d, trimmed.d, by='id', all=TRUE)
  colnames(x) <- c("id", "original.read", "trimmed.read")

  x[is.na(x)] <- ""
  # Now, get actual contamination rate per sequence
  x$n.contam <- sapply(x$id, extractContam)

  # get sequence length difference (due to trimming)
  x$n.trimmed <- with(x, nchar(original.read) - nchar(trimmed.read))
  x[, -c(2, 3)]
}


calcStats <- function(x, y) {
  stopifnot(length(x) == length(y))
  n <- length(x)
  cc <- sum(x == y)
  data.frame(n.correct=cc, n.incorrect=n-cc)
}

## Now compare trimmed reads with original.
sim.read.files <- list.files(sim.read.dir, recursive=TRUE)
trimmed.read.files <- local({
  f <- list.files(trimmers, recursive=TRUE, full.names=TRUE, pattern="trimmed-.*") 
  tmp <- gsub("([a-z]+)/results/([\\.\\d]+)/trimmed-(\\d+)-([\\.\\d]+)\\.fastq", "\\1;;;\\2;;;\\3;;;\\4",
              f, perl=TRUE)
  chunks <- strsplit(tmp, ";;;")
  x <- do.call(rbind, chunks)
  data.frame(file=f, trimmer=x[, 1], contam.rate=x[, 2], rep=x[, 3], parameter=x[, 4], stringsAsFactors=FALSE)
})


results <- with(trimmed.read.files,
                mapply(function(f, trimmer, contam, rep, parameter) {
                  ## get corresponding simulated reads
                  sr <- sim.reads[[contam]][[rep]]

                  ## get corresponding trimmed reads
                  tr <- trimmed.reads[[trimmer]][[contam]][[paste(parameter, rep, sep=";;;;")]]
                  
                  tmp <- compareFastqTrim(sr, tr)
                  out <- with(tmp, calcStats(n.trimmed, n.contam))
                  out$trimmer <- trimmer
                  out$contam.rate <- contam
                  out$rep <- rep
                  out$parameter <- parameter
                  out
                }, file, trimmer, contam.rate, rep, parameter, SIMPLIFY=FALSE))


### Analysis of results
# Currently the data is at the read-level. We need to summarize the
# results at the trimmer, contamination, and parameter level.
d <- do.call(rbind, results)
write.table(d, file="testing-results-table.txt", sep="\t", quote=FALSE, row.names=FALSE)
