## generate-reads.R -- generate reads with different contamination
## rates.

source("lib/read-sim.R")

contams <- seq(0, 0.95, 0.05)
n.reps <- 10

if (!file.exists("simulated-reads"))
  dir.create("simulated-reads")

for (rate in contams) {
  for (rep in 1:n.reps) {
    dir <- file.path("simulated-reads", rate)
    if (!file.exists(dir)) {
      dir.create(dir)
    }
    message(sprintf("creating file %d for contamination rate %f.\n", rep, rate))
    outfile <- file.path(dir, sprintf("%02d.fastq", rep))
    message(outfile)
    ff <- file("reads.fastq", open="r")
    contaminateFASTQEntry(ff, outfile, rate, verbose=FALSE)
    close(ff)
  }
}
