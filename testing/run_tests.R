## run tests and get results
library(lattice)
## source read simulation functions
source('read-sim.R')

if (!file.exists("sim-reads"))
  dir.create("sim-reads")

contam.rates <- seq(0, 0.9, 0.1)


prior.data <- lapply(contam.rates, function(contam.rate) {
  priors <- seq(0, 1, 0.1)

  ## simulate reads
  ff <- file("test-sample.fastq", open="r")
  contaminateFASTQEntry(ff, outfile=sprintf("sim-reads/sim-reads-%f.fastq", contam.rate), contam.rate, adapter)

  out <- lapply(priors, function(prior) {
    system(sprintf("../scythe -t -n 0 -a ../solexa_adapters.fa -p %f -m results/matches-prior-%f.txt -o results/out-prior-%f.fastq sim-reads/sim-reads-%f.fastq", prior, prior, prior, contam.rate), intern=TRUE)

    cmd <- sprintf("python results_parse.py -m results/matches-prior-%f.txt -r sim-reads/sim-reads-%f.fastq -t results/out-prior-%f.fastq > results/results-prior-%f.txt", prior, contam.rate, prior, prior)
    ok <- system(cmd)
    stopifnot(ok == 0)
    tmp <- t(read.table(sprintf("results/results-prior-%f.txt", prior), header=FALSE, sep="\t"))
    browser()
    if (ncol(tmp) < 5)
      return(NULL)

    d <- as.numeric(t(tmp[2, ]))
    dim(d) <- c(1, length(d))
    colnames(d) <- tmp[1, ]
    as.data.frame(d)
  })
  include <- !unlist(lapply(out, is.null))
  priors <- seq(0, 1, 0.1)[include]
  d <- cbind(contam.rate=contam.rate, prior=priors, do.call(rbind, out[include]))

  return(d)
})

## Munge data together for all contamination rates and plot with lattice
prior.data <- do.call(rbind, prior.data)
prior.data$y <- prior.data[, 'sensitivity']
prior.data$x <- 1-prior.data[, 'specificity']
png("roc-curves.png", width=1000, height=1000, res=100)
p <- xyplot(y ~ x | factor(contam.rate), data=prior.data, ylim=c(0, 1), xlim=c(0, 1),
       panel=function(x, y, subscript, ...) {
         panel.text(x, y, label=prior.data[subscript, 'prior'], cex=0.6)
         panel.abline(a=0, b=1, col="purple")
       }, xlab="1 - specificity (FPR)", ylab="sensitivity (TPR)")
print(p)
dev.off()
