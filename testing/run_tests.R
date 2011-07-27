## run tests and get results

prior.data <- lapply(seq(0, 1, 0.1), function(prior) {
  system(sprintf("../scythe -n 0 -a ../solexa_adapters.fa -p %f -m results/matches-prior-%f.txt -o /dev/null test.fastq", prior, prior), intern=TRUE)
  system(sprintf("python results_parse.py -m results/matches-prior-%f.txt -r test.fastq > results/results-prior-%f.txt", prior, prior), intern=TRUE)
  tmp <- t(read.table(sprintf("results/results-prior-%f.txt", prior), header=FALSE, sep="\t"))
  if (ncol(tmp) < 5)
    return(NULL)

  d <- as.numeric(t(tmp[2, ]))
  dim(d) <- c(1, length(d))
  colnames(d) <- tmp[1, ]
  as.data.frame(d)
})


## combine, ignoring incomplete entries
include <- !unlist(lapply(prior.data, is.null))
priors <- seq(0, 1, 0.1)[include]
prior.data <- cbind(prior=priors, do.call(rbind, prior.data[include]))

plot(x=1-prior.data[, 15], y=prior.data[, 14], ylim=c(0, 1), xlim=c(0, 1), xlab="1 - specificity (FPR)", ylab="sensitivity (TPR)", type='n')
text(x=1-prior.data[, 15], y=prior.data[, 14], ylim=c(0, 1), xlim=c(0, 1), xlab="1 - specificity (FPR)", labels=prior.data$prior, cex=0.6)
