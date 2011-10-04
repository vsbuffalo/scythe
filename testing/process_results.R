## process_results.R - process all results

## For the sim-reads file at the 0.5 contamination level
source("cutadapt/parse_results.R")
source("scythe/parse_results.R")

sens.scythe <- d.scythe[d.scythe$field == "sensitivity", 'value']
spec.scythe <- d.scythe[d.scythe$field == "specificity", 'value']
priors.scythe <- d.scythe[d.scythe$field == "specificity", 'error']

plot(1-spec.scythe, sens.scythe, xlim=c(0, 1), ylim=c(0, 1), type="n")
text(1-spec.scythe, sens.scythe, labels=priors.scythe, col="blue")

sens.cutadapt <- d.cutadapt[d.cutadapt$field == "sensitivity", 'value']
spec.cutadapt <- d.cutadapt[d.cutadapt$field == "specificity", 'value']
error.cutadapt <- d.cutadapt[d.cutadapt$field == "specificity", 'error']
text(1-spec.cutadapt, sens.cutadapt, labels=error.cutadapt, col="green")

abline(a=0, b=1, col="red")
