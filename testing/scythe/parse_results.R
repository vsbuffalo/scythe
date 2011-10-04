## results_parse.R - parse the multiple files from compare.py for scythe
dir <- "scythe/results/diff-contam/"
compare.files <- list.files(dir, pattern="compare-0")
priors <- as.numeric(sub("compare-(0\\.[0-9]+)\\.txt", "\\1", compare.files))

d <- mapply(function(f, e) {
  tmp <- read.table(file.path(dir, f), header=FALSE, sep="\t", col.names=c("field", "value"))
  cbind(tmp, error=e)
}, compare.files, priors, SIMPLIFY=FALSE)

d.scythe <- do.call(rbind, d)

## sens <- d[d$field == "sensitivity", 'value']
## spec <- d[d$field == "specificity", 'value']

## plot(1-spec, sens, xlim=c(0, 1), ylim=c(0, 1))
