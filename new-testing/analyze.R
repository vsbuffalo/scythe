## analyze.R -- analyze the results from trimming of simulated reads.


## In the cross-tool comparison, we are using 0.4 contamination
scythe.dir <- "scythe/results/0.4/"
cutadapt.dir <- "cutadapt/results/0.4/"
btrim.dir <- "btrim/results/0.4/"

scythe.result.files <- list.files(scythe.dir, pattern="results-.*\\.txt")
cutadapt.result.files <- list.files(cutadapt.dir, pattern="results-.*\\.txt")
btrim.result.files <- list.files(btrim.dir, pattern="results-.*\\.txt")

tmp.d <- lapply(scythe.result.files, function(f) {
  stopifnot(file.exists(file.path(cutadapt.dir, f))) ## should have corresponding files
  tmp <- gsub("results-(\\d+)-(0.\\d+).txt", "\\1;;;\\2", f)
  tmp <- unlist(strsplit(tmp, ";;;"))
  rep <- as.numeric(tmp[1])
  contam <- as.numeric(tmp[2])
  d.scythe <- read.table(file.path(scythe.dir, f), sep="\t", col.names=c("name", "value"), stringsAsFactors=FALSE)
  d.cutadapt <- read.table(file.path(cutadapt.dir, f), sep="\t", col.names=c("name", "value"), stringsAsFactors=FALSE)
  tmp <- rbind(cbind(d.cutadapt, rep, contam, trimmer="cutadapt"), cbind(d.scythe, rep, contam, trimmer="scythe"))
  tmp
})

# btrim has diff input parameter - num mismatches, so it needs to be handled seperately
tmp.btrim.d <- lapply(btrim.result.files, function(f) {
  tmp <- gsub("results-(\\d+)-(\\d+).txt", "\\1;;;\\2", f)
  tmp <- unlist(strsplit(tmp, ";;;"))
  rep <- as.numeric(tmp[1])
  num.mm <- as.numeric(tmp[2])
  d.btrim <- read.table(file.path(btrim.dir, f), sep="\t", col.names=c("name", "value"), stringsAsFactors=FALSE)
  if (nrow(d.btrim))
    tmp <- rbind(cbind(d.btrim, rep, contam=num.mm, trimmer="btrim"))
  else
    tmp <- NULL
  tmp
})

d <- rbind(do.call(rbind, tmp.d), do.call(rbind, tmp.btrim.d))


## just get TPR/FPR info
roc.d <- local({
  tmp <- d[d$name == "sensitivity", c("trimmer", "value", "rep", "contam")]
  tmp <- cbind(d[d$name == "sensitivity", c("trimmer", "value", "rep", "contam")], d[d$name == "specificity", c("trimmer", "value", "rep", "contam")])
  colnames(tmp) <- c("trimmer", "tpr", "rep", "contam", "trimmer", "specificity", "rep", "contam")
  tmp$fpr <- 1 - tmp$specificity
  tmp
})


# with all sim data with contam = 0.4
trimmer.colors <- c(scythe="blue", cutadapt="green", btrim="purple")
plot(tpr ~ fpr, data=roc.d, xlim=c(0, 1), ylim=c(0, 1), type="n")
for (t in levels(roc.d$trimmer)) {
  for (r in unique(roc.d$rep)) {
    # remove useless 0, 0 points
    with(roc.d[-which(roc.d$tpr == 0 & roc.d$fpr == 0), ], {
      lines(fpr[trimmer==t & rep==r], tpr[trimmer==t & rep==r], col=trimmer.colors[t])
      ## text(fpr[trimmer==t & rep==r], tpr[trimmer==t & rep==r], col=trimmer.colors[t], labels=contam[trimmer==t & rep==r])
    })
  }
}

# with just one rep, and priors/errors
trimmer.colors <- c(scythe="blue", cutadapt="green")
plot(tpr ~ fpr, data=roc.d, xlim=c(0, 1), ylim=c(0, 1), type="n")
for (t in levels(roc.d$trimmer)) {
  with(roc.d[roc.d$trimmer==t & roc.d$rep==1 & roc.d$fpr != 0 & roc.d$tpr != 0, ], {
    text(fpr[trimmer==t], tpr[trimmer==t], col=trimmer.colors[t], labels=contam[trimmer==t])
  })
}
stop()

## analyze offset data
scythe.offset.files <- list.files(scythe.dir, pattern="compare-offsets-")
cutadapt.offset.files <- list.files(cutadapt.dir, pattern="compare-offsets-")

offset.list <- lapply(scythe.offset.files, function(f) {
  stopifnot(file.exists(file.path(cutadapt.dir, f))) ## should have corresponding files
  tmp <- gsub("compare-offsets-(\\d+)-(0.\\d+).txt", "\\1;;;\\2", f)
  tmp <- unlist(strsplit(tmp, ";;;"))
  rep <- as.numeric(tmp[1])
  contam <- as.numeric(tmp[2])
  d.scythe <- read.table(file.path(scythe.dir, f), sep="\t", col.names=c("offset", "count"), stringsAsFactors=FALSE)
  d.cutadapt <- read.table(file.path(cutadapt.dir, f), sep="\t", col.names=c("offset", "count"), stringsAsFactors=FALSE)
  if (nrow(d.scythe) & nrow(d.cutadapt))
    rbind(cbind(d.scythe, rep, contam, total=sum(d.scythe$count), trimmer="scythe"), cbind(d.cutadapt, rep, contam, total=sum(d.cutadapt$count), trimmer="cutadapt"))
  else
    NULL
})

d.offsets <- do.call(rbind, offset.list)
