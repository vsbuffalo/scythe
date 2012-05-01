## run.sh -- run tests
if [ -z $1 ]
then
    contam_rate=0.4
else 
    contam_rate=$1
fi

mkdir -p scythe/results/$contam_rate
## scythe trimming - multiple priors for 40% contamination
cfile=simulated-reads/$contam_rate
for prior in 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95
do
    for rep in 01 02 03 04 05 06 07 08 09 10
    do
        ../scythe -n 0 -p $prior -a ../illumina_adapters.fa $cfile/$rep.fastq > scythe/results/$contam_rate/trimmed-$rep-$prior.fastq 2> /dev/null
    done
done

## cutadapt trimming - multiple error rates for 40% contamination
mkdir -p cutadapt/results/$contam_rate/
for error in 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95
do
    for rep in 01 02 03 04 05 06 07 08 09 10
    do
        cutadapt -e $error -a AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT $cfile/$rep.fastq > cutadapt/results/$contam_rate/trimmed-$rep-$error.fastq 2> /dev/null
    done
done


## btrim trimming - run in a mode closest to other trimmers:
## 
##  - no quality trimming
##  - 3'-end only, and the varying error is the
##    number of mismatches
mkdir -p btrim/results/$contam_rate/
for error in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19  ## Up to max error
do
    for rep in 01 02 03 04 05 06 07 08 09 10
    do
        btrim_mac -v $error -3 -a 0 -P -p btrim_adapters.fa -t $cfile/$rep.fastq -o btrim/results/$contam_rate/trimmed-$rep-$error.fastq > /dev/null
    done
done

