rfile=sim-reads/sim-reads-0.500000.fastq

for error in 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95
do
    cutadapt -e $error -a AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT $rfile > cutadapt/results/diff-error/trimmed-$error.fastq 2> /dev/null
    python compare.py -r $rfile -t cutadapt/results/diff-error/trimmed-$error.fastq > cutadapt/results/diff-error/compare-$error.txt 2> cutadapt/results/diff-error/compare-offsets-$error.txt
done