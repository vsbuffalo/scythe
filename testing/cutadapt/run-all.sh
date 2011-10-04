simdir=./sim-reads/
readsfiles=$(ls ./sim-reads/)

if [ ! -d cutadapt/results ];
then
    echo "creating directory for cutadapt results"
    mkdir -p cutadapt/results/all
fi

for rfile in $readsfiles
do
    contam=$(echo $rfile | sed 's/sim-reads-\(.*\)\.fastq/\1/')
    echo $contam
    cutadapt -e 0.2 -a AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT $simdir/$rfile > cutadapt/results/trimmed-$contam.fastq 2> /dev/null
    python compare.py -r $simdir/$rfile -t cutadapt/results/trimmed-$contam.fastq > cutadapt/results/compare-$contam.txt
done




