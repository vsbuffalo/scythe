rfile=sim-reads/sim-reads-0.500000.fastq

if [ ! -d scythe/results/diff-contam ];
then
    echo "creating directory for scythe results"
    mkdir -p scythe/results/diff-contam
fi

for contam in 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95
do
    ../scythe -n 0 -p $contam -a ../solexa_adapters.fa $rfile > scythe/results/diff-contam/trimmed-$contam.fastq 2> /dev/null
    python compare.py -r $rfile -t scythe/results/diff-contam/trimmed-$contam.fastq > scythe/results/diff-contam/compare-$contam.txt 2> scythe/results/diff-contam/compare-offsets-$contam.txt
done