clean:
	rm -rf scythe/results/*
	rm -rf cutadapt/results/*
	rm -rf btrim/results/*

clean-simulated-reads:
	rm -rf simulated-reads/*

run-trimmers:
	bash run.sh 0.4
	bash run.sh 0.7

generate-reads:
	R CMD BATCH generate-reads.R