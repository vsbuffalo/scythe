PROGRAM_NAME = scythe
VERSION = 0.93
CC = clang
CFLAGS = -Wall -pedantic -DVERSION=$(VERSION)
DEBUG = -g
#OPT = -O3
ARCHIVE = $(PROGRAM_NAME)_$(VERSION)

# Mac OS X - Linux may need different linking
LDFLAGS = -lz -lm

default: build

match.o: src/match.c src/scythe.h
	$(CC) $(CFLAGS) -c $?
scythe.o: src/scythe.c src/kseq.h src/scythe.h
	$(CC) $(CFLAGS) -c $?
util.o: src/util.c src/kseq.h src/scythe.h
	$(CC) $(CFLAGS) -c $?
prob.o: src/prob.c src/scythe.h
	$(CC) $(CFLAGS) -c $?
tests.o: src/tests/tests.c src/scythe.h
	$(CC) $(CFLAGS) -c $?

valgrind: build
	valgrind --leak-check=full --show-reachable=yes ./scythe -a solexa_adapters.fa test.fastq

test: match.o util.o prob.o tests.o
	$(CC) $(CFLAGS) $(LDFLAGS) $? -o tests && ./tests

testclean:
	rm -rf ./tests

clean:
	rm -rf *.o ./scythe *.dSYM

distclean: clean
	rm -rf *.tar.gz

dist:
	tar -zcf $(ARCHIVE).tar.gz src Makefile

build: match.o scythe.o util.o prob.o 
	$(CC) $(CFLAGS) $(LDFLAGS) $? -o scythe

debug:
	$(MAKE) build "CFLAGS=-Wall -pedantic -g -DDEBUG"


## targets for simulations
sim: simclean build
	mkdir sims
	python read_sim.py -a solexa_adapters.fa -f test.fastq -r > sims/random.fastq
	python read_sim.py -a solexa_adapters.fa -f test.fastq > sims/contam.fastq
	./scythe -n 4 -a solexa_adapters.fa -o sims/contam_out.fastq -m sims/contam_matches.txt sims/contam.fastq
	./scythe -n 4 -a solexa_adapters.fa -o sims/random_out.fastq -m sims/random_matches.txt sims/random.fastq
	perl sim_stat.pl > sims/summary_stats.txt

prior: simclean build 
	mkdir sims
	python read_sim.py -a solexa_adapters.fa -f test.fastq -r > sims/random.fastq
	python read_sim.py -a solexa_adapters.fa -f test.fastq > sims/contam.fastq
	./scythe -p 0.01 -n 4 -a solexa_adapters.fa -o sims/contam_01_out.fastq -m sims/contam_01_matches.txt sims/contam.fastq
	./scythe -p 0.05 -n 4 -a solexa_adapters.fa -o sims/contam_05_out.fastq -m sims/contam_05_matches.txt sims/contam.fastq
	./scythe -p 0.10 -n 4 -a solexa_adapters.fa -o sims/contam_10_out.fastq -m sims/contam_10_matches.txt sims/contam.fastq
	./scythe -p 0.15 -n 4 -a solexa_adapters.fa -o sims/contam_15_out.fastq -m sims/contam_15_matches.txt sims/contam.fastq
	./scythe -p 0.20 -n 4 -a solexa_adapters.fa -o sims/contam_20_out.fastq -m sims/contam_20_matches.txt sims/contam.fastq
	./scythe -p 0.25 -n 4 -a solexa_adapters.fa -o sims/contam_25_out.fastq -m sims/contam_25_matches.txt sims/contam.fastq
	./scythe -p 0.30 -n 4 -a solexa_adapters.fa -o sims/contam_30_out.fastq -m sims/contam_30_matches.txt sims/contam.fastq
	./scythe -p 0.35 -n 4 -a solexa_adapters.fa -o sims/contam_35_out.fastq -m sims/contam_35_matches.txt sims/contam.fastq
	./scythe -p 0.40 -n 4 -a solexa_adapters.fa -o sims/contam_40_out.fastq -m sims/contam_40_matches.txt sims/contam.fastq
	./scythe -p 0.45 -n 4 -a solexa_adapters.fa -o sims/contam_45_out.fastq -m sims/contam_45_matches.txt sims/contam.fastq
	./scythe -p 0.50 -n 4 -a solexa_adapters.fa -o sims/contam_50_out.fastq -m sims/contam_50_matches.txt sims/contam.fastq
	./scythe -p 0.55 -n 4 -a solexa_adapters.fa -o sims/contam_55_out.fastq -m sims/contam_55_matches.txt sims/contam.fastq
	./scythe -p 0.90 -n 4 -a solexa_adapters.fa -o sims/contam_90_out.fastq -m sims/contam_90_matches.txt sims/contam.fastq

	perl sim_stat.pl contam_01 > sims/summary_stats_01.txt
	perl sim_stat.pl contam_05 > sims/summary_stats_05.txt
	perl sim_stat.pl contam_10 > sims/summary_stats_10.txt
	perl sim_stat.pl contam_15 > sims/summary_stats_15.txt
	perl sim_stat.pl contam_20 > sims/summary_stats_20.txt
	perl sim_stat.pl contam_25 > sims/summary_stats_25.txt
	perl sim_stat.pl contam_30 > sims/summary_stats_30.txt
	perl sim_stat.pl contam_35 > sims/summary_stats_35.txt
	perl sim_stat.pl contam_40 > sims/summary_stats_40.txt
	perl sim_stat.pl contam_45 > sims/summary_stats_45.txt	
	perl sim_stat.pl contam_50 > sims/summary_stats_50.txt	
	perl sim_stat.pl contam_90 > sims/summary_stats_90.txt	

simclean:
	rm -rf sims
