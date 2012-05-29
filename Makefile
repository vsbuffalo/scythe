PROGRAM_NAME = scythe
VERSION = 0.981
CC = gcc
CFLAGS = -Wall -pedantic -DVERSION=$(VERSION) -std=gnu99
DEBUG = -g
OPT = -O3
ARCHIVE = $(PROGRAM_NAME)_$(VERSION)
LDFLAGS = -lz -lm
SDIR = src

.PHONY: clean default build distclean dist debug

default: build

match.o: $(SDIR)/match.c $(SDIR)/scythe.h
	$(CC) $(CFLAGS) -c $(SDIR)/$*.c
scythe.o: $(SDIR)/scythe.c $(SDIR)/kseq.h $(SDIR)/scythe.h
	$(CC) $(CFLAGS) -c $(SDIR)/$*.c
util.o: $(SDIR)/util.c $(SDIR)/kseq.h $(SDIR)/scythe.h
	$(CC) $(CFLAGS) -c $(SDIR)/$*.c
prob.o: $(SDIR)/prob.c $(SDIR)/scythe.h
	$(CC) $(CFLAGS) -c $(SDIR)/$*.c
tests.o: $(SDIR)/tests/tests.c $(SDIR)/scythe.h
	$(CC) $(CFLAGS) -c $(SDIR)/$*.c

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
	$(CC) $(LDFLAGS) $(DEBUG) -o scythe src/*.c