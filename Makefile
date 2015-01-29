PROGRAM_NAME = scythe
VERSION = 0.994
CC = gcc
DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CFLAGS += -g -O0
else
	CFLAGS += -O3
endif

ARCHIVE = $(PROGRAM_NAME)_$(VERSION)
override CFLAGS += -Wall -pedantic -DVERSION=$(VERSION) -std=gnu99 -fPIC
override LDFLAGS += -lz -lm
LDTESTFLAGS = -lcheck
SDIR = src
LOBJS = match.o util.o prob.o
OBJS = $(LOBJS) scythe.o

.PHONY: clean distclean dist testclean lib test all debian debian-clean

all: scythe test-scythe libscythe.so

# Executables
scythe: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o scythe $(LDFLAGS)

test-scythe: $(LOBJS) tests.o
	$(CC) $(CFLAGS) $(LOBJS) tests.o -o test-scythe $(LDFLAGS) $(LDTESTFLAGS)

%.o: $(SDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

match.o: $(SDIR)/scythe.h
scythe.o: $(SDIR)/kseq.h $(SDIR)/scythe.h
util.o: $(SDIR)/kseq.h $(SDIR)/scythe.h
prob.o: $(SDIR)/scythe.h

# special case, source is not ./src/tests.c
tests.o: $(SDIR)/scythe.h src/tests/tests.c
	$(CC) $(CFLAGS) -c src/tests/tests.c -o $@

valgrind: scythe
	valgrind --leak-check=full --show-reachable=yes ./scythe -a solexa_adapters.fa test.fastq

test: test-scythe
	./test-scythe

testclean:
	rm -rf ./tests ./test-scythe

clean:
	rm -rf *.o ./scythe ./test-scythe ./libscythe.so *.dSYM

distclean: clean
	rm -rf *.tar.gz

dist:
	tar -zcf $(ARCHIVE).tar.gz src Makefile illumina_adapters.fa


lib: libscythe.so

libscythe.so: $(LOBJS)
	$(CC) $(CFLAGS) -shared -o $@ $^ $(LDFLAGS)

debian:
	mkdir -p scythe-debian
	cp -r debian src Makefile illumina_adapters.fa scythe-debian
	tar -zcf $(ARCHIVE).orig.tar.gz src Makefile illumina_adapters.fa

debian-clean:
	rm -f scythe_*.deb scythe*.dsc scythe_*.build scythe_*.changes  scythe_*.debian.tar.*  scythe_*.orig.tar.gz
	rm -rf scythe-debian
	$(MAKE) clean distclean
