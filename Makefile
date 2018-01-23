PROGRAM_NAME = scythe
VERSION = 0.991
CC = gcc
DEBUG ?= 0
CFLAGS = -Wall -pedantic -DVERSION=$(VERSION) -std=gnu99
ifeq ($(DEBUG), 1)
	CFLAGS += -g -O0
else 
	CFLAGS += -O3
endif
ARCHIVE = $(PROGRAM_NAME)_$(VERSION)
LDFLAGS = -lz -lm
LDTESTFLAGS = -lcheck
SDIR = src
OBJS = match.o scythe.o util.o prob.o 
LOBJS = match.o util.o prob.o 


.PHONY: clean default all distclean dist tests testclean lib

default: all

%.o: $(SDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

match.o: $(SDIR)/scythe.h
scythe.o: $(SDIR)/kseq.h $(SDIR)/scythe.h
util.o: $(SDIR)/kseq.h $(SDIR)/scythe.h
prob.o: $(SDIR)/scythe.h
test.o: $(SDIR)/scythe.h

valgrind: build
	valgrind --leak-check=full --show-reachable=yes ./scythe -a solexa_adapters.fa test.fastq

test: clean match.o util.o prob.o test.o
	$(CC) $(CFLAGS) $? -o test $(LDFLAGS) $(LDTESTFLAGS) && ./test

testclean:
	rm -rf ./tests

clean:
	rm -rf *.o ./scythe *.dSYM

distclean: clean
	rm -rf *.tar.gz

dist:
	tar -zcf $(ARCHIVE).tar.gz src Makefile

all: $(OBJS)
	$(CC) $(CFLAGS) $? -o scythe $(LDFLAGS)

lib: libscythe.so

libscythe.so: CFLAGS += -fpic
libscythe.so: $(LOBJS)
	$(CC) $(CFLAGS) -shared -o $@ $^ $(LDFLAGS)
