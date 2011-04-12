/*
  util.c - utility functions for scythe.
*/

#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "kseq.h"
#include "scythe.h"

#define MIN(a,b) ((a)>(b)?(b):(a))

KSEQ_INIT(gzFile, gzread)

void *xmalloc(size_t size) {
  void *ret = malloc(size);
  if (!ret) {
    fprintf(stderr, "Out of memory, malloc failed");
    exit(EXIT_FAILURE);
  }
  return(ret);
}

adapter_array *load_adapters(gzFile fp) {
  adapter_array *aa = xmalloc(sizeof(adapter_array));
  adapter *adapters = xmalloc(sizeof(adapter)*MAX_ADAPTERS);
  kseq_t *aseq;  
  int i=0, j, seq_l, header_l;

  aseq = kseq_init(fp);
  /* 
     Note: the string lengths from kseq.h structs do *not* include the
     null byte, so these have to be taken into account.
  */
  while ((seq_l = kseq_read(aseq)) >= 0) {
    adapters[i].seq = (char *) xmalloc(seq_l*sizeof(char) + 1);
    strncpy(adapters[i].seq, aseq->seq.s, seq_l);
    header_l = aseq->name.l + aseq->comment.l + 1;
    adapters[i].name = (char *) xmalloc(header_l*sizeof(char) + 1);
    strncpy(adapters[i].name, aseq->name.s, aseq->name.l+1);
    strncat(adapters[i].name, " ", 1);
    strncat(adapters[i].name, aseq->comment.s, aseq->comment.l+1);
    adapters[i].length = seq_l;
    
    /* occurences - for recording where adapters are found */
    adapters[i].occurrences = xmalloc(seq_l*sizeof(unsigned int));
    for (j=0; j < seq_l; j++)
      adapters[i].occurrences[j] = 0;
    i++;
  }

  aa->adapters = adapters;
  aa->n = i;

  kseq_destroy(aseq);
  return(aa);
}

void destroy_adapters(adapter_array *aa, int n) {
  int i;
  for (i = 0; i < aa->n; i++) {
    free((aa->adapters)[i].seq);
    free((aa->adapters)[i].name);
    free((aa->adapters)[i].occurrences);
  }
  free(aa->adapters);
  free(aa);
}


char *fmt_matches(const char *seqa, const char *seqb, const int *matches, const int n) {
  /* 
     Make a string of matches. xmalloc() needs three times the
     sequence lengths, four new line characters, and one null byte.
  */
  char *out = xmalloc(3*(n + 1)*sizeof(char) + 2), *ptr=out;
  int i;
  /* if (strlen(seqa) == n) { */
  /*   printf("seqa: %s\nseqb: %s\n", seqa, seqb); */
  /* } */
  sprintf(out, "%.*s\n", n, seqa);
  ptr += n + 1;
  for (i = 0; i < n; i++) {
    if (matches[i] == MATCH_SCORE) {
      sprintf(ptr, "|");
    } else {
      sprintf(ptr, " ");
    }
    ptr++;
  }
  sprintf(ptr, "\n%s", seqb);
  return out;
}

void print_float_array(const float *array, int n) {
  int i;
  printf("[");
  for (i = 0; i < n; i++) {
    if (i != n-1)
      printf("%f, ", array[i]);
    else
    printf("%f]", array[i]);
  }
}

void print_int_array(const int *array, int n) {
  int i;
  printf("[");
  for (i = 0; i < n; i++) {
    if (i != n-1)
      printf("%d, ", array[i]);
    else
    printf("%d]", array[i]);
  }
}


void print_uint_array(const unsigned int *array, int n) {
  int i;
  printf("[");
  for (i = 0; i < n; i++) {
    if (i != n-1)
      printf("%d, ", array[i]);
    else
    printf("%d]", array[i]);
  }
}

int sum(const int *x, int n) {
  int i;
  int s = 0;
  for (i = 0; i < n; i++)
    s += x[i];
  return s;
}
