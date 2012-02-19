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

__KS_GETC(gzread, BUFFER_SIZE)
__KS_GETUNTIL(gzread, BUFFER_SIZE)
__KSEQ_READ


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
    if (!aseq->seq.l || !aseq->name.l) {
      fprintf(stderr, "Blank FASTA header or sequence in adapters file.\n");
      exit(EXIT_FAILURE);
    }

    adapters[i].seq = (char *) xmalloc(seq_l*sizeof(char) + 1);
    strncpy(adapters[i].seq, aseq->seq.s, seq_l);
    header_l = aseq->name.l + aseq->comment.l + 1;
    adapters[i].name = (char *) xmalloc(header_l*sizeof(char) + 1);
    strncpy(adapters[i].name, aseq->name.s, aseq->name.l+1);

    if (aseq->comment.s) {
      strncat(adapters[i].name, " ", 1);
      strncat(adapters[i].name, aseq->comment.s, aseq->comment.l+1);
    }

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

void fprint_float_array(FILE *fp, const float *array, int n) {
  int i;
  fprintf(fp, "[");
  for (i = 0; i < n; i++) {
    if (i != n-1)
      fprintf(fp, "%.3f, ", array[i]);
    else
      fprintf(fp, "%.3f]", array[i]);
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

void fprint_uint_array(FILE *fp, const unsigned int *array, int n) {
  int i;
  fprintf(fp, "[");
  for (i = 0; i < n; i++) {
    if (i != n-1)
      fprintf(fp, "%d, ", array[i]);
    else
      fprintf(fp, "%d]", array[i]);
  }
}

int sum(const int *x, int n) {
  int i;
  int s = 0;
  for (i = 0; i < n; i++)
    s += x[i];
  return s;
}

void write_fastq(gzFile output_fp, kseq_t *seq, int add_tag, char *tag, int match_n) {
  /* Heng Li's kseq.h handles FASTQ headers as such: anything after
     the first space is put in the field "comment" of the kseq_t
     struct. This function writes a single FASTQ block and wraps
     simple fprintf such that we don't have to worry about whether
     there's a comment or not. This has a variety of arguments for
     different output options,
  */
  if (match_n > 0) {
    /* If match_n is the number of matches we have to trim by, so
       output trimmed FASTQ seq (and possible header if add_tag is
       true. */
    if (add_tag) {
      if (seq->comment.s)
        fprintf(output_fp, 
                "@%s %s%s-%d\n%.*s\n+%s %s%s-%d\n%.*s\n", seq->name.s, seq->comment.s, tag, match_n, 
                (int) seq->seq.l-match_n, seq->seq.s, seq->name.s, seq->comment.s, tag, match_n, 
                (int) seq->seq.l-match_n, seq->qual.s);
      else 
        fprintf(output_fp, 
                "@%s%s-%d\n%.*s\n+%s%s-%d\n%.*s\n", seq->name.s, tag, match_n, 
                (int) seq->seq.l-match_n, seq->seq.s, seq->name.s, tag, match_n, 
                (int) seq->seq.l-match_n, seq->qual.s);

    } else {
      if (seq->comment.s)
        fprintf(output_fp, 
                "@%s %s\n%.*s\n+%s %s\n%.*s\n", seq->name.s, seq->comment.s,
                (int) seq->seq.l-match_n, seq->seq.s, seq->name.s, seq->comment.s,
                (int) seq->seq.l-match_n, seq->qual.s);
      else
        fprintf(output_fp, 
                "@%s\n%.*s\n+%s\n%.*s\n", seq->name.s,
                (int) seq->seq.l-match_n, seq->seq.s, seq->name.s,
                (int) seq->seq.l-match_n, seq->qual.s);

    }
  } else { 
    if (seq->comment.s)
      fprintf(output_fp, 
              "@%s %s\n%s\n+%s %s\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->name.s, seq->comment.s, seq->qual.s);
  else
      fprintf(output_fp, 
              "@%s\n%s\n+%s\n%s\n", seq->name.s, seq->seq.s, seq->name.s, seq->qual.s);
  }
}
