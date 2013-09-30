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

__KSEQ_BASIC(static, gzFile)
__KS_GETC(gzread, BUFFER_SIZE)
__KS_GETUNTIL(gzread, BUFFER_SIZE)
__KSEQ_READ(static)


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
  sprintf(out, "%.*s\n", n, seqa);
  ptr += n + 1;
  for (i = 0; i < n; i++) {
    if (matches[i]) {
      sprintf(ptr, "|");
    } else {
      sprintf(ptr, " ");
    }
    ptr++;
  }
  sprintf(ptr, "\n%.*s", (int) n, seqb);
  return out;
}

void print_float_array(const float *array, int n) {
  int i;
  printf("[");
  for (i = 0; i < n; i++) {
    if (i != n-1)
      printf("%.2f, ", array[i]);
    else
      printf("%.2f]", array[i]);
  }
}

void fprint_float_array(FILE *fp, const float *array, int n) {
  int i;
  fprintf(fp, "[");
  for (i = 0; i < n; i++) {
    if (i != n-1)
      fprintf(fp, "%.2f, ", array[i]);
    else
      fprintf(fp, "%.2f]", array[i]);
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

void write_fastq(FILE *output_fp, kseq_t *seq, int add_tag, int shift, int min_keep) {
  char tag[] = ";;cut_scythe";
  char *sequence, *qual;

  /* check if we have fewer than min-match */ 
  if (shift != -1 && shift <= min_keep) {
    sequence = "N";
    qual = "B";
    shift = -1;
  } else {
    sequence = seq->seq.s;
    qual = seq->qual.s;
  }
  if (shift >= 0) {
    
    if (add_tag) {
      fprintf(output_fp, 
              "@%s %.*s%s-%d\n%.*s\n+\n%.*s\n", seq->name.s, (int) seq->comment.l, 
              seq->comment.s, tag, shift, 
              (int) shift, sequence,
              (int) shift, qual);
    } else {
      fprintf(output_fp, 
              "@%s %.*s\n%.*s\n+\n%.*s\n", seq->name.s, (int) seq->comment.l,
              seq->comment.s,
              (int) shift, sequence,
              (int) shift, qual);
    }
  } else
    fprintf(output_fp, 
            "@%s %.*s\n%s\n+\n%s\n", seq->name.s, (int) seq->comment.l, 
            seq->comment.s, sequence, qual);
}

void print_summary(adapter_array *aa, float prior, int uncontaminated, 
                   int contaminated, int total) {
  /* int i; */
  fprintf(stderr, "prior: %0.3f\n", prior);
  fprintf(stderr, "\nAdapter Trimming Complete\ncontaminated: %d, uncontaminated: %d, total: %d\n", 
          contaminated, total-contaminated, total);
  fprintf(stderr, "contamination rate: %f\n\n", contaminated/(float) total);
  /* for (i = 0; i < aa->n; i++) { */
  /*   fprintf(stderr, "\nAdapter %d '%s' contamination occurences:\n", i+1, aa->adapters[i].name); */
  /*   fprint_uint_array(stderr, aa->adapters[i].occurrences, aa->adapters[i].length); */
  /*   fprintf(stderr, "\n"); */
  /* } */
}

void print_match(kseq_t *seq, match *match, FILE *matches_fp, 
                 const adapter_array *aa, quality_type qual_type) {
  /* Make a string that indicates the position of the matches with "|"s. */
  char *match_string;
  match_string = fmt_matches((aa->adapters[match->adapter_index]).seq,
                             &(seq->seq.s)[match->shift], 
                             match->match_pos, match->length);
  
  fprintf(matches_fp, "p(c|s): %f; p(!c|s): %f; adapter: %s\n%s\n%s\n%.*s\n", 
          match->ps->contam, match->ps->random,
          aa->adapters[match->adapter_index].name,
          seq->name.s, match_string, 
          (int) match->length,
          &(seq->qual.s)[match->shift]);
  fprint_float_array(matches_fp, qual_to_probs(&(seq->qual.s)[match->shift], qual_type), match->length);
  fprintf(matches_fp, "\n\n");
  free(match_string);
}
