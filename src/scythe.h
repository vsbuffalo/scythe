#ifndef SCYTHE_H
#define SCYTHE_H

#include <zlib.h>


#define MAX_ADAPTERS 20
#define MATCH_SCORE 1
#define MISMATCH_SCORE -1

#define IS_FASTQ(quality_type) INTEGER(quality_type)[0] >= 0

typedef enum {
  PHRED, 
  SANGER,
  SOLEXA,
  ILLUMINA
} quality_type;

#define Q_OFFSET 0
#define Q_MIN 1
#define Q_MAX 2

/* The quality constants come from:
   http://www.biopython.org/DIST/docs/api/Bio.SeqIO.QualityIO-module.html*/

static const int quality_contants[4][3] = {
  /* offset, min, max */
  {0, 4, 60}, /* PHRED */
  {33, 0, 93}, /* SANGER */
  {64, -5, 62}, /* SOLEXA, early Illumina (pre-pipeline 1.3) */
  {64, 0, 62} /* ILLUMINA (post-pipeline 1.3) */
};

typedef struct adapter {
  char *name;
  char *seq;
  int length;
  unsigned int *occurrences;
} adapter;

typedef struct adapter_array {
  adapter *adapters;
  int n;
} adapter_array;

enum contam {
  NOT_CONTAM,
  CONTAM
};

typedef struct likelihood_set {
  double random;
  double contam;
} likelihood_set;

typedef struct posterior_set {
  enum contam is_contam;
  double random;
  double contam;
} posterior_set;

typedef struct match {
  int *match;
  int n;
  int score;
  int *match_pos;
  float *p_quals;
  posterior_set *ps;
  int adapter_index;
} match;

/* prob.c prototypes */
float *qual_to_probs(const char *qual, quality_type q_type);
double p_data(const int *matches, float *p_quals, float p_prior_contam, float p_match, int n);
posterior_set *posterior(const int *matches, float *p_quals, float p_prior, float p_match, int n);
likelihood_set *likelihood(const int *matches, float *p_quals, float p_match, int n);

/* util.c prototypes */
void *xmalloc(size_t size);
adapter_array *load_adapters(gzFile fp);
void destroy_adapters(adapter_array *aa, int n);
char *fmt_matches(const char *seqa, const char *seqb, const int *matches, const int n);
void print_float_array(const float *array, int n);
void print_int_array(const int *array, int n);
void print_uint_array(const unsigned int *array, int n);
int sum(const int *x, int n);

/* match.c prototypes */
int *score_sequence(const char *seqa, const char *seqb, int n);
match *find_best_match(const adapter_array *aa, const char *read, float *p_quals, float prior, float p_match, int min);
void destroy_match(match *m);


#endif /* SCYTHE_H */
