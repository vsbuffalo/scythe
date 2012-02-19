#ifndef SCYTHE_H
#define SCYTHE_H

#include <zlib.h>
#include "kseq.h"

/* KSEQ_INIT() cannot be called here, because we only need the types
   defined. Calling KSEQ_INIT() would also define functions, leading
   to an unused function warning with GCC. So, the basic typedefs
   kseq.h has are included here, and each file that reads needs:
   
   This was a fix also used in Nik Joshi's Sickle (which I also helped
   with) and is the only way I know of dealing with this.

   __KS_GETC(gzread, BUFFER_SIZE)
   __KS_GETUNTIL(gzread, BUFFER_SIZE)
   __KSEQ_READ
*/

#define BUFFER_SIZE 4096
__KS_TYPE(gzFile)
__KS_BASIC(gzFile, BUFFER_SIZE)
__KSEQ_TYPE(gzFile)
__KSEQ_BASIC(gzFile)

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
void fprint_float_array(FILE *fp, const float *array, int n);
void print_int_array(const int *array, int n);
void print_uint_array(const unsigned int *array, int n);
void fprint_uint_array(FILE *fp, const unsigned int *array, int n);
int sum(const int *x, int n);
void write_fastq(gzFile output_fp, kseq_t *seq, int add_tag, char *tag, int match_n);

/* match.c prototypes */
int *score_sequence(const char *seqa, const char *seqb, int n);
match *find_best_match(const adapter_array *aa, const char *read, float *p_quals, float prior, float p_match, int min);
void destroy_match(match *m);


#endif /* SCYTHE_H */
