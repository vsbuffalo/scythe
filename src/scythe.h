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

#define MAX_ADAPTERS 1000
#define MATCH_SCORE 1
#define MISMATCH_SCORE -1

#define IS_FASTQ(quality_type) INTEGER(quality_type)[0] >= 0
#define min(a,b) ((a)>(b)?(b):(a))

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

typedef struct likelihood_set {
  double random;
  double contam;
} likelihood_set;

typedef struct posterior_set {
  int is_contam;
  double random;
  double contam;
} posterior_set;

typedef struct match {
  int *match, *match_pos; // match positions; TODO simplify
  int shift; // position adapter found, -1 if is_contam is 0
  int length; // length of adapter contamination, -1 if is_contam is 0
  int score; // match score
  int adapter_index; // index of which adapter was found
  float *p_quals; // array of qualities convered to probs
  posterior_set *ps; // posterior set
} match;

/* prob.c prototypes */
float *qual_to_probs(const char *, quality_type);
double p_data(const int *, float *, float, float, int);
posterior_set *posterior(const int *, float *, float, float, int);
likelihood_set *likelihood(const int *, float *, float, int);

/* util.c prototypes */
void *xmalloc(size_t);
adapter_array *load_adapters(gzFile);
void destroy_adapters(adapter_array *, int);
char *fmt_matches(const char *, const char *, const int *, const int);
void print_float_array(const float *, int);
void fprint_float_array(FILE *, const float *, int);
void print_int_array(const int *, int);
void print_uint_array(const unsigned int *, int);
void fprint_uint_array(FILE *, const unsigned int *, int);
extern int sum(const int *, size_t);
void write_fastq(FILE *, kseq_t *, int, int, int);
void print_summary(adapter_array *, float, int, int, int);

/* match.c prototypes */
int *score_sequence(int *, const char *, const char *, size_t);
match *find_best_match(const adapter_array *, const char *, float *, float, float, int);
void print_match(kseq_t *, match *, FILE *, const adapter_array *, quality_type);
void destroy_match(match *);


#endif /* SCYTHE_H */
