/* 
   prob.c - probability calculations used in scythe.
*/ 

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "scythe.h"

float *qual_to_probs(const char *qual, quality_type q_type) {
  /* 
     Converts qualities to the probability the base is correct.
     There are two quality-to-probability conversion formulas.
      - Solexa (pre-1.3 pipeline):
        q = -10 log10 (p/(1 - p))
        p = 1/(1 + 10^(q/10))

      - Illumina/Sanger/PHRED:
        q = -10 log10 p
        p = 1/(10^(q/10))

     From http://en.wikipedia.org/wiki/FASTQ.
  */
  int i, n = strlen(qual), q;
  float *probs = xmalloc(sizeof(float)*n);
  for (i = 0; i < n; i++) {
    q = (char) qual[i]-quality_contants[q_type][Q_OFFSET];
    if (q < quality_contants[q_type][Q_MIN] || q > quality_contants[q_type][Q_MAX]) {
      fprintf(stderr, "Base quality out of range for specified quality type (%d): %d\n", q_type, q);
      exit(EXIT_FAILURE);
    }
    if (q_type == SOLEXA) {
      probs[i] = 1 - 1/(1 + powf(10, (q/10.0)));
    } else {
      probs[i] = 1 - 1/(powf(10, (q/10.0)));
    }
  }
  return probs;
}

double p_data(const int *matches, float *p_quals, float p_prior_contam, float p_match, int n) {
  /* P(D), using the Law of Total Probability */
  double p_not_contam = 1.0, p_contam = 1.0;
  int i;
  for (i = 0; i < n; i++) {
    if (matches[i] == 1) {
      p_not_contam *= p_match;
      p_contam *= p_quals[i];
    } else {
      p_not_contam *= 1 - p_match;
      p_contam *= 1 - p_quals[i];
    }
  }
  return p_prior_contam*p_contam + (1-p_prior_contam)*p_not_contam;
}

posterior_set *posterior(const int *matches, float *p_quals, float p_prior, float p_match, int n) {
  /* The posterior for both the random model and the non-random
     model.
     
     `matches`: a int array of 0s and 1s indicating matching and
     mismatch positions.
     `p_quals`: quality values convered to probabilities.
     `p_prior`: prior probability
     `p_match`: probability of a match (for DNA: default is 0.25)
     `n`: match length.     
  */

  posterior_set *ps = xmalloc(sizeof(posterior_set));
  likelihood_set *ls = likelihood(matches, p_quals, p_match, n);
  double p_denom = p_data(matches, p_quals, p_prior, p_match, n);
  ps->contam = (ls->contam * p_prior)/p_denom;
  ps->random = (ls->random * (1-p_prior))/p_denom;
  ps->is_contam = ps->contam > ps->random;
  free(ls);
  return ps;  
}

likelihood_set *likelihood(const int *matches, float *p_quals, float p_match, int n) {
  /*
    Return a likelihood_set, which contains likelihoods of competing
    models: contamination and random matches.
  */
  likelihood_set *ls = xmalloc(sizeof(likelihood_set));
  int i;
  ls->random = 1;
  ls->contam = 1;
  for (i = 0; i < n; i++) {
    if (matches[i] == 1) {
      ls->random *= p_match; /* prob match happened by chance */
      ls->contam *= p_quals[i]; /* given contamination, prob seeing
                                   match despite base errors */
    } else {
      ls->random *= 1 - p_match;
      ls->contam *= 1 - p_quals[i];
    }
  }
  return ls;
}

