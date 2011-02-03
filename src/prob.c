/* 
   prob.c - probability calculations used in scythe.
*/ 

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "scythe.h"

float *qual_to_probs(const char *qual, quality_type q_type) {
  int i, n = strlen(qual);
  float *probs = xmalloc(sizeof(float)*n);
  for (i = 0; i < n; i++) {
    probs[i] = 1/(1 + powf(10, -((char) qual[i]-quality_contants[q_type][Q_OFFSET])/10.0));
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
  */

  posterior_set *ps = xmalloc(sizeof(posterior_set));
  likelihood_set *ls = likelihood(matches, p_quals, p_match, n);
  double p_denom = p_data(matches, p_quals, p_prior, p_match, n);
  ps->contam = (ls->contam * p_prior)/p_denom;
  ps->random = (ls->random * (1-p_prior))/p_denom;

  if (ps->contam > ps->random)
    ps->is_contam = CONTAM;
  else
    ps->is_contam = NOT_CONTAM;
  free(ls);
  return ps;  
}

likelihood_set *likelihood(const int *matches, float *p_quals, float p_match, int n) {
  likelihood_set *ls = xmalloc(sizeof(likelihood_set));
  int i;
  ls->random = 1;
  ls->contam = 1;
  for (i = 0; i < n; i++) {
    if (matches[i] == 1) {
      ls->random *= p_match;
      ls->contam *= p_quals[i];
    } else {
      ls->random *= 1 - p_match;
      ls->contam *= 1 - p_quals[i];
    }
  }
  return ls;
}

