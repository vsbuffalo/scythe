/*
  matches.c - sequence match scoring for scythe.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "scythe.h"


int *score_sequence(const char *seqa, const char *seqb, int n) {
  int i;
  int *matches;
  /* printf("find_matches(): %s\n%s\n", seqa, seqb); */
  /* if (strlen(seqa) != strlen(seqb)) */
  /*   die("matches() comparing strings of different lengths; unexpected error."); */
  matches = xmalloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    if (seqa[i] == seqb[i])
      matches[i] = MATCH_SCORE;
    else
      matches[i] = MISMATCH_SCORE;
  }
  return matches;
}

match *find_best_match(const adapter_array *aa, const char *read, 
                       float *p_quals, float prior, float p_match, int min) {
  /*
    Take an adapter_array, and check the read against all adapters,
    minding 3' and 5' ends.

    Returns NULL on no match with score > 0.
  */
  match *best_match;
  int i, l, first=1, rl=strlen(read), *max=NULL, *m=NULL, max_score=0, best_n=0, current_score, best_adapter;
  float *best_p_quals; /* for the subset of read qualities of the matching sequence */
  for (i = 0; i < aa->n; i++) {
    if (min >= aa->adapters[i].length) {
      fprintf(stderr, "error: Minimum match length (option -n) greater than or equal to length of adapter.\n");
      exit(EXIT_FAILURE);
    }
    for (l = (aa->adapters[i]).length; l > min; l--) {
      m = score_sequence(&(read)[rl-l], (aa->adapters[i]).seq, l);
      current_score = sum(m, l);
      
      /* the first sequence comparison is always the max_score */ 
      if (first) {
        max = m;
        max_score = current_score;
        best_p_quals = &(p_quals)[rl-l];
        best_n = l;
        first = 0;
        best_adapter = i;
        continue;
      }
      
      if (current_score > max_score) {
        if (max)
          free(max); /* free last max array */
        else
          max = xmalloc(l*sizeof(int));
        max = m;
        max_score = current_score;
        best_p_quals = &(p_quals)[rl-l];
        best_n = l;
        best_adapter = i;
      } else {
        free(m);
      }
      
      /* early exit when it's no longer possible to score higher */
      if (max_score >= l-1)
        break;
    }
  }

  if (max_score <= 0) {
    free(max);
    return NULL;
  }
  /* save this match */ 
  best_match = xmalloc(sizeof(match));
  best_match->match = max;
  best_match->n = best_n;
  best_match->score = max_score;
  best_match->adapter_index = best_adapter;
  best_match->p_quals = best_p_quals;
  best_match->match_pos = xmalloc(best_n*sizeof(int));
  for (i = 0; i < best_n; i++) {
    if (best_match->match[i] == MATCH_SCORE)
      best_match->match_pos[i] = 1;
    else 
      best_match->match_pos[i] = 0; 
  }

  /* calculate posterior qualities */
  best_match->ps = posterior(best_match->match_pos, best_p_quals, prior, 0.25, best_match->n);
  return best_match;
}

void destroy_match(match *m) {
  /* 
     Free all memory taken by this match object. m->p_quals is not
     freed, as this is points to an array allocated within main().
  */
  free(m->match);
  free(m->match_pos);
  free(m->ps);
  free(m);
}
