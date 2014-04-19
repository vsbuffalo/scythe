/*
  matches.c - sequence match scoring for scythe.
*/
#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "scythe.h"

int *score_sequence(int *matches, const char *seq, const char *pattern, size_t n) {
  /* 
     Score a string using constant match and mismatch scores. Assumes
     seq is longer than or of equal length to pattern. Only the first
     n characters of pattern and seq are scored. We explicitly error
     out of in the case of a null byte.


     [----------seq-------------------]
     |----pattern-----|

  */
  size_t i;
/* These are addressed by the assert below
  assert(strlen(seq) >= n);
  assert(strlen(pattern) >= n);
 */
  for (i = 0; i < n; i++) {
    assert(seq[i] && pattern[i]); /* no string termination */
    if (seq[i] == pattern[i])
      matches[i] = MATCH_SCORE;
    else
      matches[i] = MISMATCH_SCORE;
  }
  return matches;
}

match *find_best_match(const adapter_array *aa, const char *read,
                       float *p_quals, float prior, float p_match, int min_l) {
  /*
   Take an adapter array, and check the read against all
   adapters. Brute force string matching is used. This is to avoid
   approximate matching algorithms which required an a priori
   specified number mismatches.

  */

  match *best_match=NULL;
  int i, shift, max_shift, found_contam=0;
  size_t al, rl = strlen(read);
  int *best_arr=xmalloc(sizeof(int) * rl), best_adapter=0, best_length=0, best_shift=0, best_score=INT_MIN;
  int curr_score, *curr_arr=xmalloc(sizeof(int) * rl);
  posterior_set *ps=NULL;
  float *best_p_quals=NULL;

  max_shift = rl - min_l;
  for (shift = 0; shift < max_shift; shift++) {
    for (i = 0; i < aa->n; i++) {
      if (min_l >= aa->adapters[i].length) {
        fprintf(stderr, "Minimum match length (option -n) greater than or " \
                "equal to length of adapter.\n");
        exit(EXIT_FAILURE);
      }
      al = min(aa->adapters[i].length, rl - shift);
      score_sequence(curr_arr, &(read)[shift], (aa->adapters[i]).seq, al);
      curr_score = sum(curr_arr, al);
      if (curr_score > best_score) {
        best_score = curr_score;
        best_length = al;
        best_shift = shift;
        best_p_quals = &(p_quals)[shift];
        for (size_t iii = 0; iii < al; iii++)
            best_arr[iii] = curr_arr[iii];
        best_adapter = i;
        ps = posterior(best_arr, best_p_quals, prior, 0.25, best_length);
        found_contam = ps->is_contam;
        if (found_contam) {
          break;
        } else {
          free(ps);
          ps=NULL;
        }
      }
    }
    if (found_contam)
      break;
  }
  free(curr_arr);
  curr_arr = NULL;

  if (!found_contam) {
    /* no match found */
    free(best_arr);
    best_arr = NULL;
    return NULL;
  }
  best_arr = realloc(best_arr, sizeof(*best_arr) * best_length);
  /* save this match */
  best_match = xmalloc(sizeof(match));
  best_match->match = best_arr;
  best_match->shift = best_shift;
  best_match->length = best_length;
  best_match->ps = ps;
  best_match->score = best_score;
  best_match->adapter_index = best_adapter;
  best_match->p_quals = best_p_quals;
  best_match->match_pos = calloc(best_length, sizeof(int));
  for (i = 0; i < best_length; i++)
    best_match->match_pos[i] = best_match->match[i] == MATCH_SCORE;
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

