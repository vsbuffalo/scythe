/* tests.c - some tests for functions in scythe */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../scythe.h"

#define TEST_INIT int iii, ppp=0, passed=0, failed=0, checked=0;
#define TEST(test, name) if (!(test)) {                 \
    fprintf(stderr, "FAILED TEST '%s'\n", name);        \
    failed++;                                         \
  } else {                                            \
    printf("PASSED TEST '%s'\n", name);                 \
    passed++;                                         \
  } checked++;                                        \

#define TEST_ARRAY(array1, array2, n, name) for (iii = 0; iii < n; iii++) { \
    ppp = 0;                                                            \
    if (array1[iii] != array2[iii]) {                                   \
      fprintf(stderr, "FAILED TEST '%s'\n", name);                      \
      ppp = 1; failed++;                                                \
      break;                                                            \
    }                                                                   \
  }                                                                     \
  if (ppp == 0) {                                                       \
    passed++;                                                           \
    printf("PASSED TEST '%s'\n", name);                                 \
  }                                                                     \
  checked++;                                                            \

#define TEST_CLOSE printf("\ntests: %d\npassed: %d\nfailed: %d\n", checked, passed, failed);
 
void test_score_sequence(void) {
  char *seqa = "ATCGATCGATCGATCGATCGATCG";
  char *seqb = "AACGATCGATCGATCGATCGATCG";
  int cmp1[] = {1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  int *cmp1_t = calloc(24, sizeof(*cmp1_t));
  score_sequence(cmp1_t, seqa, seqb, strlen(seqa));

  TEST_INIT;
  TEST_ARRAY(cmp1, cmp1_t, 24, "score_sequence");

  TEST(sum(score_sequence(cmp1_t, seqa, seqb, strlen(seqa)), strlen(seqa)) == 22, "sum of score_sequence");

  free(cmp1_t);
  TEST_CLOSE;
}

int main(int argc, char *argv[]) {
  test_score_sequence();
  return 0;
}
