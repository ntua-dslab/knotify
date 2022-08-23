#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))

// IS_PAIR checks for AU, GC, GU, UA, CG, UG pairs
#define IS_PAIR(a, b) (IS_PAIR_ONEWAY(a, b) || IS_PAIR_ONEWAY(b, a))

// IS_PAIR_ONEWAY checks for AU, GC, GU pairs
#define IS_PAIR_ONEWAY(a, b)                                                   \
  (((a) == 'a' && (b) == 'u') || ((a) == 'g' && ((b) == 'c' || (b) == 'u')))

/*
  For simplicity in the code below, the indices of the _core_ loop stems
  are needed. See how each index maps to which sequence position:

  BRACKET  = ...(.[....)......].....
  NOTATION = ...L.R....l......r.....
             01234567890123456789012
             0         1         2

  L -> left index of core stem in left loop
  R -> left index of core stem in right loop
  l -> right index of core stem in left loop
  r -> right index of core stem in right loop

  Potential left loop stems [0, L) and (l, r)    --- example [0, 2] and [11, 16]
  Potential right loop stems (L, R) and (r, len) --- example [4, 4] and [18, 22]

  Writes results to oDotBracket, oRightLoopStems, oLeftLoopStems
*/
void pairalign(char *sequence, int i, int j, int left_loop_size, int dd_size,
               char *oDotBracket, int *oLeftLoopStems, int *oRightLoopStems) {

  int L = i;
  int R = i + left_loop_size + 1;
  int l = i + left_loop_size + dd_size + 2;
  int r = i + j - 1;

  // initialize dot bracket
  int len = strlen(sequence);
  memset(oDotBracket, '.', len);
  oDotBracket[L] = '(';
  oDotBracket[l] = ')';
  oDotBracket[R] = '[';
  oDotBracket[r] = ']';

  // left loop stems
  *oLeftLoopStems = 0;
  for (int i = 1; i <= MIN(L, r - l - 1); i++) {
    if (!IS_PAIR(sequence[L - i], sequence[l + i])) {
      break;
    }
    oDotBracket[L - i] = '(';
    oDotBracket[l + i] = ')';
    (*oLeftLoopStems)++;
  }

  // right loop stems
  *oRightLoopStems = 0;
  for (int i = 1; i <= MIN(R - L - 1, len - r); i++) {
    if (!IS_PAIR(sequence[R - i], sequence[r + i])) {
      break;
    }
    oDotBracket[R - i] = '[';
    oDotBracket[r + i] = ']';
    (*oRightLoopStems)++;
  }
}
