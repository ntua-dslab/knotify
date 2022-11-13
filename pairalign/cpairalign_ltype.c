/*
 * Copyright © 2022 Christos Pavlatos, George Rassias, Christos Andrikos,
 *                  Evangelos Makris, Aggelos Kolaitis
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the “Software”), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

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

  BRACKET  = ...(...{..[....)...}..].....
  NOTATION = ...L...M..R....l...m..r.....
             0123456789012345678901234567
             0         1         2

  L -> left index of core stem in left loop
  M -> left index of core stem in middle loop
  R -> left index of core stem in right loop
  l -> right index of core stem in left loop
  M -> right index of core stem in middle loop
  r -> right index of core stem in right loop

  a -> left index of loop stem
  b -> right index of loop stem

  For potential left loop stems, these conditions hold for a and b:
    0 <= a <= L-1         -- in example above a in [0, 2]
    l+1 <= b <= m-1       -- in example above b in [16, 18]

  For potential middle loop stems, these conditions hold for a and b:
    L+1 <= a <= M-1       -- in example above a in [4, 6]
    m+1 <= b <= r-1       -- in example above b in [20, 21]

  For potential right loop stems, these conditions hold for a and b:
    M+1 <= a <= R-1       -- in example above a in [8, 9]
    r+1 >= b <= LEN-1     -- in example above b in [23, 27]

*/
void pairalign(char *sequence, int i, int j, int left_left_loop_size,
               int left_right_loop_size, int right_left_loop_size, int dd_size,
               void (*cb)(char *, int, int, int)) {

  int L = i;
  int M = i + left_left_loop_size + 1;
  int R = i + left_left_loop_size + left_right_loop_size + 2;
  int l = i + left_left_loop_size + left_right_loop_size + dd_size + 3;
  int m = l + right_left_loop_size + 1;
  int r = i + j - 1;

  char *dot_bracket = strdup(sequence);
  int left_loop_stems = 0, middle_loop_stems = 0, right_loop_stems = 0;

  // initialize dot bracket
  int len = strlen(sequence);
  memset(dot_bracket, '.', len);
  dot_bracket[L] = '(';
  dot_bracket[l] = ')';
  dot_bracket[M] = '{';
  dot_bracket[m] = '}';
  dot_bracket[R] = '[';
  dot_bracket[r] = ']';

  // left loop stems
  left_loop_stems = 0;
  for (int a = L - 1, b = l + 1; a >= 0 && b <= m - 1; a--, b++) {
    if (!IS_PAIR(sequence[a], sequence[b])) {
      break;
    }
    dot_bracket[a] = '(';
    dot_bracket[b] = ')';
    left_loop_stems++;
  }

  middle_loop_stems = 0;
  for (int a = M - 1, b = m + 1; a >= L + 1 && b <= r - 1; a--, b++) {
    if (!IS_PAIR(sequence[a], sequence[b])) {
      break;
    }
    dot_bracket[a] = '{';
    dot_bracket[b] = '}';
    middle_loop_stems++;
  }

  // right loop stems
  right_loop_stems = 0;
  for (int a = R - 1, b = r + 1; a >= M + 1 && b <= len - 1; a--, b++) {
    if (!IS_PAIR(sequence[a], sequence[b])) {
      break;
    }
    dot_bracket[a] = '[';
    dot_bracket[b] = ']';
    right_loop_stems++;
  }

  cb(dot_bracket, left_loop_stems, middle_loop_stems, right_loop_stems);

  free(dot_bracket);
}
