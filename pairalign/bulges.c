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

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

bool gSymmetricBulges;
bool gCountStemsFromBulges;
int gMinStemsAfterBulge;
int gMaxBulgeSize;

#define MIN(a, b) ((a) < (b) ? (a) : (b))

// IS_PAIR checks for AU, GC, GU, UA, CG, UG pairs
#define IS_PAIR(a, b) (IS_PAIR_ONEWAY(a, b) || IS_PAIR_ONEWAY(b, a))

// IS_PAIR_ONEWAY checks for AU, GC, GU pairs
#define IS_PAIR_ONEWAY(a, b)                                                   \
  (((a) == 'a' && (b) == 'u') || ((a) == 'g' && ((b) == 'c' || (b) == 'u')))

/**
 * \brief given a sequence, attempt to find bulges and extra loop stems
 *
 * The left loop stems should be in the range [L, R] inclusive.
 * The right loop stems should be in the range [l, r] inclusive.
 *
 * Only a single bulge (of variable size) is allowed on each side of the loop.
 *
 * Sets the following values:
 * - o_has_bulge set to true if a bulge fitting the criteria is found
 * - o_bulge_size_left set to the size of the left bulge (starting from R-1)
 * - o_bulge_size_right set to the size of the right bulge (starting from l+1)
 * - o_stems set to the number of stems after the bulges
 *
 * \example given the following sequence and the L,R,l,r indices as shown below,
 * the result of this function would be:
 *
 * 012345678901234567890123456789012345678901
 * AUGUGCUGUAUAUGUACUUCUCUUGUUUCUCUUUUAAAAGCC
 *    I            J          i          j
 * .........(((((.................)))))......
 *
 * o_has_bulge == true
 * o_bulge_size_left == 3      (note: includes R)
 * o_bulges_size_right == 4    (note: includes l)
 * o_stems == 5
 */
void find_bulge(char *sequence, int I, int J, int i, int j, bool *o_has_bulge,
                int *o_bulge_size_left, int *o_bulge_size_right, int *o_stems) {

  *o_has_bulge = false;
  *o_stems = 0;

  for (int leftBulgeSize = 1; leftBulgeSize < MIN(J - I - 1, gMaxBulgeSize + 1);
       leftBulgeSize++) {
    for (int rightBulgeSize = 1;
         rightBulgeSize < MIN(j - i - 1, gMaxBulgeSize + 1); rightBulgeSize++) {

      if (gSymmetricBulges && leftBulgeSize != rightBulgeSize) {
        continue;
      }
      int stems = 0;

      // align stems after bulges
      for (int a = J - leftBulgeSize, b = i + rightBulgeSize;
           a >= I && b <= j && IS_PAIR(sequence[a], sequence[b]); a--, b++) {
        stems++;
      }

      // found better alignment. criteria is loop size
      if (stems >= gMinStemsAfterBulge && stems > *o_stems) {
        *o_has_bulge = true;
        *o_bulge_size_left = leftBulgeSize;
        *o_bulge_size_right = rightBulgeSize;
        *o_stems = stems;
      }
    }
  }
}

/*
  For simplicity in the code below, the indices of the _core_ loop stems
  are needed. See how each index maps to which sequence position:

  BRACKET  = ...(......[....)......].....
  NOTATION = ...L......R....l......r.....
             0123456789012345678901234567
             0         1         2

  L -> left index of core stem in left loop
  R -> left index of core stem in right loop
  l -> right index of core stem in left loop
  r -> right index of core stem in right loop

  a -> left index of loop stem
  b -> right index of loop stem

  For potential left loop stems, these conditions hold for a and b:
    0 <= a <= L-1         -- in example above a in [0, 2]
    l+1 <= b <= r-1       -- in example above b in [16, 21]

  For potential right loop stems, these conditions hold for a and b:
    L+1 <= a <= R-1       -- in example above a in [4, 9]
    r+1 >= b <= LEN-1     -- in example above b in [23, 27]

*/
void pairalign(char *sequence, int i, int j, int left_loop_size, int dd_size,
               void (*cb)(char *, int, int)) {

  int L = i;
  int R = i + left_loop_size + 1;
  int l = i + left_loop_size + dd_size + 2;
  int r = i + j - 1;

  char *dot_bracket = strdup(sequence);
  int left_loop_stems = 0, right_loop_stems = 0;

  // initialize dot bracket
  int len = strlen(sequence);
  memset(dot_bracket, '.', len);
  dot_bracket[L] = '(';
  dot_bracket[l] = ')';
  dot_bracket[R] = '[';
  dot_bracket[r] = ']';

  // left loop stems
  left_loop_stems = 0;
  for (int a = L - 1, b = l + 1; a >= 0 && b <= r - 1; a--, b++) {
    if (!IS_PAIR(sequence[a], sequence[b])) {
      break;
    }
    dot_bracket[a] = '(';
    dot_bracket[b] = ')';
    left_loop_stems++;
  }

  // right loop stems
  right_loop_stems = 0;
  for (int a = R - 1, b = r + 1; a >= L + 1 && b <= len - 1; a--, b++) {
    if (!IS_PAIR(sequence[a], sequence[b])) {
      break;
    }
    dot_bracket[a] = '[';
    dot_bracket[b] = ']';
    right_loop_stems++;
  }

  cb(dot_bracket, left_loop_stems, right_loop_stems);

  bool rBulge, lBulge;
  int rBulgeSizeRight, rBulgeSizeLeft, rStems;
  int lBulgeSizeRight, lBulgeSizeLeft, lStems;

  /*
   * 012345678901234567890123456789012345678901234567890
   * AUGUGCUGUAUAUGUACUUCUUGAAACUUGUUUCUCUUUUAAAAGCCUCGA
   * .................(((......[[)))............]]......
   * ...................L.......Rl..............r.......
   * I...............J..............i..........j........        # for left loop
   * ....................I....J...................i....j        # for right loop
   *
   * Look for bulges in left loop:
   * I = 0
   * J = L - left_loop_stems - 1
   * i = l + left_loop_stems + 1
   * j = r - 1
   *
   * Look for bulges in right loop:
   * I = L + 1
   * J = R - right_loop_stems - 1
   * i = r + right_loop_stems + 1
   * j = len(sequence) - 1
   */
  find_bulge(sequence, 0, L - left_loop_stems - 1, l + left_loop_stems + 1,
             r - 1, &lBulge, &lBulgeSizeLeft, &lBulgeSizeRight, &lStems);

  find_bulge(sequence, L + 1, R - right_loop_stems - 1,
             r + right_loop_stems + 1, len - 1, &rBulge, &rBulgeSizeLeft,
             &rBulgeSizeRight, &rStems);

  if (lBulge) {
    memset(&dot_bracket[L - left_loop_stems - lBulgeSizeLeft - lStems], '(',
           lStems);
    memset(&dot_bracket[l + left_loop_stems + 1 + lBulgeSizeRight], ')',
           lStems);

    cb(dot_bracket, left_loop_stems + (gCountStemsFromBulges ? lStems : 0),
       right_loop_stems);

    memset(&dot_bracket[L - left_loop_stems - lBulgeSizeLeft - lStems], '.',
           lStems);
    memset(&dot_bracket[l + left_loop_stems + 1 + lBulgeSizeRight], '.',
           lStems);
  }

  if (rBulge) {
    memset(&dot_bracket[R - right_loop_stems - rBulgeSizeLeft - rStems], '[',
           rStems);
    memset(&dot_bracket[r + right_loop_stems + 1 + rBulgeSizeRight], ']',
           rStems);

    cb(dot_bracket, left_loop_stems,
       right_loop_stems + (gCountStemsFromBulges ? rStems : 0));
  }

  if (lBulge && rBulge) {
    memset(&dot_bracket[L - left_loop_stems - lBulgeSizeLeft - lStems], '(',
           lStems);
    memset(&dot_bracket[l + left_loop_stems + 1 + lBulgeSizeRight], ')',
           lStems);

    cb(dot_bracket, left_loop_stems + (gCountStemsFromBulges ? lStems : 0),
       right_loop_stems + (gCountStemsFromBulges ? rStems : 0));
  }

  free(dot_bracket);
}

void initialize(int max_bulge_size, int min_stems_after_bulge,
                bool symmetric_bulges, bool count_stems_from_bulges) {
  gSymmetricBulges = symmetric_bulges;
  gMaxBulgeSize = max_bulge_size;
  gMinStemsAfterBulge = min_stems_after_bulge;
  gCountStemsFromBulges = count_stems_from_bulges;
}
