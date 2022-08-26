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

#define IS_AU(a, b) ((a) == 'a' && (b) == 'u' || (a) == 'u' && (b) == 'a')

/*
  For simplicity in the code below, the indices of the _last_ loop stems
  are needed. See how each index maps to which sequence position:

  BRACKET  = ...((([[[....)))......]]].....
  NOTATION = ...L..R........l........r.....

  L -> left index of last stem in left loop
  R -> left index of last stem in right loop
  l -> right index of last stem in left loop
  r -> right index of last stem in right loop
*/
void skip_final_au(char *sequence, char *dot_bracket, int left_loop_stems,
                   int right_loop_stems, void (*cb)(char *, int, int)) {
  int L = -1, l = -1, R = -1, r = -1;

  char *bracket = strdup(dot_bracket);

  // find indices of last loop stems
  for (int i = 0; i < strlen(sequence); i++) {
    if (bracket[i] == '(' && L == -1) {
      L = i;
    } else if (bracket[i] == '[' && R == -1) {
      R = i;
    } else if (bracket[i] == ')') {
      l = i;
    } else if (bracket[i] == ']') {
      r = i;
    }
  }

  // left loop stem ends with an AU pair
  int left_is_au = left_loop_stems > 0 && IS_AU(sequence[L], sequence[l]);
  int right_is_au = right_loop_stems > 0 && IS_AU(sequence[R], sequence[r]);

  if (left_is_au) {
    bracket[L] = bracket[l] = '.';
    cb(bracket, left_loop_stems - 1, right_loop_stems);
    bracket[L] = '(';
    bracket[l] = ')';
  }
  if (right_is_au) {
    bracket[R] = bracket[r] = '.';
    cb(bracket, left_loop_stems, right_loop_stems - 1);
    bracket[R] = '[';
    bracket[r] = ']';
  }
  if (left_is_au && right_is_au) {
    bracket[L] = bracket[l] = bracket[R] = bracket[r] = '.';
    cb(bracket, left_loop_stems - 1, right_loop_stems - 1);
  }

  free(bracket);
}
