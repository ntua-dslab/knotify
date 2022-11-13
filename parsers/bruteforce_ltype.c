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

// Description: Brute force library for detecting positions of core stems
// of pseudoknot in an RNA sequence.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct stem {
  unsigned char left, right; // [0..255]
} stem;

typedef struct core_stem {
  stem cstem[3];
} core_stem;

static int s_allow_ug;
static int s_min_dd_size;
static int s_max_dd_size;
static int s_min_window_size;
static int s_max_window_size;
static float s_min_window_size_ratio;
static float s_max_window_size_ratio;

void initialize(char *_grammar_unused, int _allow_ug, int _min_dd_size,
                int _max_dd_size, int _min_window_size, int _max_window_size,
                float _min_window_size_ratio, float _max_window_size_ratio) {
  s_allow_ug = _allow_ug;
  s_min_dd_size = _min_dd_size;
  s_max_dd_size = _max_dd_size;
  s_min_window_size = _min_window_size;
  s_max_window_size = _max_window_size;

  s_min_window_size_ratio = _min_window_size_ratio;
  s_max_window_size_ratio = _max_window_size_ratio;

  // printf("%d/%d/%d/%d/%d\n", allow_ug, min_dd_size, max_dd_size,
  //        min_window_size, max_window_size);
}

void detect_pseudoknots(char *sequence,
                        void (*cb)(int, int, int, int, int, int)) {
  int n = strlen(sequence);

  // window size is static or ratio of sequence length
  int min_window_size =
      s_min_window_size ? s_min_window_size : (n * s_min_window_size_ratio);
  int max_window_size =
      s_max_window_size ? s_max_window_size : (n * s_max_window_size_ratio);

  stem *all_stems = (stem *)malloc(n * n * sizeof(stem));
  stem *cs_position = all_stems;
  stem *cs_position2 = all_stems;
  stem *cs_position3 = all_stems;

  long int n_stems = 0;
  long int n_cstems = 0;

  for (int i = 0; i < n; i++) {
    switch (sequence[i]) {
    case 'a':
      for (int j = i + 1; j < n; j++)
        if (sequence[j] == 'u') {
          cs_position->left = i;
          cs_position->right = j;
          cs_position++;
          n_stems++;
        }
      break;
    case 'u':
      for (int j = i + 1; j < n; j++)
        if ((sequence[j] == 'a') || ((sequence[j] == 'g') && s_allow_ug)) {
          cs_position->left = i;
          cs_position->right = j;
          cs_position++;
          n_stems++;
        }
      break;
    case 'c':
      for (int j = i + 1; j < n; j++)
        if (sequence[j] == 'g') {
          cs_position->left = i;
          cs_position->right = j;
          cs_position++;
          n_stems++;
        }
      break;
    case 'g':
      for (int j = i + 1; j < n; j++)
        if ((sequence[j] == 'c') || ((sequence[j] == 'u') && s_allow_ug)) {
          cs_position->left = i;
          cs_position->right = j;
          cs_position++;
          n_stems++;
        }
      break;
    default:
      printf("Invalid character\n");
    }
  }

  // printf("Stems (size: %ld): \n", n_stems);
  // cs_position = all_stems;
  for (int i = 0; i < n_stems; i++) {
    cs_position = all_stems + i;
    // printf("%d) \t%d %d\n", i, cs_position->left, cs_position->right);
  }

  // printf("eftasa edo\n");
  core_stem *all_cstems =
      malloc(n_stems * n_stems * n_stems * sizeof(core_stem)); // zong
  core_stem *ccs_position = all_cstems;

  cs_position = all_stems;
  cs_position2 = all_stems;
  cs_position3 = all_stems;
  ccs_position = all_cstems;

  /*for (int i = 0; i < n_stems; i++) {
    cs_position = all_stems + i;
    for (int j = i + 1; j < n_stems; j++) {
      cs_position2 = all_stems + j;
      if ((cs_position->left < cs_position2->left) &&
          (cs_position2->left < cs_position->right) &&
          (cs_position->right < cs_position2->right)) {
        // printf("petixe %ld) %d %d\n",n_cstems,i,j);
        (ccs_position->cstem[0]).left = cs_position->left;
        (ccs_position->cstem[0]).right = cs_position->right;
        (ccs_position->cstem[1]).left = cs_position2->left;
        (ccs_position->cstem[1]).right = cs_position2->right;

        ccs_position++;
        n_cstems++;
      }
    }
  }*/

  for (int i = 0; i < n_stems; i++) {
    cs_position = all_stems + i;
    for (int j = i + 1; j < n_stems; j++) {
      cs_position2 = all_stems + j;
      for (int k = j + 1; k < n_stems; k++) {
        cs_position3 = all_stems + k;
        if ((cs_position->left < cs_position2->left) &&
            (cs_position2->left < cs_position3->left) &&
            (cs_position3->left < cs_position->right) &&
            (cs_position->right < cs_position2->right) &&
            (cs_position2->right < cs_position3->right)) {
          ccs_position->cstem[0].left = cs_position->left;
          ccs_position->cstem[0].right = cs_position->right;
          ccs_position->cstem[1].left = cs_position2->left;
          ccs_position->cstem[1].right = cs_position2->right;
          ccs_position->cstem[2].left = cs_position3->left;
          ccs_position->cstem[2].right = cs_position3->right;
          ccs_position++;
          n_cstems++;
        }
      }
    }
  }

  // printf("Core_Stems (size: %ld): \n", n_cstems);
  ccs_position = all_cstems;
  for (int i = 0; i < n_cstems; i++) {
    int left = (ccs_position->cstem[0]).left;
    int size = (ccs_position->cstem[2]).right - ccs_position->cstem[0].left + 1;
    int left_left_loop_size =
        ccs_position->cstem[1].left - ccs_position->cstem[0].left - 1;
    int left_right_loop_size =
        ccs_position->cstem[2].left - ccs_position->cstem[1].left - 1;
    int right_left_loop_size =
        ccs_position->cstem[1].right - ccs_position->cstem[0].right - 1;
    int dd_size =
        ccs_position->cstem[0].right - ccs_position->cstem[2].left - 1;

    // TODO (akolaitis): improve this
    if (size < min_window_size || size > max_window_size ||
        dd_size < s_min_dd_size || dd_size > s_max_dd_size ||
        left_left_loop_size == 0 || left_right_loop_size == 0 ||
        right_left_loop_size == 0 ||
        (ccs_position->cstem[1].right == ccs_position->cstem[2].right - 1)) {
      // printf("RIP! %d,%d,%d,%d\n", left, size, left_loop_size, dd_size);

      ccs_position++;
      continue;
    }

    cb(left, size, left_left_loop_size, left_right_loop_size,
       right_left_loop_size, dd_size);
    // printf("(%d %d) - (%d %d) \t\t (%c, %c) - (%c, %c) \n",
    //        (ccs_position->cstem[0]).left, (ccs_position->cstem[0]).right,
    //        (ccs_position->cstem[1]).left, (ccs_position->cstem[1]).right,
    //        sequence[(ccs_position->cstem[0]).left],
    //        sequence[(ccs_position->cstem[0]).right],
    //        sequence[(ccs_position->cstem[1]).left],
    //        sequence[(ccs_position->cstem[1]).right]);
    ccs_position++;
  }
}
