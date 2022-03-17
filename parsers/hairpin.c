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

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hashtab.h"
#include "objstack.h"
#include "yaep.h"

#define TRUE 1
#define FALSE 0

static char *input;
static int ntok;

static int read_token_func(void **attr) {
  *attr = NULL;
  if (input[ntok]) {
    return input[ntok++];
  } else {
    return -1;
  }
}

static void syntax_error_func(int err_tok_num, void *err_tok_attr,
                              int start_ignored_tok_num,
                              void *start_ignored_tok_attr,
                              int start_recovered_tok_num,
                              void *start_recovered_tok_attr) {
  return; // performance trick
  if (start_ignored_tok_num < 0) {
    fprintf(stderr, "Syntax error on token %d\n", err_tok_num);
  } else {
    fprintf(
        stderr,
        "Syntax error on token %d:ignore %d tokens starting with token = %d\n",
        err_tok_num, start_recovered_tok_num - start_ignored_tok_num,
        start_ignored_tok_num);
  }
}

static void *parse_alloc_func(int size) { return malloc(size); }

struct yaep_tree_node *parse(const char *description) {
  struct grammar *g;
  struct yaep_tree_node *root;
  int ambiguous_p;

  if ((g = yaep_create_grammar()) == NULL) {
    fprintf(stderr, "yaep_create_grammar: No memory\n");
    exit(1);
  }

  // setup some configuration flags
  yaep_set_lookahead_level(g, 2);
  yaep_set_debug_level(g, 0);
  yaep_set_one_parse_flag(g, 0);
  yaep_set_cost_flag(g, 0);

  if (yaep_parse_grammar(g, TRUE, description) != 0) {
    fprintf(stderr, "%s\n", yaep_error_message(g));
    exit(1);
  }

  int parsed = yaep_parse(g, read_token_func, syntax_error_func,
                          parse_alloc_func, NULL, &root, &ambiguous_p);
  if (parsed) {
    printf("this should not be accepted\n");
    fprintf(stderr, "yaep_parse: %s\n", yaep_error_message(g));
  }

  return root;
}

// Helper structs and function follow
struct hairpin {
  int start;
  int stems;
  int size;
  struct hairpin *next;
};

struct hairpin *initialize_new_hairpin_node() {
  struct hairpin *new_hairpin =
      (struct hairpin *)malloc(sizeof(struct hairpin));
  new_hairpin->start = 0;
  new_hairpin->stems = 0;
  new_hairpin->size = 0;
  new_hairpin->next = NULL;
  return new_hairpin;
}

struct hairpin *append_to_hairpin_list(struct hairpin *l1, struct hairpin *l2) {
  if (l1 == NULL) {
    return l2;
  }
  struct hairpin *iter = l1, *last;
  for (iter = l1; iter != NULL; last = iter, iter = iter->next) {
  }
  last->next = l2;
  return l1;
}

/**
 * Returns a hairpin list with the corresponding values, i.e. `stems` and
 * `size`, *equal or above* the given thresholds
 */
struct hairpin *filter_hairpin_list(struct hairpin *list, int stems, int size) {
  for (struct hairpin *h = list, *prev = NULL; h != NULL;) {
    if (h->stems < stems || h->size < size) {
      if (prev == NULL) {
        prev = h;
        h = h->next;
        free(prev);
        prev = NULL;
        list = h;
      } else {
        prev->next = h->next;
        free(h);
        h = prev->next;
      }
    } else {
      if (prev == NULL) {
        prev = h;
      } else {
        prev = prev->next;
      }
      h = h->next;
    }
  }
  return list;
}

void parse_left_subtree(struct yaep_tree_node *node, struct hairpin *hp) {
  if (node->type == YAEP_NIL) {
    return;
  }
  if (strstr(node->val.anode.name, "K") != NULL) {
    hp->start++;
  }
  parse_left_subtree(node->val.anode.children[0], hp);
}

/**
 * This (recursive) function should be called from `parse_tree` with the
 * *middle* child of an `S1` YAEP_ANODE. That child can either contain
 * a single tree (YAEP_ANODE) or multiple (YAEP_ALT).
 *
 * The purpose of this function is to find the size of the inner hairpin
 * loop, so it also leverages the `struct hairpin`, but make use only of
 * the `size` member.
 */
struct hairpin *parse_middle_subtree(struct yaep_tree_node *node,
                                     struct hairpin *hp, int start) {
  if (node == NULL) {
    return hp;
  }
  switch (node->type) {
  case YAEP_ANODE: {
    if (hp == NULL) {
      // first call - single tree
      hp = initialize_new_hairpin_node();
      hp->start = start;
    }
    if (strstr(node->val.anode.name, "K") != NULL) {
      hp->size++;
    }
    // we continue with a recursive call
    return parse_middle_subtree(node->val.anode.children[0], hp, start);
  }
  case YAEP_ALT: {
    struct yaep_tree_node *alt_child_node = node->val.alt.node;
    struct hairpin *h0 = NULL;
    int when = ((alt_child_node->type == YAEP_ANODE)
                // && (strcmp(alt_child_node->val.anode.name, "P5") != 0)
    );
    if (when) {
      // we don't want to consider alternative cases that do not start from
      // one of {`P1`, `P2`, `P3`, `P4`} rules
      struct hairpin *new0 = initialize_new_hairpin_node();
      new0->start = start;
      if (hp != NULL) {
        new0->size = hp->size;
      }
      h0 = parse_middle_subtree(node->val.alt.node, new0, start);
    }

    struct hairpin *new1 = initialize_new_hairpin_node();
    new1->start = start;
    if (node->val.alt.next != NULL) {
      struct hairpin *h1 =
          parse_middle_subtree(node->val.alt.next, new1, start);
      h0 = append_to_hairpin_list(h0, h1);
    }
    return h0;
  }
  case YAEP_NIL:
  default:
    return hp;
  }
}

int get_outer_subtree_length(struct yaep_tree_node *node) {
  if (node->type == YAEP_NIL) {
    return 0;
  }
  if (strstr(node->val.anode.name, "K") != NULL) {
    return 1 + get_outer_subtree_length(node->val.anode.children[0]);
  }
  return get_outer_subtree_length(node->val.anode.children[0]);
}

void calculate_stems(struct hairpin *node, int right_outer_length,
                     int seq_length) {
  while (node != NULL) {
    node->stems =
        (seq_length - (node->start + node->size + right_outer_length)) / 2;
    node = node->next;
  }
}

/**
 * Receives a top-level `YAEP_ANODE` (`S1`) tree representation and
 * calculates a list of hairpin structures. The list can contain a
 * single hairpin or multiple, depending on P ambigiuties of the `S1`
 * rule.
 */
struct hairpin *parse_full_tree(struct yaep_tree_node *node, int seq_length) {
  if (strcmp(node->val.anode.name, "S1") != 0) {
    // something's wrong
    return NULL;
  } else {
    // L
    int left_outer_length =
        get_outer_subtree_length(node->val.anode.children[0]);
    // P
    struct hairpin *hairpins = parse_middle_subtree(node->val.anode.children[1],
                                                    NULL, left_outer_length);
    // R
    int right_outer_length =
        get_outer_subtree_length(node->val.anode.children[2]);
    calculate_stems(hairpins, right_outer_length, seq_length);
    return hairpins;
  }
}

/**
 * Receives the "root" of the DAG representation produced by yaep
 * containing all the different parsed trees. Basically, there can
 * be two cases for the root node:
 * * `YAEP_ANODE` (`S1`) for a single tree
 * * `YAEP_ALT` for an arbitrary list of different trees (`YAEP_ANODE`)
 * Calls `parse_full_tree` with a `S1` `YAEP_ANODE`.
 */
struct hairpin *traverse_yaep_solution(struct yaep_tree_node *node,
                                       int seq_length) {
  if (node == NULL) {
    return NULL;
  }
  struct hairpin *hairpins = NULL;
  switch (node->type) {
  case YAEP_ANODE: {
    hairpins = parse_full_tree(node, seq_length);
  }
  case YAEP_ALT: {
    struct yaep_tree_node *alt;
    struct hairpin *alt_hairpins = NULL;
    for (alt = node; alt != NULL; alt = alt->val.alt.next) {
      alt_hairpins = parse_full_tree(alt->val.alt.node, seq_length);
      if (hairpins == NULL) {
        hairpins = alt_hairpins;
        continue;
      }
      hairpins = append_to_hairpin_list(hairpins, alt_hairpins);
    }
  }
  default:
    break;
  }
  return hairpins;
}

int hairpin_equal(struct hairpin *h1, struct hairpin *h2) {
  if (h1 == NULL || h2 == NULL) {
    return FALSE;
  }
  if (h1->start == h2->start && h1->stems == h2->stems &&
      h1->size == h2->size) {
    return TRUE;
  } else {
    return FALSE;
  }
}

struct hairpin *get_lhairpin(struct hairpin *head, int l) {
  for (; head != NULL && l > 0; head = head->next, l--) {
  }
  return (head != NULL) ? head : NULL;
}

/**
 * Swaps the i-th and j-th hairpin node in hairpin node list
 */
void swap(struct hairpin *head, int i, int j) {
  struct hairpin *tmp, *i_tmp, *j_tmp;
  int iter;

  for (tmp = head, iter = i; tmp != NULL && iter > 0; tmp = tmp->next, iter--) {
  }
  i_tmp = tmp;

  for (tmp = head, iter = j; tmp != NULL && iter > 0; tmp = tmp->next, iter--) {
  }
  j_tmp = tmp;

  int tmp_size = i_tmp->size;
  int tmp_stems = i_tmp->stems;
  int tmp_start = i_tmp->start;

  i_tmp->size = j_tmp->size;
  i_tmp->stems = j_tmp->stems;
  i_tmp->start = j_tmp->start;

  j_tmp->size = tmp_size;
  j_tmp->stems = tmp_stems;
  j_tmp->start = tmp_start;
}

int compare(struct hairpin *a, struct hairpin *b) {
  int comp = a->stems - b->stems;
  if (comp < 0) {
    return -1;
  }
  if (comp > 0) {
    return 1;
  }
  comp = a->size - b->size;
  if (comp < 0) {
    return -1;
  }
  if (comp > 0) {
    return 1;
  }
  comp = a->start - b->start;
  if (comp < 0) {
    return -1;
  }
  if (comp > 0) {
    return 1;
  }
  return comp;
}

/**
 * Sorts a hairpin list based on (with order of priority):
 * * stems
 * * size
 * * start
 */
struct hairpin *quicksort_hairpin_list(struct hairpin *head, int l, int r) {
  int i, j;
  struct hairpin *pivot, *jhairpin;
  i = l + 1;
  if (l + 1 < r) {
    pivot = get_lhairpin(head, l);
    for (j = l + 1; j <= r; j++) {
      jhairpin = get_lhairpin(head, j);
      if (jhairpin != NULL && compare(jhairpin, pivot) > 0) {
        swap(head, i, j);
        i++;
      }
    }
    swap(head, i - 1, l);
    quicksort_hairpin_list(head, l, i);
    quicksort_hairpin_list(head, i, r);
  }
  return head;
}

struct hairpin *sort_hairpin_list(struct hairpin *head) {
  struct hairpin *tmp;
  int n;
  for (tmp = head, n = 0; tmp != NULL; tmp = tmp->next, n++) {
  }
  head = quicksort_hairpin_list(head, 0, n);
  return head;
}

char *hairpin_to_string(struct hairpin *h, int str_length) {
  char *str = (char *)malloc((str_length + 1) * sizeof(char));
  int i, ind;
  for (i = 0, ind = 0; i < h->start; i++) {
    str[i] = '.';
  }
  for (ind += i, i = 0; i < h->stems; i++) {
    str[i + ind] = '(';
  }
  for (ind += i, i = 0; i < h->size; i++) {
    str[i + ind] = '.';
  }
  for (ind += i, i = 0; i < h->stems; i++) {
    str[i + ind] = ')';
  }
  for (ind += i; ind < str_length; ind++) {
    str[ind] = '.';
  }
  str[ind] = '\0';
  return str;
}

char *initialize_bracket(int loop_length) {
  char *str = malloc((loop_length + 1) * sizeof(char));
  memset(str, '.', loop_length);
  str[loop_length] = '\0';
  return str;
}

void annotate_hairpin(char *buf, struct hairpin *h) {
  memset(buf + h->start, '(', h->stems);
  memset(buf + h->start + h->stems + h->size, ')', h->stems);
}

// Return TRUE if [x1, x2] & [y1, y2] = empty.
int empty_subset(int x1, int x2, int y1, int y2) {
  return !(y1 <= x1 && x1 <= y2) && !(y1 <= x2 && x2 <= y2);
}

// Macros for hairpin indices
// 012345678901234567890
// ....((((......))))...
// LSTART=4, LEND=7, RSTART=14, REND=17
#define LSTART(h) ((h)->start)
#define LEND(h) ((h)->start + (h)->stems - 1)
#define RSTART(h) ((h)->start + (h)->stems + (h)->size)
#define REND(h) (RSTART((h)) + (h)->stems - 1)

// Check whether
//      012345678901234567890
// A == ......(((....))).....
// B == ..(((...........)))..
int hairpin_check_bulge(struct hairpin *A, struct hairpin *B, int max_bulge) {
  int left = LSTART(A) - LEND(B) - 1;
  int right = RSTART(B) - REND(A) - 1;

  // reject case where left==right==0, because this is a hairpin
  // that the grammar has already found
  return left >= 0 && right >= 0 && (left <= max_bulge && right <= max_bulge);
}

// Check whether
//      012345678901234567890
// A == ...........(((....)))
// B == ..(((..)))...........
int hairpin_check_no_overlap(struct hairpin *A, struct hairpin *B) {
  return REND(B) < LSTART(A);
}

// TODO(akolaitis): use recursion to support arbitrary max_per_loop size.
void write_hairpin_trees(struct hairpin *current, struct hairpin *next,
                         FILE *fp, int loop_length, int max_per_loop,
                         int max_bulge, int min_stems) {
  struct hairpin *iter = current;

  iter = current;
  char *buf = initialize_bracket(loop_length);
  while (iter != NULL) {
    // add to results
    if (iter->stems >= min_stems) {
      memset(buf, '.', loop_length);
      annotate_hairpin(buf, iter);
      fprintf(fp, "%s\n", buf);
    }

    if (max_per_loop > 1 || max_bulge > 0) {
      for (struct hairpin *nest = next; nest != NULL; nest = nest->next) {
        if (iter->stems + nest->stems < min_stems) {
          continue;
        }

        int merge = max_per_loop > 1 && (hairpin_check_no_overlap(iter, nest) ||
                                         hairpin_check_no_overlap(nest, iter));

        if (!merge && max_bulge > 0) {
          merge = hairpin_check_bulge(iter, nest, max_bulge) ||
                  hairpin_check_bulge(nest, iter, max_bulge);
        }

        if (!merge) {
          continue;
        }

        memset(buf, '.', loop_length);
        annotate_hairpin(buf, iter);
        annotate_hairpin(buf, nest);
        fprintf(fp, "%s\n", buf);
      }
    }

    iter = iter->next;
  }
}

char *detect_hairpins(char *grammar, char *sequence, int min_stems,
                      int min_size, int max_per_loop, int max_bulge) {
  ntok = 0;
  input = sequence;
  int input_len = strlen(input);

  char *buffer;
  size_t size;
  FILE *fp = open_memstream(&buffer, &size);

  struct yaep_tree_node *root = parse(grammar);
  struct hairpin *list = traverse_yaep_solution(root, input_len);

  list = filter_hairpin_list(list, 1, min_size);
  list = sort_hairpin_list(list);

  // drop duplicates from sorted list
  struct hairpin *iter = list, *prev = NULL;
  while (iter != NULL) {
    struct hairpin *next = iter->next;
    if (hairpin_equal(iter, prev)) {
      prev->next = next;
      free(iter);
    } else {
      prev = iter;
    }
    iter = next;
  }

  // recursively print hairpin trees with up to max_per_loop hairpins for the
  // same loop
  write_hairpin_trees(list, list ? list->next : NULL, fp, input_len,
                      max_per_loop, max_bulge, min_stems);

  fclose(fp);
  return buffer;
}
