// standard libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// custom imports
#include "hashtab.h"
#include "objstack.h"
#include "yaep.h"

// include any libs for testing ops
#include <assert.h>

#define TRUE 1
#define FALSE 0

static char *input;
static int ntok;
static int max_dd_size;

/*************************************************************************
 *                      GRAMMAR DEFINITION                             *
 *************************************************************************/

/* Grammar is defined in Python code and passed as an argument */

/*************************************************************************
                Struct definitions:
**************************************************************************/

struct pseudoknot {
  int left_loop_size;
  int dd_size;
  struct pseudoknot *next;
};

// struct to keep a long list of integers
struct size {
  int value;
  struct size *next;
};

/*************************************************************************
 *                Struct management:                                      *
 **************************************************************************/

struct pseudoknot *create_pseudoknot(int left_loop_size, int dd_size) {
  struct pseudoknot *ps =
      (struct pseudoknot *)malloc(sizeof(struct pseudoknot));
  ps->left_loop_size = left_loop_size;
  ps->dd_size = dd_size;
  ps->next = NULL;
  return ps;
};

struct size *create_size(int size) {
  struct size *size_node = (struct size *)malloc(sizeof(struct size));
  size_node->value = size;
  size_node->next = NULL;
  return size_node;
};

struct pseudoknot *
append_pseudoknot_if_not_exists(struct pseudoknot *list,
                                struct pseudoknot *upcoming_knot) {
  // FIXME: this is more like an assertion and may be removed
  // if (upcoming_knot->left_loop_size<1){
  //     return list;
  // }

  // TODO(akolaitis): This currently iterates over the whole list when
  // inserting a new pseudoknot, making insert an O(n) operation, and
  // the list creation an O(n^2) in total. We can do much better than
  // that.
  for (struct pseudoknot *it = list; it != NULL; it = it->next) {
    if (it->left_loop_size == upcoming_knot->left_loop_size &&
        it->dd_size == upcoming_knot->dd_size) {
      // printf(" -- skipping size: %d, upcoming dd %d \n",
      // upcoming_knot->left_loop_size, upcoming_knot->dd_size);
      return list;
    }
  }
  // printf(" -- Upcoming size: %d, upcoming dd %d \n",
  // upcoming_knot->left_loop_size, upcoming_knot->dd_size);
  upcoming_knot->next = list;
  return upcoming_knot;
}

struct size *append_size_if_not_exists(struct size *list, int size) {
  struct size *size_node;
  struct size *iter = list;

  while (iter != NULL) {
    if (iter->value == size) {
      return list;
    }
    iter = iter->next;
  }
  size_node = create_size(size);
  size_node->next = list;
  return size_node;
}

/**
 * concatenate two lists into one, this is a cut
 * meaning that elements are going to be maintained uniquely
 **/
struct pseudoknot *concatenate_punique(struct pseudoknot *list1,
                                       struct pseudoknot *list2) {
  struct pseudoknot *iter2 = list2, *next = NULL;
  struct pseudoknot *new_list_head = list1;

  if (list1 == NULL) {
    return list2;
  }
  if (list2 == NULL) {
    return list1;
  }

  while (iter2 != NULL) {
    next = iter2->next;
    new_list_head = append_pseudoknot_if_not_exists(new_list_head, iter2);
    iter2 = next;
  }

  return new_list_head;
}

struct size *concatenate_sunique(struct size *list1, struct size *list2) {
  struct size *iter2 = list2;
  struct size *new_list_head = list1;

  if (list1 == NULL) {
    return list2;
  }
  if (list2 == NULL) {
    return list1;
  }
  while (iter2 != NULL) {
    new_list_head = append_size_if_not_exists(new_list_head, iter2->value);

    iter2 = iter2->next;
  }

  return new_list_head;
}

/************************************************************************
 *                           Helper Functions                            *
 *************************************************************************/
void print_node_type(struct yaep_tree_node *node) {
  switch (node->type) {
  case YAEP_ANODE:
    printf("YAEP_NODE\n");
    break;
  case YAEP_NIL:
    printf("YAEP_NIL\n");
    break;
  case YAEP_ERROR:
    printf("YAEP_ERROR\n");
    break;
  case YAEP_TERM:
    printf("YAEP_TERM\n");
    break;
  case YAEP_ALT:
    printf("YAEP_ALT\n");
    break;
  case _yaep_VISITED:
    printf("_yaep_VISITED\n");
    break;
  case _yaep_MAX:
    printf("_yaep_MAX\n");
    break;
  default:
    printf("NULL\n");
    break;
  }
}

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

static void *parse_alloc_func(int size) {
  void *result;
  result = malloc(size);
  return result;
}

struct yaep_tree_node *parse(const char *description) {
  struct grammar *g;
  struct yaep_tree_node *root;
  int ambiguous_p;

  if ((g = yaep_create_grammar()) == NULL) {
    fprintf(stderr, "yaep_create_grammar: No memory\n");
    exit(1);
  }
  yaep_set_debug_level(g, 0);
  yaep_set_one_parse_flag(g, 0);
  yaep_set_cost_flag(g, 0);
  yaep_set_lookahead_level(g, 2);
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

char *terminal_codes_to_char(int code) {
  switch (code) {
  case 97:
    return "a";
  case 99:
    return "c";
  case 103:
    return "g";
  case 117:
    return "u";
  default:
    return "ERRORRRRRR";
  }
}

/*************************************************************************
 *                YAEP graph traversal functions:                         *
 **************************************************************************/

/**
 * Returns the length of a left_loop size for a particular R rule
 **/
int traverse_parse_tree_for_loop(struct yaep_tree_node *node) {
  switch (node->type) {
  case YAEP_ANODE:
    // FIXME: this requires extra investigation here ... 0 or 1?
    return traverse_parse_tree_for_loop(node->val.anode.children[0]) + 1;
    break;
  case YAEP_TERM:
    return 1;
    break;
  default:
    return -1;
    break;
  }
}

// (akolaitis) debugging
void dump(struct size *s) {
  for (struct size *a = s; a != NULL; a = a->next) {
    printf("%d ->", a->value);
  }
  printf("\n");
}

// return [A x B]
struct size *cartesianProduct(struct size *A, struct size *B) {
  struct size *sizes = NULL;
  if (A == NULL)
    return B;
  if (B == NULL)
    return A;
  // printf("iterA:"); dump(A);
  // printf("iterB:"); dump(B);
  for (struct size *iterA = A; iterA != NULL; iterA = iterA->next) {
    for (struct size *iterB = B; iterB != NULL; iterB = iterB->next) {
      sizes = append_size_if_not_exists(sizes, iterA->value + iterB->value);
    }
  }
  // printf("result:"); dump(sizes);
  return sizes;
}

// [A[0] x A[1] x ... x A[count-1]]
struct size *multiCartesianProduct(struct size **A, int count) {
  struct size *sizes = NULL;
  for (int i = 0; i < count; i++) {
    sizes = cartesianProduct(sizes, A[i]);
  }
  return sizes;
}

/**
 * Returns a list of all the dd size variations of a particular R rule
 **/
struct size *traverse_parse_tree_for_dd(struct yaep_tree_node *node) {
  switch (node->type) {
  case YAEP_ANODE:
    if ((node->val.anode.name)[0] == 'M') {
      struct size **childSizes = malloc(max_dd_size * sizeof(struct size *));
      for (int i = 0; i < max_dd_size; i++) {
        childSizes[i] = traverse_parse_tree_for_dd(node->val.anode.children[i]);
      }
      return multiCartesianProduct(childSizes, max_dd_size);
    } else if ((node->val.anode.name)[0] == 'N') {
      return traverse_parse_tree_for_dd(node->val.anode.children[0]);
    }
    break;
  case YAEP_TERM:
    return create_size(1); // should return a size node of value 1
    break;
  case YAEP_NIL:
    return create_size(0); // should return NULL node
    break;
  case YAEP_ALT:;
    struct size *anode_dds = traverse_parse_tree_for_dd(node->val.alt.node);
    struct size *alt_dds = NULL;
    if (node->val.alt.next != NULL) {
      alt_dds = traverse_parse_tree_for_dd(node->val.alt.next);
    }
    struct size *new_list = concatenate_sunique(anode_dds, alt_dds);
    return new_list;
    break;
  default:
    printf("I think something went really wrong with node type \n");
  }
  return NULL;
}

/**
 * Traverses the high level graph in the means of alternative R rules and only.
 * Returns a list of all identified pesudoknots.
 **/
struct pseudoknot *traverse_parse_tree(struct yaep_tree_node *node) {
  int loop_size;
  struct size *dd_sizes;
  struct pseudoknot *pknots = NULL;
  struct pseudoknot *pseudoknots = NULL;
  if (!node) {
    return pseudoknots;
  }
  switch (node->type) {
  case YAEP_ERROR:
    return NULL;
  case YAEP_NIL:
    return NULL;
  case YAEP_TERM:
    return NULL;
  case YAEP_ANODE:
    if ((node->val.anode.name)[0] != 'R') {
      // This should never resolve ...
      return NULL;
    }
    loop_size = traverse_parse_tree_for_loop(node->val.anode.children[0]);
    dd_sizes = traverse_parse_tree_for_dd(node->val.anode.children[1]);

    while (dd_sizes != NULL) {
      pseudoknots = append_pseudoknot_if_not_exists(
          pseudoknots, create_pseudoknot(loop_size, dd_sizes->value));
      dd_sizes = dd_sizes->next;
    }
    return pseudoknots;
  case YAEP_ALT:
    pknots = traverse_parse_tree(
        node->val.alt.node); // the list comming from the anode

    // a list comming from the alterantive node, in case it exists
    pseudoknots = NULL;
    if (node->val.alt.next != NULL) {
      pseudoknots = traverse_parse_tree(node->val.alt.next);
      // print_node_type(node->val.alt.next);
    }

    pseudoknots = concatenate_punique(pseudoknots, pknots);
    return pseudoknots;
  default:
    printf("I don't give a shit! \n");
    return pseudoknots;
  }
}

char *detect_pseudoknots(char *grammar, char *sequence, int dd_size) {
  ntok = 0;
  input = sequence;
  max_dd_size = dd_size;

  char *buffer;
  size_t size;
  FILE *fp = open_memstream(&buffer, &size);

  struct yaep_tree_node *root = parse(grammar);
  struct pseudoknot *ps = traverse_parse_tree(root);

  for (struct pseudoknot *i = ps; i != NULL; i = i->next) {
    fprintf(fp, "%d, %d", i->left_loop_size, i->dd_size);
    if (i->next != NULL)
      fprintf(fp, "\n");
  }

  fclose(fp);
  return buffer;
}
