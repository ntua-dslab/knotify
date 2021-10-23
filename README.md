RNA secondary structure prediction engine
-----------------------------------------

The monorepo provides a grammar based algorithm to  predict the   pseudoknots patters of the secondary structure of any RNA sequence.

Authors: Andrikos Christos, Makris Evaggelos, Pavlatos Christos, Rassias Georgios, Aggelos Kolaitis

Implementation details
----------------------

The algorighm consits of the two following steps:

    1.  Parse mulitple RNA subsequences to define the potenital pseudoknot structures

    2.  Choose the pseudoknot structure that is considered to be the most stable one. This steps is based on the concept of energy minimization concept.


The core algorithm was initially implemented in python based on the wide-known NLTK package. Due to serious performance issues we moved the parsing into c utilizing the `yaep` parser which is able to parse ambient grammars. The operation goes as it follows:

    1.
    2.
    2.
    2.
    2.
    2.


## Scoring

Compare prediction dot bracket with ground truth. Create confusion matrix

| Definition     | Description                                                                                                 |
| -------------- | ----------------------------------------------------------------------------------------------------------- |
| true positive  | Ground truth has a stem here, and I have correctly found that stem (matching the pair, TODO: distance +- 1) |
| true negative  | Ground truth does not have a stem, and I do not have a stem                                                 |
| false positive | I have predicted a stem here, but there is no stem in ground truth                                          |
| false negative | I have predicted no stem, but ground truth has a stem here                                                  |

## Development Environment

Building the code consists of 2 parts: Setting up a Python 3 virtual environment and building the C parser library.

```bash
$ make deps  # install package dependencies
$ make       # build the parser and setup the virtual environment at ./.venv
```

Run for all cases and save result to `result.yaml`:

```bash
$ ./.venv/bin/rna_benchmark --cases cases.yaml --parser bruteforce --max-dd-size 2 --max-stem-allow-smaller 1 --allow-ug --prune-early > result.json
$ ./.venv/bin/rna_benchmark --cases cases.yaml --parser yaep --max-dd-size 2 --max-stem-allow-smaller 1 --allow-ug --prune-early > result.json
```

## Execute

For a single sequence. See `--help` for a complete list of options:

```bash
$ ./.venv/bin/rna_analysis --parser yaep AAAAAACUAAUAGAGGGGGGACUUAGCGCCCCCCAAACCGUAACCCC
```

Run benchmark for a number of cases, print output in YAML file. See [`cases.yaml`](./cases.yaml) for an example YAML file:

```bash
$ ./.venv/bin/rna_benchmark --parser yaep [OPTIONS] --cases cases.yaml > results.yaml
```

## Unit Tests

Activate virtual environment and run with Pytest:

```bash
$ ./.venv/bin/pytest
```

## Loop sizes and indices

#### Unused right loop:
- inclusive_start_index = left_core_inner + left_loop_stems + 1
- unused_right_loop_size = right_core_outer - inclusive_start_index
- inclusive_end_index = right_core_outer - 1

### Unused left loop:
- unused_left_loop_size = left_loop_size - right_core_stems
- inclusive_start_index = left_core_outer + 1
- inclusive_end_index = inclusive_start_index + left_loop_size - right_core_stems - 1
