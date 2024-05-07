# knotify

RNA secondary structure prediction engine.

The monorepo provides a grammar based algorithm to predict the pseudoknots patterns of the secondary structure of any RNA sequence.

Authors: Andrikos Christos, Makris Evaggelos, Pavlatos Christos, Rassias Georgios, Angelos Kolaitis

## Run with Docker

`knotify` is written in a mix of Python and C code. The easiest way to run it is to build the Docker image and then run as follows. We do not currently provide pre-built Docker images, but that might change in the future.

```bash
# build docker image and tag as `knotify:dev`
docker build -t knotify:dev .

# run an example 'rna_analysis'
docker run --rm -it knotify:dev /knotify/bin/rna_analysis --sequence AAAAAACUAAUAGAGGGGGGACUUAGCGCCCCCCAAACCGUAACCCC
```

## Development Environment

`knotify` has binary depedendencies that require running on a Ubuntu 20.04 system. Building the code consists of 2 parts: Setting up a Python 3 virtual environment and building the C parser libraries.

Note that the instructions below will probably not work on newer Ubuntu systems or different Linux distributions (e.g. CentOS). This is because of the [wheel-requirements.txt](./wheel-requirements.txt) pinning the ViennaRNA C library to the ubuntu2004 version, which requires a specific glibc version. In the future, we will build ViennaRNA from source to avoid this restriction.

```bash
$ make deps  # install package dependencies
$ make       # build the parser and setup the virtual environment at ./.venv
```

After installation, make sure to activate the virtual environment before running any of the commands below:

```bash
$ . ./.venv/bin/activate
```

## Execute

### rna_analysis

For a single sequence. See `--help` for a complete list of options:

```bash
$ rna_analysis --sequence AAAAAACUAAUAGAGGGGGGACUUAGCGCCCCCCAAACCGUAACCCC
```

### rna_benchmark

Run a benchmark for a number of cases from a YAML file, saving results in JSON format in `result.json`. See [`cases.yaml`](./cases/cases.yaml) for an example YAML file. See `--help` for a list of available options.

```bash
$ rna_benchmark --cases cases/cases.yaml --max-dd-size 2 --max-stem-allow-smaller 1 --allow-ug --prune-early --parser bruteforce > result.json
$ rna_benchmark --cases cases/cases.yaml --max-dd-size 2 --max-stem-allow-smaller 1 --allow-ug --prune-early --parser yaep > result.json
```

### Calling directly from Python code

See [`scripts/00-example.py`](./scripts/00-example.py) for using `knotify` directly from Python code.

## Unit Tests

We use `pytest` for all unit tests. After enabling the virtual environment, you can run them with:

```bash
$ pytest
```

## Implementation details

The algorighm consits of the two following steps:

    1.  Parse mulitple RNA subsequences to define the potenital pseudoknot structures

    2.  Choose the pseudoknot structure that is considered to be the most stable one. This steps is based on the concept of energy minimization concept.

The core algorithm was initially implemented in python based on the wide-known NLTK package. Due to serious performance issues we moved the parsing into c utilizing the `yaep` parser which is able to parse ambient grammars.

### Scoring

Compare prediction dot bracket with ground truth. Create confusion matrix

| Definition     | Description                                                                                                 |
| -------------- | ----------------------------------------------------------------------------------------------------------- |
| true positive  | Ground truth has a stem here, and I have correctly found that stem (matching the pair, TODO: distance +- 1) |
| true negative  | Ground truth does not have a stem, and I do not have a stem                                                 |
| false positive | I have predicted a stem here, but there is no stem in ground truth                                          |
| false negative | I have predicted no stem, but ground truth has a stem here                                                  |

### Loop sizes and indices

#### Unused right loop:
- inclusive_start_index = left_core_inner + left_loop_stems + 1
- unused_right_loop_size = right_core_outer - inclusive_start_index
- inclusive_end_index = right_core_outer - 1

#### Unused left loop:
- unused_left_loop_size = left_loop_size - right_core_stems
- inclusive_start_index = left_core_outer + 1
- inclusive_end_index = inclusive_start_index + left_loop_size - right_core_stems - 1
