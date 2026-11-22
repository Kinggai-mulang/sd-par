# Sequential decorrelation for GNSS partial ambiguity resolution

MATLAB code for transformation inheritance between nested ambiguity subsets.
When one ambiguity is excluded, SD constructs a unimodular transformation for
the retained subset from the preceding reduction. LAMBDA reduction then starts
from the inherited covariance. The transformed problem is equivalent to the
integer least-squares (ILS) problem for that specified subset.

The repository accompanies the manuscript *Sequential Decorrelation Based on
Transformation Inheritance: Exploiting Covariance Relationships in GNSS Partial
Ambiguity Resolution*, by Kinggai Mula, Xiaohong Zhang, Wanke Liu, Ying Liu,
and Liguo Lu.

## Requirements

- MATLAB R2024a, the version used for the experiments and release checks.
- LAMBDA 4.0 for MATLAB, downloaded separately from its authors.
- Write access to the checkout for generated inputs and results.

[Download LAMBDA 4.0](https://www.tudelft.nl/en/ceg/about-faculty/departments/geoscience-remote-sensing/research/lambda/lambda/)

Set `config.lambdaDirectory` near the top of `matlab/sd_experiment_config.m`
to the toolbox root containing `LAMBDA_toolbox`. Alternatively, set the path
for the current MATLAB session before running an experiment.

```matlab
setenv('SD_LAMBDA_DIRECTORY', 'D:/software/LAMBDA4');
```

An explicitly configured `config.lambdaDirectory` takes precedence over the
environment variable. Keep the toolbox's copyright notices and follow its terms.
Attribution and the counted reference implementation are described in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Quick start

Open this repository as the MATLAB current folder, configure the toolbox path,
and run the default small-sample dimension experiment.

```matlab
run_experiments
```

The default configuration selects experiment 2 in `pilot` mode with three
independent covariance matrices per setting. The script generates inputs,
measures CD and SD, writes tables, and exports figures. Results are placed under
`output/numerical/pilot/random_order` and `output/figures/pilot/random_order`.

For a quick check of every primary experiment, edit these configuration entries
before running `run_experiments`.

```matlab
config.experimentSwitches = logical([1 1 1 1 1 1 1]);
config.experimentMode = "pilot";
config.stageDimensionValues = [10 12];
config.dimensionValues = [10 12];
config.conditionDimension = 10;
config.conditionExponents = [2 4];
config.gb.numSatellites = [5 6];
```

## Experiment configuration

Edit `matlab/sd_experiment_config.m`. The seven entries of `experimentSwitches`
select the following comparisons.

| Entry | Comparison | Main output table |
| --- | --- | --- |
| 1 | LAMBDA reduction and search stages | `exp01_lambda4_stage_dimension.tsv` |
| 2 | Single-step reduction over random-LDLT dimensions | `exp02_cd_sd_dimension.tsv` |
| 3 | Single-step reduction over prescribed condition numbers | `exp03_cd_sd_condition.tsv` |
| 4 | Single-step reduction and search with GPS geometry covariances | `exp04_gb.tsv` |
| 5 | IGT and swap counts for random-LDLT inputs | `exp05_reduction_operations_dimension.tsv` |
| 6 | Repeated ILS from the starting dimension down to one | `exp06_random_dimension_full_par.tsv` |
| 7 | Random, maximum-variance, and minimum-variance exclusions | `removal_comparison/removal_policy_dimension.tsv` |

`experimentMode = "full"` uses 1,000 matrices per setting. Restore the dimension
grids shown in the configuration file for the manuscript's parameter ranges.
High-dimensional random-LDLT inputs can require long integer searches.
`useCachedData = false` regenerates inputs for the selected families. Use this
setting after changing covariance, seed, sample-count, or ordering parameters.
The script validates cached inputs before reuse. Experiment 7 has a separate
`useCachedResults` switch for its timing checkpoints.

`defaultOrdering = "random-order"` randomly permutes each covariance and
excludes its final ambiguity. `"variance-ordered"` sorts original marginal
variances in ascending order before excluding the last ambiguity. Experiment 7
compares its three choices on the same random-order inputs.

### Additional comparisons

The following functions generate their required covariance inputs and use the
same configuration grids and ordering as the primary experiments.

```matlab
run_gb_chain_comparison("pilot")
run_single_step_counts("geometry-based", "pilot")
run_single_step_counts("random-condition", "pilot")
run_single_step_counts("random-dimension", "pilot")
```

GB sequence timings are written beside the primary tables as
`exp06_gb_full_par.tsv`. Count summaries and per-matrix counts are written under
`output/numerical/<mode>/<ordering>/operation_counts`.
After generating experiment 6 and the GB sequences in the same mode, run
`plot_chain_comparison("pilot")` to plot total-runtime ratios and search shares.

`plot_experiments` redraws the primary figures from saved tables.
`export_paper_results` copies available tables and figures to `paper_results`.
The exported collection reflects the experiments that have actually been run.

## Using the SD algorithm

The original-coordinate index of the ambiguity to exclude is passed explicitly.

```matlab
addpath('matlab');
setup_sd_experiment_paths();
Q = [4 1 0.5; 1 3 0.2; 0.5 0.2 2];
parent = lambda4_reduce_covariance(Q);
child = sd_lambda4_reduce_child(Q, parent, 2);
% child.Z * child.Qsubset * child.Z' equals child.Qz up to roundoff.
```

`sd_inherit_child_transform` performs inheritance. `sd_lambda4_reduce_child`
adds subsequent LAMBDA reduction. `compare_lambda4_child_reduction` provides
a paired timing comparison for one covariance and exclusion index.
The transformed float vector is `child.Z * aFloatSubset`; map an integer
candidate back using `child.inverseZ`.

## Measurement scope

Each single-step timing starts after the parent covariance has been reduced.
CD reduces the retained covariance directly. SD performs inheritance and then
LAMBDA reduction. Separate counted implementations record FLOPs, nonzero IGTs,
and swaps, including the inheritance operations. FLOPs count additions,
subtractions, multiplications, and divisions. The explicit inverse is assigned
an LU-plus-solves cost of `8*n^3/3`; this is an arithmetic cost model.

Complete sequences start with full-set reduction and solve every successive ILS
problem down to dimension one. Total timing covers reduction, search, ambiguity
transformations, and subset updates. The best back-transformed candidate and
its objective value are compared between CD and SD. Search requests two
candidates, while the sequence consistency summary records the best candidate.
These benchmarks use a given exclusion sequence without an acceptance test.

The covariance generators cover random-LDLT, condition-controlled SPD, and an
idealized triple-frequency GPS short-baseline model. GPS skies are sampled
randomly with two receivers and single-epoch observations. The largest satellite
counts extend the dimension range as an idealized setting.

The public package covers the MATLAB algorithms and simulation comparisons.
The manuscript's field benchmark uses float solutions from KinPOS and a C
implementation maintained separately. Field input data are available through
the corresponding author under the manuscript's data availability statement.
Runtimes vary with MATLAB version, machine, and inputs.

## Repository structure

```text
matlab/algorithms/   Transformation inheritance and reduction
matlab/generators/   Covariance models and input caching
matlab/metrics/      Timings and counted implementations
matlab/runners/      Experiment workflows and table output
matlab/plotting/     Shared plotting functions
data/               Generated input caches
output/             Generated measurements and figures
```

## Validation

The release was checked in MATLAB R2024a using a clean Git clone and a reset
MATLAB path. Small-sample checks exercised all seven primary experiments,
GB sequences, all three operation-count families, plotting, exports, cache
reuse, and variance-ordered sequences. The checks used three matrices per
setting with random dimensions 10 and 12, condition numbers 4 and 16 at
dimension 10, and GPS satellite counts 5 and 6. They verify the execution
workflow; the manuscript's 1,000-matrix performance runs are separate.

## License and citation

Original contributions are distributed under the [MIT License](LICENSE),
copyright 2026 Kinggai Mula and contributors. Cite the SD manuscript when using
this implementation in research. Cite the LAMBDA 4.0 publication when using its
toolbox.

Massarweh, L., Verhagen, S., and Teunissen, P. J. G. (2025).
New LAMBDA toolbox for mixed-integer models: estimation and evaluation.
GPS Solutions, 29, 14. https://doi.org/10.1007/s10291-024-01738-z
