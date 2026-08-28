# Human Skin Permeability QSPR Manuscript Restart

This repository contains the simplified manuscript workflow for an interpretable QSPR analysis of human skin permeability.

The current branch is a clean restart of a previously expanded analysis. The expanded version was archived separately on:

```text
archive/expanded-analysis
```

This branch is intended to rebuild the manuscript around a smaller, clearer analysis pipeline.

## Project aim

The manuscript asks:

> Can an interpretable, compound-aware QSPR model predict human skin permeability better than classical baselines while remaining auditable?

The project is not mainly about producing the most complex possible machine-learning model. The focus is:

- interpretable QSPR modeling;
- compound-aware validation;
- transparent descriptor screening;
- benchmark comparison;
- applicability-domain assessment;
- manuscript-ready tables and figures.

## Current manuscript framing

The manuscript should be written around the following core narrative:

1. Literature-curated human skin permeability datasets are useful but heterogeneous.
2. Repeated compound measurements can make ordinary row-wise validation optimistic.
3. Leave-one-compound-out cross-validation is used as the primary model-selection design.
4. Descriptor redundancy is screened before model selection to avoid unstable combinations in interpretable parametric models.
5. The selected interpretable QSPR model is compared with classical and flexible benchmark models.
6. Applicability-domain and validation-sensitivity analyses define where the model is most defensible.

## Current simplified workflow

The intended active script structure is:

```text
R/
├── 00_config.R
├── 01_prepare_data.R
├── 02_screen_descriptors.R
├── 03_select_qspr_model.R
├── 04_benchmark_models.R
├── 05_define_applicability_domain.R
├── 06_validation_sensitivity.R
├── 07_model_interpretation.R
├── 08_supplementary_diagnostics.R
├── 20_make_manuscript_tables.R
├── 30_make_manuscript_figures.R
└── 99_run_all.R
```

Scripts should map directly to manuscript Methods sections.

## Script roles

### `00_config.R`

Defines common paths, candidate descriptors, endpoint names, model-selection settings, output paths, and shared helper settings.

### `01_prepare_data.R`

Imports and cleans the human skin permeability dataset.

Expected outputs include:

- cleaned modeling dataset;
- cleaning-flow summary;
- removed-row logs;
- data-cleaning table source.

### `02_screen_descriptors.R`

Screens candidate descriptors before model selection.

This script should generate:

- candidate descriptor correlation heatmap;
- full candidate descriptor VIF table;
- red-flag predictor combinations to exclude from candidate model search;
- soft-warning descriptor combinations for interpretation only.

Important principle:

```text
Do not remove descriptors globally.
Only skip candidate formulas containing highly redundant descriptor combinations.
```

Expected main figure:

```text
Figure 2. Candidate descriptor correlation structure.
```

Expected output files:

```text
results/descriptor_screening/
figures/figure2_descriptor_correlation.png
tables/tableS1_descriptor_redundancy_summary.csv
```

### `03_select_qspr_model.R`

Performs compound-aware interpretable model selection.

Primary validation design:

```text
leave-one-compound-out cross-validation
```

This script should read the red-flag predictor combinations from descriptor screening and skip candidate formulas containing those combinations.

Expected outputs include:

- selected model formula;
- LOCO-CV predictions;
- model-selection summary;
- final model coefficient table source;
- observed-versus-predicted figure source.

### `04_benchmark_models.R`

Compares the selected interpretable QSPR model with benchmark models.

Benchmarks may include:

- null mean model;
- Potts-Guy-style model;
- linear selected-predictor model;
- selected-predictor random forest;
- RDKit descriptor benchmark, if retained.

The benchmark section should not overclaim. The interpretable model only needs to improve over classical baselines and approach flexible models.

### `05_define_applicability_domain.R`

Defines and evaluates the descriptor-space applicability domain.

Expected outputs include:

- domain membership table;
- coverage summary;
- error summary inside and outside the domain;
- applicability-domain figure source.

### `06_validation_sensitivity.R`

Evaluates sensitivity to validation design.

Core planned comparison:

```text
LOCO-CV vs repeated row-wise CV
```

Leave-one-reference-out validation should only be included if it is fully scripted and exported.

### `07_model_interpretation.R`

Contains model-interpretation outputs that are useful for manuscript figures or supplement.

Possible outputs:

- partial-effect plots;
- interaction plots;
- selected ablation results.

This script should stay focused. Do not let interpretation diagnostics overwhelm the main manuscript.

### `08_supplementary_diagnostics.R`

Contains optional analyses that support, but do not drive, the main manuscript.

Examples:

- coefficient stability;
- residual diagnostics;
- high-error compounds;
- heteroscedastic uncertainty analysis.

These should be supplementary unless directly needed in the main Results.

### `20_make_manuscript_tables.R`

Collects analysis outputs and writes final manuscript-ready tables.

Main tables should include:

1. Data cleaning and attrition.
2. Candidate model selection summary.
3. Final selected model coefficients.
4. Benchmark model performance.
5. Applicability-domain summary.
6. Validation-design sensitivity summary.

Supplementary tables should include:

- descriptor redundancy summary;
- red-flag predictor combinations;
- additional diagnostics.

### `30_make_manuscript_figures.R`

Creates or collects manuscript-ready figures.

Main figures should include:

1. Workflow figure.
2. Descriptor correlation structure.
3. Observed versus LOCO-CV predicted values.
4. Benchmark comparison.
5. Applicability-domain/error figure.
6. Model interpretation figure, if needed.

Supplementary figures should be kept separate.

### `99_run_all.R`

Runs the simplified workflow in order.

Example intended order:

```r
source("R/00_config.R")
source("R/01_prepare_data.R")
source("R/02_screen_descriptors.R")
source("R/03_select_qspr_model.R")
source("R/04_benchmark_models.R")
source("R/05_define_applicability_domain.R")
source("R/06_validation_sensitivity.R")
source("R/07_model_interpretation.R")
source("R/08_supplementary_diagnostics.R")
source("R/20_make_manuscript_tables.R")
source("R/30_make_manuscript_figures.R")
```

Unfinished optional scripts can be commented out.

## Current simplification rule

The previous analysis became too broad. The restart should prioritize manuscript completion.

Keep in the main manuscript:

- data cleaning;
- descriptor screening;
- LOCO-CV model selection;
- benchmark comparison;
- applicability-domain analysis;
- validation-design sensitivity.

Move to supplement or archive:

- extensive coefficient stability analysis;
- extensive residual diagnostics;
- heteroscedastic uncertainty modeling;
- excessive ablation variants;
- excessive partial-effect figures;
- leave-one-reference-out validation unless fully implemented;
- redundant RDKit benchmarking details.

## Descriptor screening plan

Descriptor screening should be simple and auditable.

Use pairwise Pearson correlations to identify direct redundancy.

Use VIF to identify multivariable redundancy.

Planned logic:

```text
Hard red flags:
- pairwise |r| >= 0.85
- extreme VIF group, e.g. VIF >= 100

Soft warnings:
- 0.70 <= pairwise |r| < 0.85
- 5 <= VIF < 100
```

Hard red-flag combinations are excluded from candidate model search.

Soft-warning combinations are retained but interpreted cautiously.

Current example from earlier screening:

```text
MWa + MVh
```

was strongly correlated and is expected to be a hard red-flag pair.

The solubility/partitioning descriptors:

```text
logKowb
LogSaqd
LogSoce
```

showed extreme VIF behavior in the full descriptor pool and should be handled as a group-level redundancy issue.

## Manuscript structure

The manuscript should use this simplified structure:

```text
Introduction
Methods
Results
Discussion
Conclusion
```

Suggested Results subsections:

```text
1. Cleaned modeling dataset
2. Descriptor redundancy and candidate predictor structure
3. Compound-level model selection
4. Benchmark comparison
5. Applicability-domain analysis
6. Validation-design sensitivity
```

Avoid adding many extra Results subsections unless they directly support the central manuscript claim.

## Provisional manuscript claim

A cautious version of the central claim:

> This study develops an interpretable, compound-aware QSPR workflow for human skin permeability and evaluates it against classical and flexible benchmark models. The workflow emphasizes leakage-aware validation, descriptor redundancy control, and applicability-domain assessment rather than maximum black-box predictive complexity.

## Data availability note

The raw dataset and some derived/private files may not be included in the public repository.

The README should not imply that the full analysis is reproducible from public files unless all required data inputs are included.

Use language such as:

```text
Analysis scripts and generated non-sensitive outputs are included. Raw data files are omitted where redistribution is not appropriate.
```

## Current priority

Immediate work should proceed in this order:

1. Finalize simplified script names and workflow.
2. Rebuild `00_config.R`.
3. Rebuild `01_prepare_data.R`.
4. Rebuild `02_screen_descriptors.R`.
5. Rebuild `03_select_qspr_model.R`.
6. Re-run the core workflow.
7. Check whether the selected model changes.
8. Rebuild manuscript tables and figures.
9. Rewrite Results using the final outputs.
10. Finalize Discussion, Abstract, and data/code availability statement.

## Notes for continuation

This repository is currently in a cleanup/restart phase. Do not infer that all old scripts are still part of the active manuscript workflow. The current goal is not to expand the analysis. The goal is to produce a clean, defensible manuscript with a small number of well-aligned scripts, tables, and figures.

When working with this repository, prioritize:

- simplifying;
- aligning scripts with Methods sections;
- avoiding overclaiming;
- keeping diagnostics supplementary;
- producing manuscript-ready outputs.

Do not re-expand the workflow unless explicitly required.
