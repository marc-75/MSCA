## Version: 1.0.0
### 2025-02-19

Contained very basic information and only the make_state_matrix that can be useful

## Version: 1.1.0
### 2025-06-02

### Changes

- Added:

    - `fast_clara_jaccard`	Fast CLARA-like clustering using Jaccard dissimilarity
    - `get_cluster_sequences`	Extract sequences of length k within clusters
    - `sequence_stats()` to compute frequency, conditional probability, and relative risk of sequences by cluster

- A vignette displaying a basic workflow is available 

### Bug Fixes

- `fast_jaccard_dist`	were corrected 

### Comments on current version 

This version will hopeful permit to run basic analyses of electronic health records. Further examples and functions are expected soon.  


## Version: 1.2.0
### 2025-06-15

### Changes

Corrected, modified and integrated `sequence_stats()` and `get_cluster_sequences`. The sequence frenquencies are computed by patient not on the total number of sequences

Corrected some spelling in the vignette.

### To do

Plot methods for sequences


## Version: 1.2.0
### 2025-06-15

### Changes

Corrected, modified and integrated `sequence_stats()` and `get_cluster_sequences`. The sequence frenquencies are computed by patient not on the total number of sequences

Corrected some spelling in the vignette.

### To do

Plot methods for sequences

## Version: 1.3.0
### Started 2025-06-15

### Changes

Implement CLARANS as an option for the `fast_clara_jaccard` function

## Version: 1.4.0
### 2026-09-27

### Changes

- `fast_clara_jaccard()` now returns the function call and all argument values
  (including defaults) in `$call`, to make analyses easier to reproduce.
- `get_cluster_sequences()`: clusters are now returned in sorted order.
- Vignette: typos corrected.

### Bug Fixes

- `fast_clara_jaccard()`: fixed medoid indexing (fastkmedoids returns 0-based
  indices), which could drop one medoid or select the neighbouring patient.

### Comments on current version

MSCA was archived from CRAN on 2026-03-24 because its dependency
`fastkmedoids` was temporarily failing CRAN checks. This version restores
MSCA on CRAN and makes it independent of `fastkmedoids` for its default use.
