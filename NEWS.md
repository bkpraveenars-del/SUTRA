# sutra 0.1.0

First release. Modules 1 and 2 of the SUTRA architecture.

* `sutra_spatial()` fits a two-dimensional tensor-product P-spline spatial model
  and returns genotype BLUEs together with their full error variance-covariance
  matrix.
* `sutra_stagewise()` performs stage-wise integration across environments by
  spectral whitening on the genotype-contrast space, with `fully_efficient`,
  `diagonal` and `unweighted` weighting for direct comparison, and optional
  shrinkage of the variance matrix toward its diagonal.
* `sutra_condition()` reports, per environment, the condition number of the
  matrix that would be inverted, so that an ill-conditioned block is visible
  before the model is fitted.
* `sutra_simulate()` generates multi-environment data with known genotypic
  effects, known variance components and a smooth spatial surface deliberately
  misaligned with block boundaries.
* `sutra_compare_weighting()` runs all three weighting schemes on the same input.

Known limitations, documented and tested:

* The advantage of precision weighting reverses above a crossover
  genotype-by-environment variance of approximately 0.15. Do not weight a
  stage-two model that lacks an interaction term.
* The off-diagonal elements of the variance matrix are second-order and can
  harm when genotype-by-environment variance is large; use `shrink`.
