# sutra — Stage-wise Unified Trait–Reaction-norm Architecture

A variance-preserving analytical continuum for multi-environment trial (MET)
analysis. Spatially adjusted genotype means are passed to a stage-wise linear
mixed model **together with their error variance–covariance structure**, rather
than as bare point estimates.

Status: **Modules 1 and 2 implemented and validated. Modules 3 and 4 not yet
written.** Every claim below is reproducible from `inst/scripts`.

## Install

```r
# dependencies: SpATS, lme4, Matrix
R CMD INSTALL sutra
```

## Use

```r
library(sutra)

# Module 1 -- one environment at a time
s1 <- lapply(unique(dat$env), function(e)
  sutra_spatial(dat[dat$env == e, ], response = "yield", genotype = "geno",
                col = "col", row = "row", env = e))

# Look before you leap: is any block ill-conditioned?
sutra_condition(s1)

# Module 2 -- stage-wise integration
fit <- sutra_stagewise(s1, weighting = "fully_efficient")

# What is the weighting actually worth on YOUR data?
sutra_compare_weighting(s1)$summary
```

## What the validation shows — including where the method fails

| Gate | Question | Result |
|---|---|---|
| 1 | Does Module 1 reproduce SpATS? | H² = 0.77000 vs 0.77000 on `wheatdata`; variance components identical to machine precision |
| 2 | Does the rotation reproduce a single-stage fit? | r = 0.99 against `lme4` on raw plot data |
| 3 | Is variance preservation worth anything? | +0.004 under balance (nothing); **+0.062** under heterogeneous variance, unequal replication and 30% dropout |
| 4 | When does it stop working? | **The sign reverses above a crossover-GxE variance of ≈0.15**, reaching −0.228 at 1.2 |

Three findings that contradict the usual framing, and that you should know before
using this package:

1. **The raw variance–covariance matrix of genotype means must not be used as a
   weight.** Its dominant eigenvector is the constant vector — intercept
   uncertainty, 42% of total variance on `wheatdata`, correlation 0.99968 with
   the constant vector. Whitening by it makes the fully efficient analysis
   *worse* than an unweighted one. `sutra_stagewise()` projects onto the contrast
   space by default (`contrast_space = TRUE`). Once projected, the genuine
   off-diagonal correlation is 2.7%.

2. **The off-diagonals are second-order.** They help by 0.002–0.006 without GxE
   and *hurt* in four of five scenarios once GxE is present, because they must be
   estimated. Carry the standard errors always; carry the covariances only after
   checking. `shrink` interpolates between the two.

3. **Do not weight a stage-two model that has no GxE term.** Precision weighting
   narrows the effective environmental sample. With crossover GxE this makes the
   estimate of the genotype main effect worse, badly. This is a limitation of the
   current two-module implementation, not of the architecture — but it is real
   today.

## Reproduce

```sh
Rscript inst/scripts/gate2_stagewise_vs_singlestage.R
Rscript inst/scripts/gate3_definitive.R
Rscript inst/scripts/gate4_crossover.R
```

## Licence

GPL-3.
