test_that("Gate 1: Module 1 reproduces the SpATS reference heritability", {
  skip_if_not_installed("SpATS")
  data(wheatdata, package = "SpATS", envir = environment())
  d <- wheatdata[!is.na(wheatdata$yield), ]
  d$R <- as.factor(d$row); d$C <- as.factor(d$col)
  ref <- SpATS::SpATS(response = "yield",
                      spatial = ~ SpATS::PSANOVA(col, row, nseg = c(10, 20)),
                      genotype = "geno", genotype.as.random = TRUE,
                      fixed = ~ colcode + rowcode, random = ~ R + C, data = d,
                      control = list(tolerance = 1e-3, monitoring = 0))
  s <- sutra_spatial(d, response = "yield", genotype = "geno",
                     col = "col", row = "row", fixed = ~ colcode + rowcode,
                     nseg = c(10, 20), genotype.as.random = TRUE)
  expect_equal(s$H2, SpATS::getHeritability(ref), tolerance = 1e-6)
})

test_that("Sigma is symmetric, positive definite, and named", {
  sim <- sutra_simulate(n_gen = 25, n_env = 1, n_rep = 2, seed = 7)
  s <- sutra_spatial(sim$data, response = "y", genotype = "genotype",
                     col = "col", row = "row")
  expect_true(isSymmetric(unname(s$Sigma), tol = 1e-8))
  expect_true(s$diagnostics$Sigma_pd)
  expect_equal(rownames(s$Sigma), s$blues$genotype)
})

test_that("the contrast projection removes the intercept direction", {
  sim <- sutra_simulate(n_gen = 30, n_env = 1, n_rep = 2, seed = 11)
  s <- sutra_spatial(sim$data, response = "y", genotype = "genotype",
                     col = "col", row = "row")
  V <- s$Sigma; n <- nrow(V)
  e <- eigen(V, symmetric = TRUE)
  # the dominant eigenvector should BE the constant vector: this is the
  # intercept variance that must not be used as a weight
  expect_gt(abs(sum(e$vectors[, 1] * rep(1 / sqrt(n), n))), 0.9)
  P <- diag(n) - matrix(1 / n, n, n)
  ec <- eigen(P %*% V %*% P, symmetric = TRUE)$values
  # after projection the last eigenvalue is the removed null direction
  expect_lt(abs(ec[n]) / max(ec), 1e-8)
})

test_that("unweighted uses one global scalar, not per-environment variance", {
  sim <- sutra_simulate(n_gen = 25, n_env = 3, n_rep = 2,
                        var_e = c(0.2, 1.0, 4.0), seed = 13)
  s1 <- lapply(unique(sim$data$env), function(e)
    sutra_spatial(sim$data[sim$data$env == e, ], response = "y",
                  genotype = "genotype", col = "col", row = "row", env = e))
  u <- sutra_stagewise(s1, weighting = "unweighted")
  f <- sutra_stagewise(s1, weighting = "fully_efficient")
  # with strongly heterogeneous variance the two must NOT coincide;
  # if they do, the unweighted baseline is secretly weighted
  expect_lt(cor(u$effects$value, f$effects$value), 0.999)
})

test_that("weighting recovers truth better under heterogeneity than under balance", {
  gain <- function(var_e, seed) {
    # crossover MUST be 0 here: with crossover GxE present, precision weighting
    # narrows the effective environmental sample and the expected sign reverses
    # (see inst/scripts/gate4_crossover.R).  Leaving the simulator's default
    # crossover in place makes this test measure GxE, not weighting.
    sim <- sutra_simulate(n_gen = 40, n_env = 6, n_rep = 2, var_ge = 0,
                          crossover = 0, var_e = var_e, rep_range = c(1, 3),
                          seed = seed)
    s1 <- lapply(unique(sim$data$env), function(e)
      sutra_spatial(sim$data[sim$data$env == e, ], response = "y",
                    genotype = "genotype", col = "col", row = "row", env = e))
    cmp <- sutra_compare_weighting(s1)
    tr <- sim$truth$g[cmp$effects$genotype]
    cor(cmp$effects$fully_efficient, tr) - cor(cmp$effects$unweighted, tr)
  }
  het <- mean(vapply(1:4, function(i) gain(c(.2,.4,.8,1.5,2.5,4), 900 + i), numeric(1)))
  hom <- mean(vapply(1:4, function(i) gain(1.0, 900 + i), numeric(1)))
  expect_gt(het, hom)
})

test_that("Gate 4: precision weighting reverses under crossover GxE", {
  gain <- function(crossover, seed) {
    sim <- sutra_simulate(n_gen = 40, n_env = 6, n_rep = 2, var_ge = 0,
                          crossover = crossover, var_e = c(.2,.4,.8,1.5,2.5,4),
                          rep_range = c(1, 3), seed = seed)
    s1 <- lapply(unique(sim$data$env), function(e)
      sutra_spatial(sim$data[sim$data$env == e, ], response = "y",
                    genotype = "genotype", col = "col", row = "row", env = e))
    cmp <- sutra_compare_weighting(s1)
    tr <- sim$truth$g[cmp$effects$genotype]
    cor(cmp$effects$fully_efficient, tr) - cor(cmp$effects$unweighted, tr)
  }
  none <- mean(vapply(1:4, function(i) gain(0.0, 900 + i), numeric(1)))
  high <- mean(vapply(1:4, function(i) gain(1.2, 900 + i), numeric(1)))
  expect_gt(none, 0)     # weighting helps when environments agree
  expect_lt(high, 0)     # and hurts when they do not
  expect_gt(none - high, 0.1)
})

test_that("sutra_condition flags ill-conditioned blocks", {
  sim <- sutra_simulate(n_gen = 30, n_env = 2, n_rep = 2, seed = 21)
  s1 <- lapply(unique(sim$data$env), function(e)
    sutra_spatial(sim$data[sim$data$env == e, ], response = "y",
                  genotype = "genotype", col = "col", row = "row", env = e))
  cd <- sutra_condition(s1)
  expect_true(all(c("env","condition","risk") %in% names(cd)))
  expect_true(all(cd$condition > 1))
})
