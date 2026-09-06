#' Simulate a multi-environment trial with known truth
#'
#' Generates plot-level MET data in which the genotypic effects, the spatial
#' surface and every variance component are known. Real data can never validate
#' a pipeline, because the truth is not observed; simulation is the only setting
#' in which "the framework recovers the genetic signal" is a testable claim
#' rather than a hope.
#'
#' The spatial surface is a smooth non-linear trend deliberately misaligned with
#' block boundaries, which is the situation a randomised complete block design
#' cannot represent.
#'
#' @param n_gen number of genotypes.
#' @param n_env number of environments.
#' @param n_rep replicates per genotype per environment.
#' @param n_row,n_col field dimensions per environment; defaults chosen to hold
#'   \code{n_gen * n_rep} plots.
#' @param var_g genotypic variance.
#' @param var_ge genotype-by-environment interaction variance.
#' @param var_spatial variance of the smooth spatial surface.
#' @param var_row,var_col variance of discrete row/column management effects.
#' @param var_e residual (plot) variance. May be a vector of length
#'   \code{n_env} to impose heterogeneous error variance across environments,
#'   which together with imbalance is the condition under which weighting
#'   actually matters.
#' @param dropout proportion of genotype-by-environment cells removed at random
#'   (sparse testing / selection dropout).
#' @param rep_range integer vector of length 2; if given, replication per
#'   genotype is drawn uniformly in this range instead of being fixed at
#'   \code{n_rep}, creating unequal replication.
#' @param crossover if >0, adds a genotype-specific linear response to an
#'   environmental index of this magnitude, generating rank-changing GxE.
#' @param seed RNG seed.
#' @return list with \code{data} (plot records), \code{truth} (the generating
#'   genotype effects, spatial surfaces and variance components) and
#'   \code{settings}.
#' @examples
#' sim <- sutra_simulate(n_gen = 20, n_env = 2, n_rep = 2, seed = 4)
#' str(sim$settings)
#' head(sim$data)
#' # the generating genotypic effects are known, which is the point
#' head(sim$truth$g)
#' @export
sutra_simulate <- function(n_gen = 100, n_env = 6, n_rep = 2,
                           n_row = NULL, n_col = NULL,
                           var_g = 1.0, var_ge = 0.5, var_spatial = 2.0,
                           var_row = 0.10, var_col = 0.08, var_e = 0.6,
                           crossover = 0.6, dropout = 0, rep_range = NULL,
                           seed = 1) {
  set.seed(seed)
  var_e <- if (length(var_e) == 1) rep(var_e, n_env) else var_e
  if (length(var_e) != n_env)
    stop("var_e must have length 1 or n_env")
  # field must hold the largest possible plot count, allowing for unequal
  # replication where some genotypes are sown up to rep_range[2] times
  max_rep <- if (is.null(rep_range)) n_rep else max(rep_range)
  n_plot_max <- n_gen * max_rep
  if (is.null(n_row)) n_row <- ceiling(sqrt(n_plot_max * 1.25))
  if (is.null(n_col)) n_col <- ceiling(n_plot_max / n_row)
  while (n_row * n_col < n_plot_max) n_col <- n_col + 1
  n_plot <- n_gen * n_rep

  gen  <- sprintf("G%03d", seq_len(n_gen))
  g    <- stats::rnorm(n_gen, 0, sqrt(var_g)); names(g) <- gen
  eidx <- seq(-1.5, 1.5, length.out = n_env)          # environmental index
  slope <- stats::rnorm(n_gen, 0, sqrt(crossover)); names(slope) <- gen
  env_main <- stats::rnorm(n_env, 0, 1.2)
  ge <- matrix(stats::rnorm(n_gen * n_env, 0, sqrt(var_ge)), n_gen, n_env,
               dimnames = list(gen, sprintf("E%d", seq_len(n_env))))

  surfaces <- list(); out <- list()
  for (j in seq_len(n_env)) {
    u <- seq_len(n_row); v <- seq_len(n_col)
    U <- outer(u, rep(1, n_col)); V <- outer(rep(1, n_row), v)
    a <- stats::runif(4, -1, 1)
    S <- a[1] * sin(2 * pi * U / n_row * 0.9 + a[2]) +
         a[3] * cos(2 * pi * V / n_col * 0.7) +
         0.8 * exp(-(((U - n_row * 0.3) / (n_row * 0.22))^2 +
                     ((V - n_col * 0.7) / (n_col * 0.24))^2))
    S <- S / stats::sd(as.vector(S)) * sqrt(var_spatial)
    surfaces[[j]] <- S

    reff <- stats::rnorm(n_row, 0, sqrt(var_row))
    ceff <- stats::rnorm(n_col, 0, sqrt(var_col))

    # which genotypes appear in this environment, and how often
    gen_j <- gen
    if (dropout > 0) {
      keep_g <- stats::runif(n_gen) > dropout
      if (sum(keep_g) < 5) keep_g[sample(n_gen, 5)] <- TRUE
      gen_j <- gen[keep_g]
    }
    reps_j <- if (is.null(rep_range)) rep(n_rep, length(gen_j))
              else sample(rep_range[1]:rep_range[2], length(gen_j), replace = TRUE)
    gvec <- rep(gen_j, times = reps_j)
    n_pj <- length(gvec)
    if (n_pj > n_row * n_col)
      stop("environment ", j, ": more plots than field cells; increase n_row/n_col")
    gvec <- gvec[sample(n_pj)]

    cells <- expand.grid(row = u, col = v)[seq_len(n_pj), , drop = FALSE]
    cells <- cells[sample(nrow(cells)), , drop = FALSE]     # randomise layout

    y <- env_main[j] + g[gvec] + slope[gvec] * eidx[j] + ge[gvec, j] +
         S[cbind(cells$row, cells$col)] +
         reff[cells$row] + ceff[cells$col] +
         stats::rnorm(n_pj, 0, sqrt(var_e[j]))

    out[[j]] <- data.frame(
      env = sprintf("E%d", j), genotype = gvec,
      row = cells$row, col = cells$col,
      env_index = eidx[j], y = as.numeric(y),
      stringsAsFactors = FALSE)
  }
  dat <- do.call(rbind, out)
  rownames(dat) <- NULL

  list(data = dat,
       truth = list(g = g, slope = slope, ge = ge, env_main = env_main,
                    env_index = eidx, surfaces = surfaces,
                    varcomp = list(g = var_g, ge = var_ge, spatial = var_spatial,
                                   row = var_row, col = var_col, e = var_e),
                    H2_expected = var_g / (var_g + mean(var_e) / n_rep)),
       settings = list(n_gen = n_gen, n_env = n_env, n_rep = n_rep,
                       n_row = n_row, n_col = n_col, dropout = dropout,
                       rep_range = rep_range, seed = seed))
}
