#' Module 1: micro-environmental spatial deconvolution
#'
#' Fits a two-dimensional tensor-product P-spline (PS-ANOVA) spatial model to a
#' single field trial and returns spatially adjusted genotype BLUEs together
#' with their FULL, non-diagonal error variance-covariance matrix.
#'
#' The retention of the off-diagonal elements of \code{Sigma} is the whole point
#' of the module: they are what a conventional export of "adjusted means" throws
#' away, and what \code{\link{sutra_stagewise}} needs in order to reproduce a
#' single-stage analysis.
#'
#' Genotype is fitted as FIXED by default. This is deliberate: the module exists
#' to emit BLUEs. Fitting genotype as random here and again in stage two shrinks
#' the genotypic estimates twice and biases genetic variance downward
#' (Holland & Piepho 2024, "Don't BLUP twice", G3 14:jkae250).
#'
#' @param data data.frame of plot records for ONE environment.
#' @param response character; name of the response column.
#' @param genotype character; name of the genotype factor column.
#' @param col,row character; names of the numeric field column/row coordinates.
#' @param fixed,random one-sided formulae for additional fixed/random terms
#'   (e.g. \code{~ replicate}); \code{NULL} for none.
#' @param nseg integer vector of length 2, number of spline segments for
#'   (col, row). Defaults to roughly half the number of distinct positions,
#'   which is the SpATS default recommendation.
#' @param env character; a label identifying this environment.
#' @param genotype.as.random logical; keep FALSE unless you know why you want
#'   BLUPs out of stage one.
#' @param tolerance REML convergence tolerance.
#' @return An object of class \code{sutra_stage1}: a list with \code{blues},
#'   \code{Sigma}, \code{H2}, \code{varcomp}, \code{diagnostics} and \code{fit}.
#' @examples
#' # small simulated trial, one environment
#' sim <- sutra_simulate(n_gen = 20, n_env = 1, n_rep = 2, seed = 1)
#' s <- sutra_spatial(sim$data, response = "y", genotype = "genotype",
#'                    col = "col", row = "row")
#' s
#' # the FULL error variance-covariance matrix, not just its diagonal
#' dim(s$Sigma)
#' @export
sutra_spatial <- function(data, response, genotype, col, row,
                          fixed = NULL, random = NULL, nseg = NULL,
                          env = "E1", genotype.as.random = FALSE,
                          tolerance = 1e-3) {
  stopifnot(is.data.frame(data))
  for (v in c(response, genotype, col, row)) {
    if (!v %in% names(data)) stop("column not found in data: ", v)
  }
  d <- data
  d[[genotype]] <- as.factor(d[[genotype]])
  d[[col]] <- as.numeric(as.character(d[[col]]))
  d[[row]] <- as.numeric(as.character(d[[row]]))
  d <- d[!is.na(d[[response]]), , drop = FALSE]

  ncol_u <- length(unique(d[[col]]))
  nrow_u <- length(unique(d[[row]]))
  if (is.null(nseg)) nseg <- c(max(3, floor(ncol_u / 2)), max(3, floor(nrow_u / 2)))

  # discrete row/column random effects (management artefacts)
  d$.Rf <- as.factor(d[[row]])
  d$.Cf <- as.factor(d[[col]])
  rand_terms <- c(".Rf", ".Cf")
  if (!is.null(random)) rand_terms <- c(rand_terms, all.vars(random))
  rand_fml <- stats::as.formula(paste("~", paste(rand_terms, collapse = " + ")))

  spat <- stats::as.formula(
    sprintf("~ PSANOVA(%s, %s, nseg = c(%d, %d))", col, row, nseg[1], nseg[2]))
  # PSANOVA() must be resolvable when SpATS evaluates the spatial formula,
  # whether or not the user has attached the SpATS package
  environment(spat) <- environment()

  fit <- SpATS::SpATS(
    response = response, spatial = spat, genotype = genotype,
    genotype.as.random = genotype.as.random,
    fixed = fixed, random = rand_fml, data = d,
    control = list(tolerance = tolerance, monitoring = 0))

  pred <- stats::predict(fit, which = genotype, return.vcov.matrix = TRUE)
  V <- attr(pred, "vcov")
  # SpATS may return a 'spam' sparse object; coerce to a dense base matrix
  V <- matrix(as.numeric(as.matrix(V)), nrow = nrow(pred), ncol = nrow(pred))
  gl <- as.character(pred[[genotype]])
  dimnames(V) <- list(gl, gl)
  V <- (V + t(V)) / 2                      # enforce exact symmetry

  blues <- data.frame(
    env      = env,
    genotype = gl,
    value    = as.numeric(pred$predicted.values),
    se       = sqrt(diag(V)),
    stringsAsFactors = FALSE)

  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  off <- V[upper.tri(V)]
  diagnostics <- list(
    eff_dim          = fit$eff.dim,
    n_obs            = nrow(d),
    n_genotypes      = length(gl),
    psanova_nseg     = nseg,
    sigma2_residual  = fit$psi[1],
    Sigma_offdiag_ratio = mean(abs(off)) / mean(diag(V)),
    Sigma_min_eigen  = min(ev),
    Sigma_condition  = max(ev) / min(ev),
    Sigma_pd         = all(ev > 0),
    semivariogram    = tryCatch(
      sutra_semivariogram(d, response, col, row, fit), error = function(e) NULL))

  out <- list(
    env = env, blues = blues, Sigma = V,
    H2 = tryCatch(SpATS::getHeritability(fit), error = function(e) NA_real_),
    varcomp = fit$var.comp, diagnostics = diagnostics, fit = fit,
    genotype.as.random = genotype.as.random)
  class(out) <- "sutra_stage1"
  out
}

#' @export
print.sutra_stage1 <- function(x, ...) {
  cat("<sutra Module 1: spatial deconvolution>\n")
  cat("  environment      :", x$env, "\n")
  cat("  plots / genotypes:", x$diagnostics$n_obs, "/", x$diagnostics$n_genotypes, "\n")
  cat("  heritability     :", if (is.na(x$H2)) "NA" else round(x$H2, 4), "\n")
  cat("  Sigma            :", nrow(x$Sigma), "x", ncol(x$Sigma),
      "| positive definite:", x$diagnostics$Sigma_pd, "\n")
  cat("  mean|off-diag| / mean diag :",
      round(x$diagnostics$Sigma_offdiag_ratio, 4),
      "  <- discarded by unweighted two-stage analysis\n")
  cat("  condition number :", signif(x$diagnostics$Sigma_condition, 4), "\n")
  invisible(x)
}

#' Empirical semivariogram of model residuals
#'
#' Reports residual spatial autocorrelation before and after the fitted spatial
#' surface is removed. A post-correction curve that is flat to within sampling
#' error is the diagnostic that Module 1 worked; a curve that still rises means
#' spatial structure is being left in the residual.
#'
#' @param d data.frame of plot records.
#' @param response,col,row column names.
#' @param fit a fitted SpATS object (optional).
#' @param maxlag largest lag distance, in plots.
#' @return data.frame with columns \code{lag}, \code{gamma_raw}, \code{gamma_adj}.
#' @examples
#' sim <- sutra_simulate(n_gen = 20, n_env = 1, n_rep = 2, seed = 1)
#' s <- sutra_spatial(sim$data, response = "y", genotype = "genotype",
#'                    col = "col", row = "row")
#' head(s$diagnostics$semivariogram)
#' @export
sutra_semivariogram <- function(d, response, col, row, fit = NULL, maxlag = 12) {
  y <- d[[response]]
  r <- if (!is.null(fit)) fit$residuals else y - mean(y, na.rm = TRUE)
  cc <- d[[col]]; rr <- d[[row]]
  grid_of <- function(v) {
    m <- matrix(NA_real_, nrow = max(rr), ncol = max(cc))
    m[cbind(rr, cc)] <- v
    m
  }
  gr <- grid_of(y - mean(y, na.rm = TRUE)); ga <- grid_of(r)
  sv <- function(M, L) {
    a <- as.vector(M[-seq_len(L), , drop = FALSE] - M[seq_len(nrow(M) - L), , drop = FALSE])
    b <- as.vector(M[, -seq_len(L), drop = FALSE] - M[, seq_len(ncol(M) - L), drop = FALSE])
    v <- c(a, b); v <- v[is.finite(v)]
    if (!length(v)) return(NA_real_)
    0.5 * mean(v^2)
  }
  L <- seq_len(min(maxlag, max(nrow(gr), ncol(gr)) - 1))
  data.frame(lag = L,
             gamma_raw = vapply(L, function(l) sv(gr, l), numeric(1)),
             gamma_adj = vapply(L, function(l) sv(ga, l), numeric(1)))
}
