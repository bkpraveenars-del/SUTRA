#' Module 2: fully efficient stage-wise linear mixed model
#'
#' Takes the per-environment output of \code{\link{sutra_spatial}} and fits the
#' across-environment model
#' \deqn{\hat{y} = X\beta + Zu + e,\qquad \mathrm{Var}(e) = \Sigma}
#' where \eqn{\Sigma} is the block-diagonal collection of the FULL stage-one
#' error variance-covariance matrices, not their diagonals and not the identity.
#'
#' Tractability comes from whitening rather than from approximation. Each
#' environment's block is spectrally decomposed, \eqn{\Sigma_j = U\Lambda U^{\top}},
#' and the stage-one means and their design matrices are pre-multiplied by
#' \eqn{\Lambda^{-1/2}U^{\top}}. On the rotated scale the errors are independent
#' with unit variance, so an ordinary mixed-model solver returns the fully
#' efficient weighted solution (Piepho et al. 2012).
#'
#' Three weighting schemes are provided so that the cost of the usual shortcut
#' can be measured rather than assumed:
#' \describe{
#'   \item{\code{fully_efficient}}{the full \eqn{\Sigma} (this framework).}
#'   \item{\code{diagonal}}{only \eqn{\mathrm{diag}(\Sigma)}; the common
#'         "weighted two-stage" compromise.}
#'   \item{\code{unweighted}}{\eqn{\Sigma = I}; what exporting adjusted means
#'         to a spreadsheet actually assumes.}
#' }
#'
#' @param stage1 a list of \code{sutra_stage1} objects, one per environment.
#' @param weighting one of \code{"fully_efficient"}, \code{"diagonal"},
#'   \code{"unweighted"}.
#' @param genotype.as.random logical; TRUE returns genotype BLUPs (default),
#'   FALSE returns across-environment BLUEs.
#' @param contrast_space logical; project each environment's variance matrix
#'   onto the genotype-contrast space before whitening, removing the intercept
#'   variance. Keep TRUE. Setting FALSE reproduces the naive rotation and is
#'   provided only to demonstrate why it fails.
#' @param shrink numeric in [0, 1]; linear shrinkage of the stage-one variance
#'   matrix toward its own diagonal, \eqn{(1-\alpha)\Sigma + \alpha\,
#'   \mathrm{diag}(\Sigma)}. \code{0} uses the full matrix, \code{1} is
#'   equivalent to \code{weighting = "diagonal"}. Whitening by an ESTIMATED
#'   covariance amplifies that matrix's own estimation error along its smallest
#'   eigenvalues, so an ill-conditioned block can make the fully efficient
#'   analysis markedly worse than an unweighted one. Shrinkage is the standard
#'   remedy and is strongly recommended whenever any block has a condition
#'   number above roughly 30; see \code{\link{sutra_condition}}.
#' @param max_condition numeric; cap the condition number of each block by
#'   flooring its eigenvalues at \code{max(ev) / max_condition}. An alternative
#'   to \code{shrink}; \code{Inf} disables it.
#' @param ridge small value added to the eigenvalues before inversion, guarding
#'   against a numerically singular block.
#' @return An object of class \code{sutra_stage2}.
#' @examples
#' sim <- sutra_simulate(n_gen = 20, n_env = 3, n_rep = 2, seed = 2)
#' s1 <- lapply(unique(sim$data$env), function(e)
#'   sutra_spatial(sim$data[sim$data$env == e, ], response = "y",
#'                 genotype = "genotype", col = "col", row = "row", env = e))
#' fit <- sutra_stagewise(s1, weighting = "fully_efficient")
#' fit
#' @export
sutra_stagewise <- function(stage1,
                            weighting = c("fully_efficient", "diagonal", "unweighted"),
                            genotype.as.random = TRUE,
                            contrast_space = TRUE,
                            shrink = 0,
                            max_condition = Inf,
                            ridge = 1e-8) {
  weighting <- match.arg(weighting)
  if (inherits(stage1, "sutra_stage1")) stage1 <- list(stage1)
  if (!all(vapply(stage1, inherits, logical(1), "sutra_stage1")))
    stop("stage1 must be a list of sutra_stage1 objects")

  gen_levels <- sort(unique(unlist(lapply(stage1, function(s) s$blues$genotype))))
  env_levels <- vapply(stage1, function(s) s$env, character(1))
  if (anyDuplicated(env_levels)) stop("environment labels must be unique")

  yl <- Zl <- list()
  diag_scale <- numeric(0)

  for (j in seq_along(stage1)) {
    s  <- stage1[[j]]
    yj <- s$blues$value
    gj <- factor(s$blues$genotype, levels = gen_levels)
    nj <- length(yj)
    if (nj < 2) next

    Zj <- stats::model.matrix(~ gj - 1)
    colnames(Zj) <- gen_levels

    Wj <- switch(weighting,
      fully_efficient = s$Sigma,
      diagonal        = diag(diag(s$Sigma), nrow = nj),
      unweighted      = diag(1, nrow = nj))
    Wj <- (Wj + t(Wj)) / 2
    if (shrink > 0 && weighting == "fully_efficient") {
      Wj <- (1 - shrink) * Wj + shrink * diag(diag(Wj), nrow = nj)
    }

    # --- CONTRAST PROJECTION -------------------------------------------------
    # Within one environment the genotype means are identified only up to a
    # common shift, which the environment main effect absorbs.  The variance of
    # that shift (the intercept variance) is therefore NOT part of the
    # information about genotype differences.  It is also, empirically, the
    # dominant eigenvalue of the raw vcov -- its eigenvector is the constant
    # vector -- so whitening by the raw matrix crushes the mean direction by a
    # large factor and distorts the rotated genotype design.
    #
    # Projecting onto the contrast space P = I - J/n removes it.  The projected
    # matrix has rank n-1; the null direction is dropped, and no environment
    # column is then needed in X because centring has already absorbed it.
    if (contrast_space) {
      P  <- diag(nj) - matrix(1 / nj, nj, nj)
      Wc <- P %*% Wj %*% P
      e  <- eigen((Wc + t(Wc)) / 2, symmetric = TRUE)
      k  <- nj - 1L                                   # drop the null direction
      ev <- e$values[seq_len(k)]
      if (any(ev <= 0)) {
        bad <- sum(ev <= 0)
        ev  <- pmax(ev, ridge * max(e$values))
        warning("environment ", s$env, ": ", bad,
                " non-positive eigenvalue(s) in the contrast space were floored")
      }
      if (is.finite(max_condition))
        ev <- pmax(ev, max(ev) / max_condition)
      U  <- e$vectors[, seq_len(k), drop = FALSE]
      A  <- t(U) / sqrt(ev)
    } else {
      e  <- eigen(Wj, symmetric = TRUE)
      ev <- pmax(e$values, ridge * max(e$values))
      A  <- t(e$vectors) / sqrt(ev)
    }

    yl[[length(yl) + 1L]] <- A %*% yj
    Zl[[length(Zl) + 1L]] <- A %*% Zj
    diag_scale <- c(diag_scale, mean(diag(Wj)))
  }

  y <- as.numeric(do.call(rbind, yl))
  Z <- do.call(rbind, Zl)
  X <- matrix(0, nrow = length(y), ncol = 0)   # env effects absorbed by centring
  if (!contrast_space) {
    X <- matrix(1, nrow = length(y), ncol = 1, dimnames = list(NULL, "(Intercept)"))
  }

  fit <- .sutra_fit_lmm(y, X, Z, gen_levels,
                        genotype.as.random = genotype.as.random)

  out <- list(
    weighting = weighting,
    genotype.as.random = genotype.as.random,
    environments = env_levels,
    genotypes = gen_levels,
    effects = fit$effects,          # data.frame(genotype, value, se)
    varcomp = fit$varcomp,
    contrast_space = contrast_space,
    shrink = shrink, max_condition = max_condition,
    residual_scale = fit$sigma2,    # ~1 + attributable GxE when weighting is right
    n_rotated = length(y),
    fit = fit$obj)
  class(out) <- "sutra_stage2"
  out
}

# internal: fit the rotated model
.sutra_fit_lmm <- function(y, X, Z, gen_levels, genotype.as.random = TRUE) {
  if (genotype.as.random) {
    # y = X beta + Z u + e ,  u ~ N(0, sigma_g^2 I) , e ~ N(0, sigma_e^2 I)
    # fit via lme4 with an explicit numeric design: use lmer with a dummy
    # grouping factor is not possible for a rotated Z, so solve REML directly.
    est <- .sutra_reml_ridge(y, X, Z)
    eff <- data.frame(genotype = gen_levels,
                      value = as.numeric(est$u),
                      se    = as.numeric(est$u_se),
                      stringsAsFactors = FALSE)
    list(effects = eff,
         varcomp = c(genotype = est$sigma2_g, residual = est$sigma2_e),
         sigma2  = est$sigma2_e, obj = est)
  } else {
    XX <- if (ncol(X)) cbind(X, Z[, -1, drop = FALSE]) else
          cbind(1, Z[, -1, drop = FALSE])
    qrf <- qr(XX)
    b <- qr.coef(qrf, y); b[is.na(b)] <- 0
    r <- y - XX %*% b
    s2 <- sum(r^2) / (length(y) - qrf$rank)
    XtXinv <- tryCatch(chol2inv(qr.R(qrf)), error = function(e) MASS::ginv(crossprod(XX)))
    se <- sqrt(pmax(diag(XtXinv) * s2, 0))
    ng <- length(gen_levels)
    p0 <- max(ncol(X), 1L)
    val <- c(0, b[(p0 + 1):(p0 + ng - 1)])
    ses <- c(NA_real_, se[(p0 + 1):(p0 + ng - 1)])
    eff <- data.frame(genotype = gen_levels, value = val, se = ses,
                      stringsAsFactors = FALSE)
    list(effects = eff, varcomp = c(residual = s2), sigma2 = s2, obj = qrf)
  }
}

# internal: REML for y = X beta + Z u + e  with u ~ N(0, s2g I), e ~ N(0, s2e I)
# parameterised by lambda = s2e / s2g and profiled, so one-dimensional and stable.
.sutra_reml_ridge <- function(y, X, Z, interval = c(1e-6, 1e6)) {
  n <- length(y); p <- ncol(X); q <- ncol(Z)
  ZtZ <- crossprod(Z); Zty <- crossprod(Z, y); yty <- sum(y^2)
  if (p > 0) {
    XtX <- crossprod(X); XtZ <- crossprod(X, Z); Xty <- crossprod(X, y)
    build <- function(lam) rbind(cbind(XtX, XtZ), cbind(t(XtZ), ZtZ + lam * diag(q)))
    rhs   <- rbind(Xty, Zty)
  } else {
    build <- function(lam) ZtZ + lam * diag(q)
    rhs   <- Zty
  }

  nll <- function(loglam) {
    lam <- exp(loglam)
    C <- build(lam)
    ch <- tryCatch(chol(C), error = function(e) NULL)
    if (is.null(ch)) return(1e10)
    sol <- backsolve(ch, backsolve(ch, rhs, transpose = TRUE))
    rss <- yty - sum(sol * rhs)
    if (rss <= 0) return(1e10)
    logdetC <- 2 * sum(log(diag(ch)))
    # REML objective (up to an additive constant)
    0.5 * ((n - p) * log(rss) + logdetC - q * log(lam))
  }
  opt <- stats::optimize(nll, interval = log(interval))
  lam <- exp(opt$minimum)

  C <- build(lam)
  ch <- chol(C)
  sol <- backsolve(ch, backsolve(ch, rhs, transpose = TRUE))
  rss <- yty - sum(sol * rhs)
  s2e <- as.numeric(rss / (n - p))
  s2g <- s2e / lam
  Cinv <- chol2inv(ch)
  u  <- sol[(p + 1):(p + q)]
  us <- sqrt(pmax(diag(Cinv)[(p + 1):(p + q)] * s2e, 0))
  list(u = u, u_se = us, beta = sol[seq_len(p)],
       sigma2_e = s2e, sigma2_g = s2g, lambda = lam, converged = TRUE)
}

#' @export
print.sutra_stage2 <- function(x, ...) {
  cat("<sutra Module 2: stage-wise integration>\n")
  cat("  weighting        :", x$weighting, "\n")
  cat("  environments     :", length(x$environments),
      "| genotypes:", length(x$genotypes), "\n")
  cat("  variance comps   :",
      paste(names(x$varcomp), signif(x$varcomp, 4), sep = " = ", collapse = "  "), "\n")
  if (x$weighting == "fully_efficient") {
    cat("  residual scale   :", signif(x$residual_scale, 4),
        "(should be ~1; a large departure means the stage-1 variance",
        "structure is wrong)\n")
  }
  invisible(x)
}

#' Compare weighting schemes on the same stage-one input
#'
#' Runs \code{\link{sutra_stagewise}} under all three weighting schemes and
#' returns their genotype effects side by side. This is the measurement that
#' turns "unweighted two-stage analysis is biased" from an assertion into a
#' number for your own data.
#'
#' @param stage1 list of \code{sutra_stage1} objects.
#' @param ... passed to \code{sutra_stagewise}.
#' @return list with \code{effects} (data.frame) and \code{summary}.
#' @examples
#' sim <- sutra_simulate(n_gen = 20, n_env = 3, n_rep = 2, seed = 2)
#' s1 <- lapply(unique(sim$data$env), function(e)
#'   sutra_spatial(sim$data[sim$data$env == e, ], response = "y",
#'                 genotype = "genotype", col = "col", row = "row", env = e))
#' # what is the weighting actually worth on these data?
#' sutra_compare_weighting(s1)$summary
#' @export
sutra_compare_weighting <- function(stage1, ...) {
  schemes <- c("fully_efficient", "diagonal", "unweighted")
  fits <- lapply(schemes, function(w) sutra_stagewise(stage1, weighting = w, ...))
  names(fits) <- schemes
  eff <- data.frame(genotype = fits[[1]]$effects$genotype, stringsAsFactors = FALSE)
  for (w in schemes) eff[[w]] <- fits[[w]]$effects$value
  vc <- do.call(rbind, lapply(schemes, function(w) {
    v <- fits[[w]]$varcomp
    data.frame(weighting = w,
               sigma2_g = unname(v["genotype"]),
               sigma2_e = unname(v["residual"]),
               residual_scale = fits[[w]]$residual_scale,
               stringsAsFactors = FALSE)
  }))
  rk <- sapply(schemes, function(w) eff[[w]])
  vc$cor_with_fully_efficient <- round(
    apply(rk, 2, function(v) stats::cor(v, rk[, "fully_efficient"])), 5)
  vc$rank_cor_with_fully_efficient <- round(
    apply(rk, 2, function(v) stats::cor(v, rk[, "fully_efficient"], method = "spearman")), 5)
  list(effects = eff, summary = vc, fits = fits)
}


#' Condition numbers of the stage-one variance matrices
#'
#' Whitening by an estimated covariance matrix amplifies that matrix's own
#' estimation error along its smallest eigenvalues. An ill-conditioned block can
#' therefore make a fully efficient analysis WORSE than an unweighted one. This
#' function reports, per environment, the condition number of the contrast-space
#' variance matrix that \code{\link{sutra_stagewise}} would actually invert, so
#' the risk can be seen before the model is fitted rather than inferred from a
#' disappointing result afterwards.
#'
#' As a working rule, a condition number above about 30 in any environment is a
#' reason to set \code{shrink} above zero.
#'
#' @param stage1 a list of \code{sutra_stage1} objects.
#' @return data.frame with one row per environment.
#' @examples
#' sim <- sutra_simulate(n_gen = 20, n_env = 2, n_rep = 2, seed = 3)
#' s1 <- lapply(unique(sim$data$env), function(e)
#'   sutra_spatial(sim$data[sim$data$env == e, ], response = "y",
#'                 genotype = "genotype", col = "col", row = "row", env = e))
#' # a condition number above about 30 means shrinkage is advisable
#' sutra_condition(s1)
#' @export
sutra_condition <- function(stage1) {
  if (inherits(stage1, "sutra_stage1")) stage1 <- list(stage1)
  do.call(rbind, lapply(stage1, function(s) {
    n <- nrow(s$Sigma)
    P <- diag(n) - matrix(1 / n, n, n)
    ev <- eigen((P %*% s$Sigma %*% P + t(P %*% s$Sigma %*% P)) / 2,
                symmetric = TRUE, only.values = TRUE)$values[seq_len(n - 1)]
    off <- s$Sigma[upper.tri(s$Sigma)]
    cJ <- mean(off); V0 <- s$Sigma - cJ
    data.frame(env = s$env, n_genotypes = n,
               condition = max(ev) / min(ev),
               offdiag_corr_after_projection =
                 mean(abs(V0[upper.tri(V0)])) / mean(diag(V0)),
               risk = ifelse(max(ev) / min(ev) > 30, "shrink recommended", "ok"),
               stringsAsFactors = FALSE)
  }))
}
