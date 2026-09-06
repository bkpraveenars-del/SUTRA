## GATE 3 (definitive) -------------------------------------------------------
## Two metrics, two GxE settings, five imbalance scenarios.
##  var_ge = 0   isolates the weighting question (no GxE leakage into sigma^2_g)
##  var_ge = 0.4 is the realistic case; expected leakage is var_ge / n_env
suppressMessages(library(sutra))
N_SIM <- 25; N_ENV <- 6

run <- function(label, var_e, dropout, rep_range, var_ge) {
  rows <- list()
  for (i in seq_len(N_SIM)) {
    sim <- sutra_simulate(n_gen = 60, n_env = N_ENV, n_rep = 2, var_g = 1,
                          var_ge = var_ge, var_spatial = 1.5, var_row = .15,
                          var_col = .12, var_e = var_e, dropout = dropout,
                          rep_range = rep_range, crossover = 0, seed = 5000 + i)
    d <- sim$data
    s1 <- lapply(unique(d$env), function(e)
      try(sutra_spatial(d[d$env==e,], response="y", genotype="genotype",
                        col="col", row="row", env=e, genotype.as.random=FALSE),
          silent = TRUE))
    if (any(vapply(s1, inherits, logical(1), "try-error"))) next
    cmp <- try(sutra_compare_weighting(s1), silent = TRUE)
    if (inherits(cmp, "try-error")) next
    tr <- sim$truth$g[cmp$effects$genotype]
    for (w in c("fully_efficient","diagonal","unweighted"))
      rows[[length(rows)+1]] <- data.frame(w = w,
        acc = cor(cmp$effects[[w]], tr),
        vg  = cmp$summary$sigma2_g[cmp$summary$weighting == w])
  }
  a <- do.call(rbind, rows)
  o <- do.call(rbind, lapply(split(a, a$w), function(z) data.frame(
    weighting = z$w[1], n = nrow(z),
    acc = mean(z$acc), acc_se = sd(z$acc)/sqrt(nrow(z)),
    vg = mean(z$vg), vg_bias_pct = 100*(mean(z$vg) - 1))))
  o$scenario <- label; o$var_ge <- var_ge
  o$expected_leak_pct <- 100 * var_ge / N_ENV
  o[, c("var_ge","scenario","weighting","n","acc","acc_se","vg","vg_bias_pct","expected_leak_pct")]
}

scen <- list(
  list("A balanced, homogeneous",        0.8,                    0,   NULL),
  list("B balanced, heterogeneous",      c(.2,.4,.8,1.5,2.5,4),  0,   NULL),
  list("C unequal rep, homogeneous",     0.8,                    0,   c(1,4)),
  list("D unequal rep + heterogeneous",  c(.2,.4,.8,1.5,2.5,4),  0,   c(1,4)),
  list("E D + 30% dropout",              c(.2,.4,.8,1.5,2.5,4), .30,  c(1,4)))

out <- list()
for (vge in c(0, 0.4))
  for (s in scen)
    out[[length(out)+1]] <- run(s[[1]], s[[2]], s[[3]], s[[4]], vge)
res <- do.call(rbind, out)
saveRDS(res, "/home/claude/gate3_definitive.rds")
cat("\n=== GATE 3 DEFINITIVE (true sigma^2_g = 1.0) ===\n\n")
for (v in unique(res$var_ge)) {
  cat("---- var_ge =", v, " (expected sigma^2_g leakage:",
      unique(res$expected_leak_pct[res$var_ge==v]), "%) ----\n")
  print(res[res$var_ge==v, c("scenario","weighting","acc","acc_se","vg_bias_pct")],
        digits = 4, row.names = FALSE)
  cat("\n")
}
