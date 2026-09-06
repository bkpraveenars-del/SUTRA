## GATE 4 -------------------------------------------------------------------
## Precision weighting concentrates information on the most precisely measured
## environments.  That is correct when environments differ in precision but
## agree in ranking.  When they disagree -- crossover GxE -- it narrows the
## effective environmental sample.  How much crossover destroys the benefit?
suppressMessages(library(sutra))
N_SIM <- 20
run <- function(crossover, var_e, label) {
  g <- gd <- numeric(0); cond <- numeric(0)
  for (i in seq_len(N_SIM)) {
    sim <- sutra_simulate(n_gen = 50, n_env = 6, n_rep = 2, var_g = 1,
                          var_ge = 0, var_spatial = 1.5, var_e = var_e,
                          rep_range = c(1, 3), crossover = crossover,
                          seed = 900 + i)
    s1 <- lapply(unique(sim$data$env), function(e)
      try(sutra_spatial(sim$data[sim$data$env == e, ], response = "y",
            genotype = "genotype", col = "col", row = "row", env = e), silent = TRUE))
    if (any(vapply(s1, inherits, logical(1), "try-error"))) next
    cmp <- try(sutra_compare_weighting(s1), silent = TRUE)
    if (inherits(cmp, "try-error")) next
    tr <- sim$truth$g[cmp$effects$genotype]
    au <- cor(cmp$effects$unweighted, tr)
    g  <- c(g,  cor(cmp$effects$fully_efficient, tr) - au)
    gd <- c(gd, cor(cmp$effects$diagonal, tr)        - au)
    cond <- c(cond, max(sutra_condition(s1)$condition))
  }
  data.frame(scenario = label, crossover = crossover, n = length(g),
             gain_fullSigma = mean(g), se_full = sd(g)/sqrt(length(g)),
             gain_diagonal  = mean(gd),
             mean_condition = mean(cond))
}
het <- c(.2,.4,.8,1.5,2.5,4)
res <- rbind(
  run(0.0, het, "heterogeneous variance"),
  run(0.3, het, "heterogeneous variance"),
  run(0.6, het, "heterogeneous variance"),
  run(1.2, het, "heterogeneous variance"),
  run(0.0, 0.8, "homogeneous variance"),
  run(0.6, 0.8, "homogeneous variance"))
cat("\n=== GATE 4: gain over UNWEIGHTED as crossover GxE increases ===\n\n")
print(res, digits = 4, row.names = FALSE)
saveRDS(res, "/home/claude/gate4.rds")
