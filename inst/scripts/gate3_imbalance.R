## GATE 3 -------------------------------------------------------------------
## The literature's claim is NOT mainly about the accuracy of genotype point
## estimates.  Piepho et al. (2012) and Damesa et al. (2019) claim that ignoring
## the stage-one variance structure BIASES THE GENETIC VARIANCE and hence
## heritability.  So the primary metric here is bias in sigma^2_g against a known
## truth of 1.0, with accuracy reported alongside.
suppressMessages(library(sutra))
TRUE_VG <- 1.0

run_scenario <- function(label, var_e, dropout, rep_range, n_sim = 30,
                         n_gen = 60, n_env = 6) {
  rows <- list()
  for (i in seq_len(n_sim)) {
    sim <- sutra_simulate(n_gen = n_gen, n_env = n_env, n_rep = 2,
                          var_g = TRUE_VG, var_ge = 0.4, var_spatial = 1.5,
                          var_row = 0.15, var_col = 0.12, var_e = var_e,
                          dropout = dropout, rep_range = rep_range,
                          crossover = 0, seed = 1000 + i)
    d <- sim$data
    s1 <- lapply(unique(d$env), function(e)
      try(sutra_spatial(d[d$env == e, ], response = "y", genotype = "genotype",
                        col = "col", row = "row", env = e,
                        genotype.as.random = FALSE), silent = TRUE))
    if (any(vapply(s1, inherits, logical(1), "try-error"))) next
    cmp <- try(sutra_compare_weighting(s1), silent = TRUE)
    if (inherits(cmp, "try-error")) next
    truth <- sim$truth$g[cmp$effects$genotype]
    for (w in c("fully_efficient", "diagonal", "unweighted")) {
      rows[[length(rows) + 1]] <- data.frame(
        sim = i, weighting = w,
        acc = cor(cmp$effects[[w]], truth),
        vg  = cmp$summary$sigma2_g[cmp$summary$weighting == w])
    }
  }
  a <- do.call(rbind, rows)
  agg <- do.call(rbind, lapply(split(a, a$weighting), function(z)
    data.frame(weighting = z$weighting[1], n = nrow(z),
               acc = mean(z$acc),
               vg_mean = mean(z$vg),
               vg_bias_pct = 100 * (mean(z$vg) - TRUE_VG) / TRUE_VG)))
  agg$scenario <- label
  agg[, c("scenario", "weighting", "n", "acc", "vg_mean", "vg_bias_pct")]
}

scen <- rbind(
  run_scenario("A balanced, homogeneous (null)", 0.8, 0, NULL),
  run_scenario("B balanced, heterogeneous var", c(.2,.4,.8,1.5,2.5,4), 0, NULL),
  run_scenario("C unequal rep, homogeneous",    0.8, 0, c(1,4)),
  run_scenario("D unequal rep + heterogeneous", c(.2,.4,.8,1.5,2.5,4), 0, c(1,4)),
  run_scenario("E D + 30% dropout",             c(.2,.4,.8,1.5,2.5,4), .30, c(1,4)),
  run_scenario("F extreme: 40x var spread, 1-5 reps, 40% dropout",
               c(.15,.3,.8,2,4,6), .40, c(1,5))
)
cat("\n=== GATE 3 ===  true sigma^2_g = 1.0 ; bias% is what the paper claims\n\n")
print(scen, digits = 4, row.names = FALSE)
saveRDS(scen, "/home/claude/gate3.rds")
