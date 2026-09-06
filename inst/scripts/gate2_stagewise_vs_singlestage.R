## GATE 2 -------------------------------------------------------------------
## Does the spectrally-rotated stage-wise fit reproduce a genuine single-stage
## analysis of the raw plot data?  Simulated with NO smooth spatial surface so
## that a true single-stage model is fittable in lme4 and the comparison is
## like-for-like.
suppressMessages({library(sutra); library(lme4)})

sim <- sutra_simulate(n_gen = 40, n_env = 5, n_rep = 2,
                      var_g = 1.0, var_ge = 0.4, var_spatial = 0,
                      var_row = 0.15, var_col = 0.12, var_e = 0.8,
                      crossover = 0, seed = 42)
dat <- sim$data
dat$env <- factor(dat$env); dat$genotype <- factor(dat$genotype)
dat$rf <- factor(paste(dat$env, dat$row)); dat$cf <- factor(paste(dat$env, dat$col))

## ---- gold standard: single-stage LMM on raw plots -------------------------
ss <- lmer(y ~ env + (1 | genotype) + (1 | rf) + (1 | cf), data = dat, REML = TRUE)
b_ss <- ranef(ss)$genotype[, 1]
names(b_ss) <- rownames(ranef(ss)$genotype)

## ---- stage one, per environment ------------------------------------------
s1 <- lapply(levels(dat$env), function(e) {
  sutra_spatial(subset(dat, env == e), response = "y", genotype = "genotype",
                col = "col", row = "row", env = e, genotype.as.random = FALSE)
})

## ---- stage two, three weightings -----------------------------------------
cmp <- sutra_compare_weighting(s1)
eff <- cmp$effects
b_ss <- b_ss[eff$genotype]

res <- data.frame(
  weighting = c("fully_efficient", "diagonal", "unweighted"),
  cor_with_single_stage = sapply(c("fully_efficient","diagonal","unweighted"),
                                 function(w) cor(eff[[w]], b_ss)),
  spearman = sapply(c("fully_efficient","diagonal","unweighted"),
                    function(w) cor(eff[[w]], b_ss, method = "spearman")),
  rmse = sapply(c("fully_efficient","diagonal","unweighted"),
                function(w) sqrt(mean((scale(eff[[w]]) - scale(b_ss))^2))),
  row.names = NULL)

truth <- sim$truth$g[eff$genotype]
res$cor_with_TRUTH <- sapply(c("fully_efficient","diagonal","unweighted"),
                             function(w) cor(eff[[w]], truth))
res$single_stage_cor_with_TRUTH <- cor(b_ss, truth)

cat("\n=== GATE 2: stage-wise vs single-stage (n_gen=40, n_env=5) ===\n")
print(res, digits = 5)
cat("\nvariance components by weighting:\n")
print(cmp$summary, digits = 5)

pass <- res$cor_with_single_stage[res$weighting == "fully_efficient"] > 0.99
cat("\nGATE 2:", ifelse(pass, "PASS", "FAIL"),
    "(fully efficient vs single stage r > 0.99)\n")
saveRDS(list(res = res, summary = cmp$summary), "/home/claude/gate2.rds")
