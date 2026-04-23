library(glue)

base_path <- "/home/annihong/projects/multiplex-social-networks/multilevelp2/code"
out_file  <- "/home/annihong/projects/multip2/inst/replication_code/applied_multilevel_bully/comparison.txt"

fit1 <- readRDS(file.path(base_path, "bully1_2000_full.Rds"))
fit2 <- readRDS(file.path(base_path, "bully1_2000_full2.Rds"))
fit3 <- readRDS(file.path(base_path, "bully1_2000_full3.Rds"))

fits <- list(full = fit1, full2 = fit2, full3 = fit3)

sink(out_file)

# 1. Sampling settings
cat("=== Sampling settings ===\n")
for (nm in names(fits)) {
  f <- fits[[nm]]
  cat(glue("{nm}: iter={f@sim$iter}, warmup={f@sim$warmup}, chains={length(f@sim$samples)}, saved={(f@sim$iter - f@sim$warmup) * length(f@sim$samples)}\n"))
}

# 2. Model structure inferred from par_dims
cat("\n=== Model structure (inferred from parameter dimensions) ===\n")
for (nm in names(fits)) {
  pd <- fits[[nm]]@par_dims
  T_  <- pd$mu                        # number of layers
  H_  <- pd$cross_mu                  # number of layer pairs
  # sigma_within dim: c(T, 2 + sum(D_within[1,]))
  within_dim <- pd$sigma_within[2]    # 2 + sum(D_within[1,])
  # sigma_cross dim: c(H, 2 + sum(D_cross[1,]))
  cross_dim  <- if (!is.null(pd$sigma_cross)) pd$sigma_cross[2] else NA
  # z dim: c(2*T, sum(n)) -> actor random effects
  n_actors   <- pd$z[2]
  # L inferred from mu_random: c(L, T)
  L_         <- pd$mu_random[1]
  # group covariates: sigma_within actor dim
  cat(glue(
    "{nm}: T={T_}, H={H_}, L={L_}, total_actors={n_actors},",
    " within_dim={within_dim} (2 + sum(D_within[1,])),",
    " cross_dim={cross_dim} (2 + sum(D_cross[1,]))\n"
  ))
}

# 4. Parameter estimates
params_to_compare <- c("mu", "rho", "cross_mu", "cross_rho", "alpha_fixed_coef", "beta_fixed_coef", "mu_group_coef")
cat("\n=== Parameter estimates (mean [2.5%, 97.5%]) ===\n")
for (param in params_to_compare) {
  cat(glue("\n-- {param} --\n"))
  for (nm in names(fits)) {
    s <- summary(fits[[nm]])$summary
    rows <- s[grep(glue("^{param}(\\[|$)"), rownames(s)), c(1, 4, 8), drop = FALSE]
    rows <- round(rows, 3)
    colnames(rows) <- c("mean", "2.5%", "97.5%")
    cat(nm, ":\n"); print(rows)
  }
}



sink()
cat("Output written to", out_file, "\n")
