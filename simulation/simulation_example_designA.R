############################################################
# Sample Simulation Code (Design A) - consistent with manuscript
# Scenarios:
#  - Clean
#  - Unidirectional contamination: contaminate ONLY X by multiplying by 10
#    on the last h = floor(n * eps) observations
#  - Bidirectional contamination: replace m = ceiling(eps * n) rows by
#    Z_out ~ N2( (8, -8)^T, diag(9, 9) )
#
# Estimators (8 total, as in manuscript):
#  rho_C (CCC), rho_R (rCCC), rho_g2, rho_g3, rho_g4, rho_g5, rho_L1, rho_Bayes
#
# RMSE computed w.r.t. true rho_C (for Design A, rho_C = rho)
############################################################

set.seed(20251230)

suppressPackageStartupMessages({
  library(MASS)  # mvrnorm
})

# ----------------------------------------------------------
# 0) Source core functions from repository (adjust path if needed)
# ----------------------------------------------------------

# ----------------------------------------------------------
# 1) Data generation: Design A (bivariate normal, mean 0,0; var 1,1; corr rho)
#    Manuscript: (X, Y)^T ~ N2( (0,0)^T, [[1, rho],[rho,1]] )
# ----------------------------------------------------------
gen_design_A <- function(n, rho) {
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  Z <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
  list(x = Z[, 1], y = Z[, 2])
}

# ----------------------------------------------------------
# 2) Contamination mechanisms (MATCHING THE MANUSCRIPT)
#    Unidirectional: contaminate ONLY X-values of last h=floor(n*eps) obs by *10
#    Bidirectional : replace m=ceiling(eps*n) rows by Z_out ~ N2((8,-8), diag(9,9))
# ----------------------------------------------------------
add_contamination <- function(x, y, eps = 0.10,
                              scenario = c("clean", "unidir", "bidir")) {
  scenario <- match.arg(scenario)
  n <- length(x)
  
  if (scenario == "clean" || eps <= 0) return(list(x = x, y = y))
  
  if (scenario == "unidir") {
    h <- floor(n * eps)
    if (h <= 0) return(list(x = x, y = y))
    idx <- (n - h + 1):n              # LAST h observations (as in manuscript)
    x[idx] <- 10 * x[idx]             # contaminate ONLY X (multiply by 10)
    return(list(x = x, y = y))
  }
  
  if (scenario == "bidir") {
    m <- ceiling(eps * n)
    if (m <= 0) return(list(x = x, y = y))
    idx <- sample.int(n, size = m, replace = FALSE)  # randomly select m rows
    
    mu_out <- c(8, -8)
    Sigma_out <- matrix(c(9, 0, 0, 9), 2, 2)
    Z_out <- MASS::mvrnorm(n = m, mu = mu_out, Sigma = Sigma_out)
    
    x[idx] <- Z_out[, 1]
    y[idx] <- Z_out[, 2]
    return(list(x = x, y = y))
  }
}

# ----------------------------------------------------------
# 3) Estimators (8 total, aligned with manuscript notation)
#    Mapping consistent with your Table 2 ordering:
#      g2 = "absolute"
#      g3 = "absolute1.5"
#      g4 = "winsor"
#      g5 = "huber"
# ----------------------------------------------------------
compute_8_estimators <- function(x, y) {
  c(
    rho_C     = ccc(x, y),
    rho_Rccc     = rccc(x, y),
    rho_g2    = robust_ccc(x, y, method = "absolute"),
    rho_g3    = robust_ccc(x, y, method = "absolute1.5"),
    rho_g4    = robust_ccc(x, y, method = "winsor"),
    rho_g5    = robust_ccc(x, y, method = "huber"),
    rho_L1    = ccc_L1(x, y),
    rho_Bayes = ccc_bayes_mvt(x, y)
  )
}

# ----------------------------------------------------------
# 4) RMSE (w.r.t. true rho_C)
#    For Design A (Table 1), rho_C = rho
# ----------------------------------------------------------
rmse <- function(est, rho_true) sqrt(mean((est - rho_true)^2, na.rm = TRUE))

# ----------------------------------------------------------
# 5) Simulation driver (Design A only)
# ----------------------------------------------------------
run_sim_designA_one_setting <- function(n = 100, rho = 0.9, 
                                        eps = 0.10, 
                                        nrep = 500) {
  
  scenarios <- c("clean", "unidir", "bidir")
  method_names <- names(compute_8_estimators(rnorm(10), rnorm(10)))
  
  # store: scenario -> nrep x 8
  store <- lapply(scenarios, function(.) {
    matrix(NA_real_, nrow = nrep, ncol = length(method_names),
           dimnames = list(NULL, method_names))
  })
  names(store) <- scenarios
  
  for (r in seq_len(nrep)) {
    dat <- gen_design_A(n, rho)
    
    for (s in scenarios) {
      dat_s <- add_contamination(dat$x, dat$y, eps = eps, scenario = s)
      store[[s]][r, ] <- compute_8_estimators(dat_s$x, dat_s$y)
    }
  }
  
  rho_true <- rho  # Design A: true rho_C = rho (Table 1)
  rmse_mat <- sapply(scenarios, function(s) {
    apply(store[[s]], 2, rmse, rho_true = rho_true)
  })
  
  round(rmse_mat, 6)
}

# ----------------------------------------------------------
# 6) Run a sample scenario (editor-friendly example)
# ----------------------------------------------------------
n_sample    <- 100
rho_sample  <- 0.9
eps_sample  <- 0.10
nrep_sample <- 500   # set 3000 to mirror manuscript (will be slower)

cat("=== Sample Simulation Example (Design A, manuscript-consistent) ===\n")
cat("n =", n_sample, "| rho =", rho_sample, "| eps =", eps_sample, "| nrep =", nrep_sample, "\n\n")

rmse_results <- run_sim_designA_one_setting(
  n    = n_sample,
  rho  = rho_sample,
  eps  = eps_sample,
  nrep = nrep_sample
)

print(rmse_results)
