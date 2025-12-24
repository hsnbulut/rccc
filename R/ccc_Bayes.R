ccc_bayes_mvt<-function(x, y, conf.level = 0.95, nmcmc = 2000) {
  ratings <- cbind(x, y)
  fit <- robust_bayesian_ccc(
    ratings,
    num_observers = 2,
    n_mcmc   = nmcmc,
    n_burnin = floor(nmcmc / 2),
    nu_fixed = 10
  )
  as.numeric(fit$point_estimate)
}

robust_bayesian_ccc<-function(Y,                # n x d matrix
                              num_observers,    # d
                              n_mcmc   = 5000,
                              n_burnin = 1000,
                              nu_fixed = 10,    # nu sabit tutuluyor
                              mu_0     = NULL,
                              Sigma_0_scale = 1000,
                              rho_prior = NULL,
                              V_prior   = NULL,
                              seed      = NULL) {
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  n <- nrow(Y)
  d <- ncol(Y)
  
  if (d != num_observers) {
    stop("num_observers ile Y'nin sütun sayısı uyumlu değil.")
  }
  if (n < 3) stop("En az 3 gözlem olmalı.")
  
  if (!is.null(seed)) set.seed(seed)
  
  # Önceler
  if (is.null(mu_0))    mu_0    <- rep(0, d)
  Sigma_0 <- diag(Sigma_0_scale, d)
  Sigma_0_inv <- solve(Sigma_0)
  
  rho <- if (is.null(rho_prior)) d else rho_prior
  V   <- if (is.null(V_prior)) diag(1, d) else V_prior
  V_inv <- solve(V)
  
  # MCMC saklama
  mu_samples    <- matrix(NA_real_, n_mcmc, d)
  Sigma_samples <- array(NA_real_, dim = c(d, d, n_mcmc))
  lambda_samples <- matrix(NA_real_, n_mcmc, n)
  ccc_samples   <- numeric(n_mcmc)
  
  # Başlangıç
  mu_current    <- colMeans(Y)
  Sigma_current <- cov(Y)
  lambda_current <- rep(1, n)
  nu_current    <- nu_fixed
  
  for (t in seq_len(n_mcmc)) {
    # A. lambda güncelleme
    Sigma_inv_current <- solve(Sigma_current)
    for (i in seq_len(n)) {
      yi_minus_mu <- Y[i, ] - mu_current
      quad <- drop(t(yi_minus_mu) %*% Sigma_inv_current %*% yi_minus_mu)
      shape <- (nu_current + d) / 2
      rate  <- (nu_current + quad) / 2
      lambda_current[i] <- rgamma(1, shape = shape, rate = rate)
    }
    lambda_samples[t, ] <- lambda_current
    
    # B. mu güncelleme
    sum_lambda <- sum(lambda_current)
    sum_lambda_Y <- colSums(lambda_current * Y)
    
    A <- Sigma_inv_current * sum_lambda + Sigma_0_inv
    b <- Sigma_inv_current %*% sum_lambda_Y + Sigma_0_inv %*% mu_0
    A_inv <- solve(A)
    mu_current <- MASS::mvrnorm(1, mu = A_inv %*% b, Sigma = A_inv)
    mu_samples[t, ] <- mu_current
    
    # C. Sigma^-1 ~ Wishart güncelleme
    S_lambda <- matrix(0, d, d)
    for (i in seq_len(n)) {
      yi_minus_mu <- Y[i, ] - mu_current
      S_lambda <- S_lambda + lambda_current[i] * (yi_minus_mu %*% t(yi_minus_mu))
    }
    
    Wishart_V_inv <- V_inv + S_lambda
    Wishart_V     <- solve(Wishart_V_inv)
    Wishart_nu    <- rho + n
    
    Sigma_inv_new <- stats::rWishart(1, df = Wishart_nu, Sigma = Wishart_V)[,,1]
    Sigma_current <- solve(Sigma_inv_new)
    Sigma_samples[,,t] <- Sigma_current
    
    # D. nu sabit tutuluyor (nu_fixed)
    # (Makaledeki tam Bayes modeli için ayrıca nu güncellemesi gerekir.)
    
    # E. CCC hesapla (Feng formülüne göre basit d=2 durumu)
    nu_term <- nu_current / (nu_current - 2)
    
    # d=2 olduğundan:
    # Numerator: 2 * nu/(nu-2) * sigma_12
    numerator <- 2 * nu_term * Sigma_current[1, 2]
    
    # Denominator:
    # (d-1) * nu/(nu-2) * (sigma_11 + sigma_22) + (mu_1 - mu_2)^2
    term1 <- (d - 1) * nu_term * (Sigma_current[1,1] + Sigma_current[2,2])
    term2 <- (mu_current[1] - mu_current[2])^2
    denominator <- term1 + term2
    
    ccc_samples[t] <- numerator / denominator
  }
  
  # Burn-in sonrası örnekler
  keep <- (n_burnin + 1):n_mcmc
  ccc_post <- ccc_samples[keep]
  
  point_estimate <- stats::median(ccc_post)
  cred_int <- stats::quantile(ccc_post, probs = c(0.025, 0.975), na.rm = TRUE)
  
  list(
    point_estimate   = point_estimate,
    credible_interval = cred_int,
    ccc_samples      = ccc_post,
    mu_samples       = mu_samples[keep, , drop = FALSE],
    Sigma_samples    = Sigma_samples[,, keep, drop = FALSE]
  )
}