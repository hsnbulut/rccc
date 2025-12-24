rccc<-function(x, y, alpha = 0.75, tol = 1e-12) {
  
  # Pairwise finite observations
  ok <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  
  # Common sample size after filtering
  n_common <- length(x)
  if (n_common < 3) return(NA_real_)
  
  # Degenerate variance check
  if (var(x) < tol || var(y) < tol) return(NA_real_)
  
  # Robust covariance (MCD)
  mcd <- tryCatch(
    CovMcd(data.frame(x, y), alpha = alpha),
    error = function(e) NULL
  )
  if (is.null(mcd)) return(NA_real_)
  
  mu <- mcd$center
  S  <- mcd$cov
  
  denom <- S[1,1] + S[2,2] + (mu[1] - mu[2])^2
  
  # Numerical stability checks
  if (!is.finite(denom) || denom < tol) return(NA_real_)
  if (!is.finite(S[1,2])) return(NA_real_)
  
  coef <- (2 * S[1,2]) / denom
  if (!is.finite(coef)) return(NA_real_)
  
  coef
}
