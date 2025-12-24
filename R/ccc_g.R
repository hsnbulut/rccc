# Robust concordance correlation suggested by King and Chinchilli
robust_ccc<-function(x, y,
                     method = c("absolute","absolute1.5","winsor","huber"),
                     delta  = 1.5, z0 = NULL) {
  method <- match.arg(method)
  if (is.null(z0)) z0 <- 1.5 * sd(x - y)
  gfun <- switch(method,
                 absolute     = function(z) abs(z),
                 absolute1.5  = function(z) abs(z)^delta,
                 winsor       = function(z) ifelse(abs(z) <= z0, z^2, z0^2),
                 huber        = function(z) ifelse(abs(z) <= z0, 0.5*z^2, z0*abs(z) - 0.5*z0^2))
  E_dep <- mean(gfun(x - y))
  E_ind <- mean(outer(x, y, function(a,b) gfun(a - b)))
  rho_g <- 1 - E_dep / E_ind
  max(-1, min(1, rho_g))
}