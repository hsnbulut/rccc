ccc_L1<-function(x, y,
                 equal.means = TRUE,
                 boots       = FALSE) {
  mat <- cbind(x, y)
  fit <- L1pack::l1ccc(mat,
                       equal.means = equal.means,
                       boots       = boots)
  fit$rho1
}