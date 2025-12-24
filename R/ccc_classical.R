# classical concordance correlation coefficient 
ccc<-function(x, y) {
  n  <- length(x)
  mx <- mean(x); my <- mean(y)
  sxx <- sum((x - mx)^2) / n
  syy <- sum((y - my)^2) / n
  sxy <- sum((x - mx) * (y - my)) / n
  2 * sxy / (sxx + syy + (mx - my)^2)
}