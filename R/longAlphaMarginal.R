#' Marginal correlation estimation for X->M
#'
#' @param X repeated measured exposure / treatment / input.
#' @param M time-invariant repeated measure of M high-dimensional mediators.
#' @param COV covariants, can be null.
#' @param wave is used to indicate the measurement occasion in repeated measures,
#' representing which time point a sample was measured.
#'
#' @return \code{Est_alpha}: a matrix of the main results of X->M margin.
#' First row is the estimated coefficients,
#' Second row is the estimated variance.
#' Third row is the estimated p-value.
#'
#' @export

longAlphaMarginal = function(X, M, COV, wave) {

  p = ncol(M)
  Est_alpha = matrix(0, 3, p)
  colnames(Est_alpha) = colnames(M)

  X = X[which(wave == min(wave))]
  M = data.frame(M[which(wave == min(wave)),])
  COV = COV[which(wave == min(wave)),]
  for(j in 1:p) {
    MX = data.frame(M = M[,j], X = X, COV = COV)
    fit = stats::lm(M ~ ., data = MX)
    Est_alpha[1,j] = summary(fit)$coef[2,1]          # coefficients of alpha
    Est_alpha[2,j] = summary(fit)$coef[2,2]^2        # variance of alpha
    Est_alpha[3,j] = summary(fit)$coef[2,4]          # p-value of alpha
  }
  return(Est_alpha = Est_alpha)
}
