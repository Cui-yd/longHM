#' Marginal correlation estimation for M->Y
#'
#' @param Y repeated measured outcomes for n individuals.
#' @param X repeated measured exposure / treatment / input.
#' @param M time-invariant repeated measure of M high-dimensional mediators.
#' @param COV covariants, can be null.
#' @param id the id of all the individuals.
#' @param wave is used to indicate the measurement occasion in repeated measures,
#' representing which time point a sample was measured.
#' @param method method used in this framework. Choices include "lmm" (linear mixed effect model),
#' "gee" (generalized estimating equation).
#'
#' @return \code{Est_beta}: a matrix of the main results of M->Y margin.
#' First row is the estimated coefficients,
#' Second row is the estimated variance.
#' Third row is the estimated p-value.
#'
#' @import lmerTest
#' @import geepack
#'
#' @export

longBetaMarginal = function(Y, X, M, COV, id, wave, method = c("lmm", "gee")) {

  p = ncol(M)
  Est_beta = matrix(0, 3, p)
  colnames(Est_beta) = colnames(M)

  for(j in 1:p) {
    MY = data.frame(M = M[, j], Y = Y, X = X, COV = COV, id = id, wave = wave)
    switch (method,
            "lmm" = {
              lmm.fit = suppressWarnings(lmer(Y ~ . -id + (1+COV.Z1|id), data = MY))
              lmm.coef = summary(lmm.fit)$coef
              Est_beta[1,j] = lmm.coef[2,1]                   # coefficients of beta
              Est_beta[2,j] = lmm.coef[2,2]^2                 # variance of beta
              Est_beta[3,j] = lmm.coef[2,5]                   # p-value of beta
            },
            "gee" = {
              gee.fit = suppressWarnings(geeglm(formula = Y ~ . -id, family = gaussian, data = MY, id = id, corstr="exchangeable"))
              gee.coef = summary(gee.fit)$coef
              Est_beta[1,j] = gee.coef[2,1]                   # coefficients of beta
              Est_beta[2,j] = gee.coef[2,2]^2                 # variance of beta
              Est_beta[3,j] = gee.coef[2,4]                   # p-value of beta
            }
    )
  }
  return(Est_beta = Est_beta)
}
