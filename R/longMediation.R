#' MAIN FUNCTION
#'
#' @param Y repeated measured outcomes for n individuals.
#' @param X repeated measured exposure / treatment / input.
#' @param M time-invariant repeated measure of M high-dimensional mediators.
#' @param COV covariants, can be null.
#' @param id the id of all the individuals.
#' @param wave is used to indicate the measurement occasion in repeated measures,
#' representing which time point a sample was measured.
#' @param topN the number of top variables selected during the SIS 
#' (Sure Independence Screening) step. If it is null, use default topN = 2*n/log(n), 
#' n is the total sample size. It can be either a single numeric value or a vector 
#' of two numeric values. If it is a single value, the top topN variables are 
#' selected in both the X->M and M->Y steps. If it is a vector of two values, 
#' the first value (topN1) indicates the number of top variables selected in the
#' X->M margin, and the second value (topN2) indicates the number of top variables 
#' selected in the M->Y margin
#' @param SIS_by specify the marginal correlation method for conducting the SIS. 
#' The user can choose between two options: "pval" (sort by p-value) or 
#' "coef" (sort by absolute value of coefficient).
#' @param method method used in this framework. Choices include 
#' "lmm" (linear mixed effect model),
#' "gee" (generalized estimating equation).
#' @param verbose logical: whether to provides additional detailed information 
#' or messages during the execution of this function (default = FALSE).
#'
#'
#' @return A list with the following members:
#' @return \code{result}: a data frame of main results.
#' sob_pval_fdr (FDR-BH adjusted Sobel test p-value).
#' sob_pval_bon (Bonferroni adjusted Sobel test p-value).
#' sob_pval_by (FDR-BY adjusted Sobel test p-value).
#' join_pval_fdr (FDR-BH adjusted Joint significant test p-value).
#' join_pval_bon (Bonferroni adjusted Joint significant test p-value).
#' join_pval_by (FDR-BY adjusted Joint significant test p-value).
#' ab_coef (Estimated indirect effect, alpha * beta).
#' ab_var (Variance of estimated indirect effect).
#' conf_low (Lower confidence bound).
#' conf_up (Upper confidence bound).
#' sob_pval(Sobel test p-value before multiple hypothesis adjustment).
#' alpha_coef (Estimated alpha coefficient).
#' alpha_var (Estimated alpha variance).
#' alpha_pval(Estimated alpha p-value).
#' beta_coef(Estimated beta coefficient).
#' beta_var (Estimated beta variance).
#' beta_pval (Estimated beta  p-value).
#' @return \code{subID_SIS}: ID of the intersect SIS selected mediators.
#' @return \code{subM_SIS}: a matrix of the selected mediators in step 1.
#' @return \code{Beta.margin}: marginal correlation matrix in M->Y margin.
#' @return \code{Alpha.margin}: marginal correlation matrix in X->M margin.
#' @return \code{elapsed}: elapsed time, units = mins.
#'
#'
#' @import lmerTest
#' @import geepack
#' @importFrom utils head
#' @importFrom stats gaussian p.adjust pnorm qnorm
#'
#' @export
#'
#' @examples
#' # Example data for demonstration, with 200 individuals and 500 mediators:
#'
#' data(example_lmm)
#' n = length(unique(example_lmm$id))
#' topN = ceiling(2 * n/log(n))
#' results_lmm = longMediation(Y = example_lmm$Y, X = example_lmm$X,
#' M = example_lmm$M, COV = example_lmm$COV, id = example_lmm$id,
#' wave = example_lmm$wave, topN = topN, SIS_by = "coef",
#' method = "lmm", verbose = TRUE)
#'
#' data(example_gee)
#' n = length(unique(example_gee$id))
#' topN = ceiling(2 * n/log(n))
#' results_gee = longMediation(Y = example_gee$Y, X = example_gee$X,
#' M = example_gee$M, COV = example_gee$COV, id = example_gee$id,
#' wave = example_gee$wave, topN = topN, SIS_by = "coef",
#' method = "gee", verbose = TRUE)



longMediation = function(Y, X, M, COV = NULL, id, wave, topN = NULL,
                         SIS_by = c("coef","pval"), method = c("lmm", "gee"),
                         verbose = F) {


  begin = Sys.time()

  ##### FIRST STEP: INTERSECT SIS #####

  if (verbose) message("Step 1: Screening...", "     (", Sys.time(), ")")

  n = length(unique(id))

  if (is.null(topN)) {
    topN = 2*ceiling(n/log(n))
  } else if (length(topN) == 1) {
    topNAlpha = topNBeta = topN
  } else {
    topNAlpha = topN[1]
    topNBeta = topN[2]
  }
  if (verbose) print(paste0("Select top ", topNAlpha, " at X->M pathway and top ", topNBeta, " at M->Y pathway"))

  Beta.margin = longBetaMarginal(Y = Y, M = M, X = X, COV = COV, id = id, wave = wave, method = method)
  Alpha.margin = longAlphaMarginal(X = X, M = M, COV = COV, wave = wave)
  switch(SIS_by,
         "coef" = {
           SIS.beta.rank = sort(abs(Beta.margin[1, ]), decreasing = T)
           SIS.alpha.rank = sort(abs(Alpha.margin[1, ]), decreasing = T)
         },
         "pval" = {
           SIS.beta.rank = sort(Beta.margin[3, ], decreasing = F)
           SIS.alpha.rank = sort(Alpha.margin[3, ], decreasing = F)
         })

  top.beta = names(head(SIS.beta.rank, topNBeta))
  top.alpha = names(head(SIS.alpha.rank, topNAlpha))
  subID_SIS = intersect(top.beta, top.alpha)

  if (length(subID_SIS) == 0) {
    s = 0
    while (s >= 0) {
      s = s + 1
      topN = s * ceiling(n/log(n))
      top.beta = names(head(SIS.beta.rank, topN))
      top.alpha = names(head(SIS.alpha.rank, topN))
      subID_SIS = intersect(top.beta, top.alpha)
      if (length(subID_SIS) > 0) {
        break
      }
    }
  }

  subM_SIS = data.frame(M[, subID_SIS])
  colnames(subM_SIS) = subID_SIS
  subdat_SIS = data.frame(id, Y, X, COV, subM_SIS)

  if(verbose) message("        Top ", length(subID_SIS), " mediators selected (by marginal correlation): ", paste(subID_SIS, sep = " "), "\n")



  ##### SECOND STEP: INDIRECT EFFECT TEST #####

  if(verbose) message("Step 2: Indirect effect testing ...", "     (", Sys.time(), ")")

  # beta estimation
  switch (method,
          "lmm" = {
            lmm.fit = suppressWarnings(lmer(Y ~ . -id + (1+Z1|id), data = subdat_SIS))
            lmm.coef = as.data.frame(summary(lmm.fit)$coef)
            lmm.M = lmm.coef[grep("M", rownames(lmm.coef), value=TRUE),]
            beta_coef = lmm.M[,1]
            beta_var = lmm.M[,2]^2
            beta_pval = lmm.M[,5]
          },
          "gee" = {
            gee.fit = suppressWarnings(geeglm(Y ~ . -id, id=id, data=subdat_SIS, corstr="exchangeable", family = gaussian))
            gee.coef = summary(gee.fit)$coef
            gee.M = gee.coef[grep("M", rownames(gee.coef), value=TRUE),]
            beta_coef = gee.M[, 1]
            beta_var = gee.M[, 2]^2
            beta_pval = gee.M[, 4]
          }
  )


  # alpha estimation
  alpha.fit = longAlphaMarginal(X = X, M = subM_SIS, COV = COV, wave = wave)
  alpha_coef = alpha.fit[1, subID_SIS]
  alpha_var = alpha.fit[2, subID_SIS]
  alpha_pval = alpha.fit[3, subID_SIS]


  # p value calculation
  ab_coef = alpha_coef * beta_coef
  ab_var = (alpha_coef^2) * (beta_var) + (beta_coef^2) * (alpha_var)
  # confidence interval
  conf_low = ab_coef - qnorm(0.975)*sqrt(ab_var);  conf_up = ab_coef + qnorm(0.975)*sqrt(ab_var)


  # sobel test for alpha and beta
  s.test = (abs(ab_coef))/(sqrt(ab_var))
  sob_pval = 2 * (1 - pnorm(s.test))
  sob_pval_fdr = p.adjust(sob_pval, "fdr", length(subID_SIS))
  sob_pval_bon = p.adjust(sob_pval, "bonferroni", length(subID_SIS))
  sob_pval_by = p.adjust(sob_pval, "BY", length(subID_SIS))

  if(verbose){
    message("        Significant FDR(BH) Sobel mediator(s): ", names(sob_pval_fdr)[which(sob_pval_fdr < 0.05)])
    message("        Significant bonforroni Sobel mediator(s): ", names(sob_pval_bon)[which(sob_pval_bon < 0.05)])
    message("        Significant FDR(BY) Sobel mediator(s): ", names(sob_pval_by)[which(sob_pval_by < 0.05)], "\n")
  }

  # joint test for alpha and beta
  alpha_pval_fdr = p.adjust(alpha_pval, "fdr", length(subID_SIS))
  beta_pval_fdr = p.adjust(beta_pval, "fdr", length(subID_SIS))
  pval_bind_fdr = rbind(alpha_pval_fdr, beta_pval_fdr)
  join_pval_fdr = apply(pval_bind_fdr, 2, max)

  alpha_pval_bon = p.adjust(alpha_pval, "bonferroni", length(subID_SIS))
  beta_pval_bon = p.adjust(beta_pval, "bonferroni", length(subID_SIS))
  pval_bind_bon = rbind(alpha_pval_bon, beta_pval_bon)
  join_pval_bon = apply(pval_bind_bon, 2, max)

  alpha_pval_by = p.adjust(alpha_pval, "BY", length(subID_SIS))
  beta_pval_by = p.adjust(beta_pval, "BY", length(subID_SIS))
  pval_bind_by = rbind(alpha_pval_by, beta_pval_by)
  join_pval_by = apply(pval_bind_by, 2, max)

  elapsed = difftime(Sys.time(), begin, units = "mins")

  if(verbose) {
    message("        Significant FDR(BH) Joint mediator(s): ", names(join_pval_fdr)[which(join_pval_fdr < 0.05)])
    message("        Significant bonforroni Joint mediator(s): ", names(join_pval_bon)[which(join_pval_bon < 0.05)])
    message("        Significant FDR(BY) Joint mediator(s): ", names(join_pval_by)[which(join_pval_by < 0.05)])
  }

  result = data.frame(sob_pval_fdr, sob_pval_bon, sob_pval_by,
                      join_pval_fdr, join_pval_bon, join_pval_by,
                      ab_coef, ab_var, conf_low, conf_up, sob_pval,
                      alpha_coef, alpha_var, alpha_pval,
                      beta_coef, beta_var, beta_pval)

  return(list(result = result, subID_SIS = subID_SIS, subM_SIS = subM_SIS,
              Beta.margin = Beta.margin, Alpha.margin = Alpha.margin, elapsed = elapsed))
}

