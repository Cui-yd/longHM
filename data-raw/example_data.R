## code to prepare `example_lmm` and `example_gee` dataset goes here

usethis::use_data(example_lmm, overwrite = TRUE)
usethis::use_data(example_gee, overwrite = TRUE)

library(MASS)

sigmaGenerate = function(visit, sigma){
  S = matrix(sigma[2], nrow = visit, ncol = visit)
  diag(S) = sigma[1]
  return(S)
}

simData = function(n, p, alpha, beta, visit.times, sigma_vec, seed = NULL, method = c("lmm", "gee")) {

  if (is.null(seed)) seed = 123456
  set.seed(seed)


  visit = sample((visit.times[1]):(visit.times[2]), n, replace = T)
  wave = as.vector(unlist(sapply(visit, function(x) (.2*(1:x-1)))))
  id = rep(1:n, times = visit)
  X = rnorm(n, 1, 0.5)


  ck = t(runif(p,0,0.5))
  M = matrix(0, n, p)
  colnames(M) = paste0("M",1:ncol(M))
  Z1 = runif(n, 0, 1)
  Z2 = rbinom(n, 1, 0.4)
  for (i in 1:n) {
    em = rnorm(p, 0, 1)
    M[i,] =  ck + X[i]*alpha + Z1[i]*0.5 + Z2[i]*0.5 + em
  }
  Z1 = rep(Z1, times = visit) + wave*0.5
  Z2 = rep(Z2, times = visit)
  COV = data.frame(Z1, Z2)
  X = rep(X, times = visit)
  M = apply(M, 2, function(x) rep(x, times = visit))

  switch (method,
          "lmm" = {
            e_ij = unlist(sapply(visit, function(x) rnorm(x, 0, 0.5)))
            sigma = matrix(c(sigma_vec[1], sigma_vec[2], sigma_vec[2], sigma_vec[1]), nrow = 2)
            bi = apply(mvrnorm(n, mu=c(0,0),Sigma=sigma), 2, function(x) rep(x, times = visit))
            Y = (0.5+bi[,1]) + 0.8*X + (0.5+bi[,2])*Z1 + 0.5*Z2 + M%*%beta + e_ij
          },
          "gee" = {
            e_ij = list()
            for (i in 1:n){
              sigma = sigmaGenerate(visit = visit[i], sigma = sigma_vec)
              e_ij[[i]] = mvrnorm(1, rep(0, visit[i]), sigma)
            }
            e_ij = unlist(e_ij)
            Y = 0.5 + 0.8*X + 0.5*Z1 + 0.5*Z2 + M%*%beta + e_ij
          }
  )

  return(list(Y = Y, X = X, M = M, COV = COV, id = id, wave = wave, visit = visit, e_ij = e_ij, alpha = alpha, beta = beta))
}


### set parameters

p = 500
n = 200
topN = ceiling(2 * n/log(n))

alpha = beta = rep(0, p)
alpha[1:6] = c(rep(0.6, 4), rep(0.4, 2))
beta[c(1:4, 7:8)] = c(rep(0.6, 4), rep(0.4, 2))
sigma_vec = c(1, 0.2)
visit.times = c(6,9)


### generate data
example_lmm = simData(n, p, alpha, beta, visit.times = visit.times, sigma_vec = sigma_vec,
                      seed = 123456, method = "lmm")

example_gee = simData(n, p, alpha, beta, visit.times = visit.times, sigma_vec = sigma_vec,
                      seed = 123456, method = "gee")

