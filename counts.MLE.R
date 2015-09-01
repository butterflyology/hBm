Counts.MLE = function() {

negLL = function(param) {
  pA = rep(expit(param[1]), nsites)
  pB = rep(expit(param[2]), nsites)
  beta = param[-(1:2)]
  phi = xMat %*% beta
  lambda = as.numeric(exp(phi))
  piMat = cbind(pA * pB, pA * (1-pB), (1-pA) * pB)
  logLike = dpois(yMat, lambda * piMat, log = TRUE)
  (-1) * sum(logLike)
}

expit = function(z) { 1 / (1 + exp(-z)) }
logit = function(z) { log(z) - log(1-z) }
  
# read data from file
  df = read.table('counts.txt', header = TRUE)
  yMat = as.matrix(df[, -1])
  nsites = dim(yMat)[1]
  hab = df[, 'habitat']
  xMat = model.matrix(~hab - 1)

# fit model to data
  ySum = apply(yMat, 1, sum)
  pGuess = 0.5
  lambdaGuess = mean(ySum) / pGuess
  paramGuess = c(rep(logit(pGuess), 2), rep(log(lambdaGuess), 3))
  fit = optim(par = paramGuess, fn = negLL, method = 'BFGS', hessian = TRUE)
  
  list(logLikelihood = -fit$value, mle = fit$par, covMat = chol2inv(chol(fit$hessian)))
}

