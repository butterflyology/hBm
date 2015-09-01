# Code based on Royale and Dorazio 2008
# Just some pretty straightforward code to get thinking about how to handle count data
# MLE function in counts.MLE.R
Counts.Bayes = function() {

# read data from file
  df = read.table('counts.txt', header = TRUE)
  yMat = as.matrix(df[,-1])
  hab = df[, 'habitat']
  xMat = model.matrix(~hab-1)


# model specification in WinBUGS
  modelFilename = "countsBayes.txt"
cat('
model {

for (k in 1:ncovs) {
    beta[k] ~ dnorm(0.0, 1.0E-2)
    lambda.mean[k] <- exp(beta[k])
}

pA ~ dunif(0,1)
pB ~ dunif(0,1)

pi.1 <- pA*pB
pi.2 <- pA*(1-pB)
pi.3 <- (1-pA)*pB

for (i in 1:n) {

    log(lambda[i]) <- inprod(x[i,], beta[])
    mu.1[i] <- lambda[i]*pi.1
    mu.2[i] <- lambda[i]*pi.2
    mu.3[i] <- lambda[i]*pi.3

    y[i,1] ~ dpois(mu.1[i])
    y[i,2] ~ dpois(mu.2[i])
    y[i,3] ~ dpois(mu.3[i])
}
}
', fill = TRUE, file = modelFilename)



# arguments for bugs()
  nsites = dim(yMat)[1]
  ncovs = dim(xMat)[2]

  data = list(n=nsites, ncovs = ncovs, y = yMat, x = xMat)
  params = list('beta', 'pA', 'pB')

  inits = function() {
    lambdaGuess = mean(apply(yMat,1,sum))
    list(pA = runif(1, 0, 1), pB = runif(1, 0, 1), beta = log(rep(lambdaGuess, ncovs)) )
  }

  
# call to bugs()
  library(R2WinBUGS)
  fit = bugs(data, inits, params,  model.file = modelFilename,
    n.chains=1, n.iter=30000, n.burnin=10000, n.thin=5, bugs.seed=sample(1:9999,size = 1), debug = FALSE, DIC = FALSE)

  fit
}

