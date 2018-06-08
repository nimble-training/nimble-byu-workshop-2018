## @knitr dbetabin

dbetabin <- nimbleFunction(
    run = function(x = double(0), alpha = double(0), beta = double(0), size = double(0), 
        log = integer(0, default = 0)) {
        
        returnType(double(0))
        logProb <- lgamma(size+1) - lgamma(x+1) - lgamma(size - x + 1) +
            lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta) +
            lgamma(x + alpha) + lgamma(size - x + beta) - lgamma(size + alpha + beta)
        if(log) return(logProb)
        else return(exp(logProb))
    })

rbetabin <- nimbleFunction(
    run = function(n = integer(0), alpha = double(0), beta = double(0), size = double(0)) {
        returnType(double(0))
        if(n != 1) print("rbetabin only allows n = 1; using n = 1.")
        p <- rbeta(1, alpha, beta)
        return(rbinom(1, size = size, prob = p))
    })

## @knitr littersMarg-code

littersMargCode <- nimbleCode({
  for (i in 1:G) {
     for (j in 1:N) {
     	 # (marginal) likelihood (data model)
        r[i,j] ~ dbetabin(a[i], b[i], n[i,j])
     }
     # prior for hyperparameters
     a[i] ~ dgamma(1, .001)
     b[i] ~ dgamma(1, .001)
   }
})

## @knitr littersMarg-model

littersMargModel <- nimbleModel(littersMargCode, 
          data = littersData, constants = littersConsts, inits = littersInits)

cLittersMargModel <- compileNimble(littersMargModel)

## @knitr sv-code

stochVolCode <- nimbleCode({
  x[1] ~ dnorm(phi * x0, sd = sigma)
  y[1] ~ dnorm(0, var = betaSquared * exp(x[1]))
  for(t in 2:T){
        x[t] ~ dnorm(phi * x[t-1], sd = sigma)
        y[t] ~ dnorm(0, var = betaSquared * exp(x[t]))
  }
  x0 ~ dnorm(1, sd = sigma)
  phi <- 2 * phiStar - 1
  phiStar ~ dbeta(18, 1)
  sigma ~ T(dt(mu = 0, sigma = 1, df = 1), 0, )
  betaSquared <- beta^2
  beta ~ T(dt(mu = 0, sigma = 1, df = 1), 0, )
})

## @knitr sv-model

library('stochvol')
data('exrates')
y <- 100 * logret(exrates$USD[exrates$date > '2012-02-01'])
stochVolModel <- nimbleModel(code = stochVolCode,
   constants = list(T = 44), data = list(y = y),
   inits = list(beta = .5992, phi = .9702,
   sigma = .178, x0 = 0))
CstochVolModel <- compileNimble(stochVolModel)
