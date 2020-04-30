
# Draws from the truncated and scaled beta distribution 

rtruncdist = function(nsamples, mean, sd, lower, upper) {

  if(lower>mean || upper<mean) stop("invalid mean and lower/upper")
  if(sd<0) stop("sd must be positive")

  # transform to mean and var of beta dist
  mu = (mean - lower)/(upper-lower)
  v = sd^2 / (upper - lower)^2

  if( v>0.25) stop("Invalid sd on (lower, upper)")

  a = ( (1-mu)/v - 1/mu ) * mu^2
  b = a * (1/mu - 1)
  res = lower + (upper - lower)*rbeta(nsamples, shape1=a, shape2=b)
  return(res)
}

