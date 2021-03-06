
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

  if(any(is.na(res))) {
    cat("mean =", mean, "\n")
    cat("lower =", lower, "\n")
    cat("upper =", upper, "\n")
    cat("sd =", sd, "\n")
    stop("NA's in res")
  }

  return(res)
}


dtruncdist = function(x, mean, sd, lower, upper) {

  if(lower>mean || upper<mean) stop("invalid mean and lower/upper")
  if(sd<0) stop("sd must be positive")

  if (x < lower || x > upper) {return(0)} else {
  # transform to mean and var of beta dist
  mu = (mean - lower)/(upper-lower)
  v = sd^2 / (upper - lower)^2

  if( v>0.25) stop("Invalid sd on (lower, upper)")

  a = ( (1-mu)/v - 1/mu ) * mu^2
  b = a * (1/mu - 1)
  
  return( ((x - lower)/(upper - lower))^(a-1) * (1-(x - lower)/(upper - lower))^(b-1) / beta(a,b) / (upper - lower) )}
}
