# Simulation for change in mean:
mean.simu <- function(n, m)
{
  changetype = 'mean'
  a = floor((n-30)/60)
  v = sample(1:a, m, replace = FALSE) - runif(m)/2
  cps = sort(round(v*60, 0))
  ngroup = diff(c(0,cps,n))
  x=unlist(lapply(ngroup, function(x) rnorm(x,rnorm(1,0,2.5),1)))
  return(list(x=x,cps=cps))
}


# Simulation for change in both mean and variance:
var.simu <- function(n, m)
{
  changetype = 'meanvar'
  a = floor((n-30)/60)
  v = sample(1:a, m, replace = FALSE) - runif(m)/2
  cps = sort(round(v*60, 0))
  ngroup = diff(c(0,cps,n))
  x=unlist(lapply(ngroup, function(x) rnorm(x,0,rlnorm(1,0,log(10)/2))))
  return(list(x=x,cps=cps))
}

