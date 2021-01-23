cost <- function(x, changetype){
  if (changetype=="mean") {
    return(sum((x-mean(x))^2))
  }
  if (changetype=="meanvar"){
    n = length(x)
    if (n==1){
      return(Inf)
    }
    else {return (n*(log(sum((x-mean(x))^2)/n)+1))}
  }
}

OP <- function(x, beta = 0.1, changetype){
  n = length(x)
  cp = rep(0, n)
  F_cost = rep(-beta, n+1)
  for (i in 2:(n+1)){
    Fcompare = rep(0, i-1)
    min = Inf
    point = NULL
    for (j in 1:(i-1)){
      Fcompare[j] = F_cost[j] + cost(x[j:(i-1)],changetype) + beta
      if (Fcompare[j] <= min) {
        min = Fcompare[j]
        point = j-1 
      }
    }
    F_cost[i] = min
    cp[i-1] = point
  }
  cps = vector(mode="numeric")
  while (cp[n]>0){
    cps = append(cps, cp[n])
    n = cp[n]
  }
  return(list(cps=cps, Q=tail(F_cost, 1)))
}



PELT <- function(x, beta = 0.1, changetype){
  n = length(x)
  cp = rep(0, n)
  R = c(0)
  F_cost = rep(-beta, n+1)
  for (i in 2:(n+1)){
    Fcompare = rep(0, i-1)
    min = Inf
    point = NULL
    for (j in R){
      Fcompare[j+1] = F_cost[j+1] + cost(x[(j+1):(i-1)], changetype)
      if ((Fcompare[j+1] + beta) <= min) {
        min = Fcompare[j+1] + beta
        point = j
      }
    }
    F_cost[i] = min
    cp[i-1] = point
    R_new = c(i-1)
    for (k in R){
      if (Fcompare[k+1] <= min) {
        R_new = c(R_new, k)
      }
    }
    R = R_new
    if ((changetype == 'meanvar')&(i>2)){
      R[1] = i-2
      }
  }
  cps = vector(mode="numeric")
  while (cp[n]>0){
    cps = append(cps, cp[n])
    n = cp[n]
  }
  return(list(cps=cps, Q=tail(F_cost, 1)))
}
