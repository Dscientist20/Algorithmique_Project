## GPL-3 License
## Copyright (c) 2020 Groupe 8

cost <- function(x){
  return(sum((x-mean(x))^2)/2)
}

#' Optimal partitioning algorithm
#'
#' @description Detection by partitioning
#' @param x an initial numeric data
#' @return the change points vector

OP<- function(x, beta = 1){
  n = length(x)
  cp = vector(mode="list")
  cp[[1]] = 0
  F_cost = rep(0, n+1)
  for (i in 2:(n+1)){
    min = Inf
    point = NULL
    for (j in 1:(i-1)){
      if ((F_cost[j] + cost(x[j:(i-1)]) + beta) < min) {
        min = F_cost[j] + cost(x[j:(i-1)]) +beta
        point = j-1 
      }
    }
    F_cost[i] = min
    cp[[i]] = append(cp[[j]], point)
  }
  cps = unique(cp[[n+1]])
  return(cps[-1])
}

#########################################################
#########################################################

#' Pruned Exact Linear Time algorithm
#'
#' @description Detecting by Pruning
#' @param x an initial numeric data
#' @return the change points vector

PELT <- function(x, beta = 1){
  n = length(x)
  cp = vector(mode="list")
  cp[[1]] = 0
  R = vector(mode="list")
  R[[1]] = 0
  F_cost = rep(-beta, n+1)
  for (i in 2:(n+1)){
    Fcompare = rep(0, i-1)
    min = Inf
    point = NULL
    for (j in R[[i-1]]){
      Fcompare[j+1] = F_cost[j+1] + cost(x[(j+1):(i-1)])
      if ((Fcompare[j+1] + beta) <= min) {
        min = Fcompare[j+1] + beta
        point = j 
      }
    }
    F_cost[i] = min
    cp[[i]] = unique(append(cp[[point+1]], point))
    R[[i]] = c(i-1)
    for ( k in R[[i-1]]){
      if (Fcompare[k+1] <= min) {
        R[[i]] = append(R[[i]], k)
      }
    }
  }
  cps = unique(cp[[n+1]])
  return(cps[-1])
}

#########################################################
one.simu <- function(n, type = "sample", algo)
{
  m = sample(n/10)
  if(type == "sample"){v <- sample(m)}else{v <- m:1}
  w = sample(m)
  x=rep(v,w*n/sum(w))+runif(length(rep(v,w*n/sum(w))))
  if (algo == 'OP'){
    t <- system.time(OP(x, 1))[[1]]
  }
  else if (algo == 'PELT'){
    t <- system.time(PELT(x, 1))[[1]]
  }
  return(t)
}
