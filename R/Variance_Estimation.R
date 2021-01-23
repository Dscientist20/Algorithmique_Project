hall <- function(x, order = 3){
  n=length(x)
  d1 = c(0.5^0.5,-0.5^0.5)
  d2 = c(0.8090, -0.5,-0.3090) 
  d3 = c(0.1942, 0.2809, 0.3832, -0.8582)
  d4 = c(0.2708, -0.0142, 0.6909, -0.4858, -0.4617)
  d5 = c(0.9064, -0.2600, -0.2167, -0.1774, -0.1420, -0.1103)
  d = get(paste("d",order,sep=""))
  s=0
  for (k in 1:(n-order)){
    s = s + sum((d*x[k:(k+order)]))^2
  }
  return(sqrt(s/(n-order)))
}

