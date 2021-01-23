CROPS <- function(x, beta, changetype, max_iter = 50){
  n=length(beta)-1
  stat=rep(1,n)
  stat=c(stat, 0)
  ncp=rep(0, n+1)
  cps=vector(mode='list')
  Q=rep(0, n+1)
  iter=1
  
  beta1 = beta[1]
  beta2 = beta[2]
  cp1 = PELT_rcpp(x, beta1, changetype)
  cp2 = PELT_rcpp(x, beta2, changetype)
  cps[[1]] = cp1$cps
  cps[[2]] = cp2$cps
  ncp[1] = length(cps[[1]])
  ncp[2] = length(cps[[2]])
  Q[1] = cp1$Q
  Q[2] = cp2$Q
  
  while (sum(stat)>0){
    if (iter>max_iter) {break}
    iter=iter+1
    n=length(beta)-1
    for (i in n:1){
      if (stat[i]==1){
        beta1 = beta[i]
        beta2 = beta[i+1]
        Q1 = Q[i]
        Q2 = Q[i+1]
        n1 = ncp[i]
        n2 = ncp[i+1]
        if (n1 > (n2+1)){
          beta_int = (Q1-Q2-n1*beta1+n2*beta2)/(n2-n1)
          int = PELT_rcpp(x, beta_int, changetype)
          cps_int = int$cps
          Q_int = int$Q
          n_int = length(cps_int)
          if ((beta_int == beta2)|(beta_int == beta1)){
            stat[i] = 0
          }
          else{
            beta = append(beta, beta_int, i)
            ncp = append(ncp, n_int, i)
            cps = append(cps, list(cps_int), i)
            Q = append(Q, Q_int , i)
            if (n_int!=n2){
              stat = append(stat, 1, i)
            }
            else {stat = append(stat, 0, i)}
          }
        }
        else {stat[i]=0}
      }
    }
  }
  return(list(beta=beta,ncp=ncp,cps=cps))
}
