#include <Rcpp.h> //to use the NumericVector object
using namespace Rcpp; //to use the NumericVector object

#include<vector> //to use std::vector<double>

// [[Rcpp::export]]
double cost_rcpp(NumericVector x, String changetype){
  int n = x.size();
  double m = mean(x);
  float c = 0;
  for(int i=0; i<n; ++i){
    c += pow(x[i]- m,2);
  }
  if (changetype == "mean"){
    return c;
  }
  if (changetype == "meanvar"){
    if (n==1){
      return std::numeric_limits<double>::infinity();
    }
    else {return n*(log(c/n)+1);}
  }
}


// [[Rcpp::export]]
List OP_rcpp(NumericVector x, double beta, String changetype){
  int n = x.size();
  IntegerVector cp (n);
  NumericVector Fcost (n+1, -beta);
  for (int i = 1; i <= n; i++){
    NumericVector Fcompare (i);
    double min = std::numeric_limits<double>::infinity();
    for (int j = 0; j < i; j++){
      Fcompare[j] = Fcost[j] + cost_rcpp(x[Range(j, i-1)], changetype) + beta;
      if (Fcompare[j] <= min) {
        min = Fcompare[j];
        cp[i-1] = j;
      }
    }
    Fcost[i] = min;
  }
  IntegerVector cps;
  while (cp[n-1]>0){
    cps.push_back(cp[n-1]);
    n = cp[n-1];
  }
  List result;
  result["cps"] = cps;
  result["Q"] = tail(Fcost, n=1);
  return(result);;
}

  
// [[Rcpp::export]]
List PELT_rcpp(NumericVector x, double beta, String changetype){
  int n = x.size();
  IntegerVector cp (n);
  IntegerVector R (1);
  NumericVector Fcost (n+1, -beta);
  NumericVector nc (n);
  for (int i = 1; i <= n; i++){
    int nR = R.size();
    NumericVector Fcompare (i);
    double min = std::numeric_limits<double>::infinity();
    int nc_min = 0;
    for (int j = 0; j < nR; j++){
      Fcompare[R[j]] = Fcost[R[j]] + cost_rcpp(x[Range(R[j], i-1)], changetype);
      if (Fcompare[R[j]] + beta < min) {
        min = Fcompare[R[j]] + beta;
        cp[i-1] = R[j];
        nc_min = nc[R[j]] + 1;
      }
      else if (Fcompare[R[j]] + beta == min) {
        if (nc[R[j]] + 1 <= nc_min){
        cp[i-1] = R[j];
        nc_min = nc[R[j]] + 1;
        }
      }
    }
    Fcost[i] = min;
    IntegerVector Rnew (1,i);
    for (int k = 0; k < nR; k++){
      if (Fcompare[R[k]] <= Fcost[i]) {
        Rnew.push_back(R[k]);
      }
    }
    R = Rnew;
    if ((changetype == "meanvar")&&(i>1)){
      R[0]=i-1;
      }
  }
  IntegerVector cps;
  while (cp[n-1]>0){
    cps.push_back(cp[n-1]);
    n = cp[n-1];
  }
  List result;
  result["cps"] = cps;
  result["Q"] = tail(Fcost, n=1);
  return(result);
}


