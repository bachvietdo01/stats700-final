#include "global.h"

// Cholesky rank one update and downdate to be used in cpp

arma::mat chol_update_arma(arma::mat LL, arma::vec xx, int D) {
  arma::mat L(LL.memptr(), D, D, true);
  arma::vec x(xx.memptr(), D, true);
  int n = x.size();
  for(int k=0; k<n; k++){
    double r = sqrt(L(k, k)*L(k, k) + x[k]*x[k]);
    double c = r / L(k, k);
    double s = x[k] / L(k, k);
    L(k, k) = r;
    for(int j=k+1; j<n; j++){
      L(k, j) = (L(k, j) + s*x[j])/c;
      x[j] = c*x[j] - s*L(k, j);
    }
  }
  return L;
}

arma::mat chol_downdate(arma::mat LL, arma::vec xx, int D) {
  arma::mat L(LL.memptr(), D, D, true);
  arma::vec x(xx.memptr(), D, true);
  int n = x.size();
  for(int k=0; k<n; k++){
    double r = sqrt(L(k, k)*L(k, k) - x[k]*x[k]);
    if(std::isnan(r)){
      printf("Error: chol downdate problem!\n");
    }
    double c = r / L(k, k);
    double s = x[k] / L(k, k);
    L(k, k) = r;
    for(int j=k+1; j<n; j++){
      L(k, j) = (L(k, j) - s*x[j])/c;
      x[j] = c*x[j] - s*L(k, j);
    }
  }
  return L;
}


// Cholesky rank one update and downdate to be used in R

//' @export
// [[Rcpp::export]]
NumericMatrix chol_update(NumericMatrix LL, NumericVector xx) {
  NumericMatrix L = Rcpp::clone(LL);
  NumericVector x = Rcpp::clone(xx);
  int n = x.size();
  for(int k=0; k<n; k++){
    double r = sqrt(L(k, k)*L(k, k) + x[k]*x[k]);
    double c = r / L(k, k);
    double s = x[k] / L(k, k);
    L(k, k) = r;
    for(int j=k+1; j<n; j++){
      L(k, j) = (L(k, j) + s*x[j])/c;
      x[j] = c*x[j] - s*L(k, j);
    }
  }
  return L;
}

//' @export
// [[Rcpp::export]]
NumericMatrix chol_downdate(NumericMatrix LL, NumericVector xx) {
  NumericMatrix L = Rcpp::clone(LL);
  NumericVector x = Rcpp::clone(xx);
  int n = x.size();
  for(int k=0; k<n; k++){
    double r = sqrt(L(k, k)*L(k, k) - x[k]*x[k]);
    double c = r / L(k, k);
    double s = x[k] / L(k, k);
    L(k, k) = r;
    for(int j=k+1; j<n; j++){
      L(k, j) = (L(k, j) - s*x[j])/c;
      x[j] = c*x[j] - s*L(k, j);
    }
  }
  return L;
}

double logsumexp(double a, double b){
  double m = max(a, b);
  if(arma::is_finite(m)){
    return(log(exp(a-m) + exp(b-m)) + m);
  } else{
    return(m);
  }
}

// Numerical apprixmation of log(V_n(t))

//' @export
// [[Rcpp::export]]
NumericVector calculate_log_V_n(double alpha, int N, int how_many){
  how_many = min(how_many, N);
  NumericVector log_V_n(how_many);
  
  double a, b, c;
  
  for(int t=1; t<=how_many; t++){
    a = 0;
    c = -1 * R_PosInf;
    int k = t;
    while(abs(a-c) > 1e-12){
      a = c;
      b = Rf_lgammafn(t+1) - Rf_lgammafn(k-t+1) - Rf_lgammafn(k*alpha+N) + Rf_lgammafn(k*alpha) - log(exp(1)-1) - Rf_lgammafn(k+1);
      c = logsumexp(a, b);
      k += 1;
    }
    log_V_n[t-1] = c;
  }
  return log_V_n;
}
