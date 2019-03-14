#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma ;


// [[Rcpp::export]]
Rcpp::List MGSPsamp(int p, int n, int k,
                    double as, double bs, double df,
                    double ad1, double bd1,
                    double ad2, double bd2,
                    double adf, double bdf,
                    double b0, double b1,
                    int sp, int nrun, int burn, int thin,
                    double prop, double epsilon,
                    arma::vec ps, arma::mat Sigma,
                    arma::mat Lambda, arma::mat meta,
                    arma::mat veta, arma::mat psijh,
                    arma::vec theta, arma::vec tauh,
                    arma::mat Plam,arma::mat Y,
                    arma::mat scaleMat, Rcpp::StringVector output){
  Rcpp::StringVector cm = "covMean";
  Rcpp::StringVector cv = "covSamples";
  Rcpp::StringVector fs = "factSamples";
  Rcpp::StringVector ss = "sigSamples";
  Rcpp::StringVector nf = "numFactors";
  
  bool covm = Rcpp::any(Rcpp::in(cm, output)).is_true();
  bool cov = Rcpp::any(Rcpp::in(cv, output)).is_true();
  bool fac = Rcpp::any(Rcpp::in(fs, output)).is_true();
  bool sig = Rcpp::any(Rcpp::in(ss, output)).is_true();
  bool nfa = Rcpp::any(Rcpp::in(nf, output)).is_true();
  mat COVMEAN;
  cube OMEGA;
  Rcpp::List LAMBDA;
  cube SIGMA(p, p, sp, fill::zeros);
  vec K(sp, fill::zeros);
  if(covm) mat COVMEAN(p, p, fill::zeros);
  if(cov) cube OMEGA(p, p, sp, fill::zeros);
  if(fac) Rcpp::List LAMBDA;
  if(sig) cube SIGMA(p, p, sp, fill::zeros);
  if(nfa) vec K(sp, fill::zeros);
  int ind = 0;
  
  for(int i=0; i<nrun; ++i){
    // UPDATE ETA
    mat Lmsg = Lambda.each_col() % ps;
    mat Veta1 = eye<mat>(k,k) + Lmsg.t() * Lambda;
    mat Tmat = chol(Veta1);
    mat Q, R;
    qr(Q, R, Tmat);
    mat S = inv(R);
    mat Veta = S * S.t();
    mat Meta = Y * Lmsg * Veta;
    mat noise(n, k, fill::randn);
    mat eta = Meta + noise * S.t();
    
    // UPDATE LAMBDA
    mat eta2 = eta.t() * eta;    // prepare eta crossproduct before the loop
    for(int j = 0; j < p; ++j) {
      mat Llamt = chol(diagmat(Plam.row(j)) + ps(j)*eta2);
      Lambda.row(j) = solve(Llamt, randn<vec>(k) + 
        solve(Llamt, solve(Llamt.t(), ps(j) * eta.t() * Y.col(j)))).t();
    }
    
    // UPDATE psihj
    mat Lambda_sq = pow(Lambda,2);
    mat shape = Lambda_sq.each_row() % tauh.t();
    for (int i=0; i < p; i++) {
      for (int j=0; j < k; j++) {
        psijh(i,j) = arma::randg(distr_param(df/2 + 0.5,1/shape(i,j)));
        
      }
    }
    
    // UPDATE THETA & TAUH
    mat matr = psijh % pow(Lambda, 2);
    double ad = ad1 + 0.5 * p * k;
    double bd = bd1 + 0.5 * theta[0] * sum(tauh.t() % sum(matr, 0));
    theta[0] = randg(distr_param(ad, 1 / bd));           
    tauh = cumprod(theta);
    
    for(int h=1; h<k; ++h) {
      double ad = ad2 + 0.5 * p *(k-h);
      vec tauh_sub = tauh.subvec(h,k-1);
      double bd = bd2 + 0.5 * sum(tauh_sub.t() % sum(matr.cols(h,k-1), 0));
      theta(h) = randg(distr_param(ad, 1 / bd)); 
      tauh = cumprod(theta);
    }
    
    // UPDATE SIGMA
    mat Ytil = Y - eta * Lambda.t();
    mat bsvec =  bs + 0.5*sum(pow(Ytil,2), 0);
    for(int i = 0; i < p ; ++i){
      ps(i) = randg(distr_param(as + 0.5*n, 1 / bsvec(0,i))); 
    } 
    Sigma = diagmat(1/ps);
    
    //UPDATE PLAM
    Plam = psijh.each_row() % tauh.t();
    
    //ADAPT K
    double prob = 1 / std::exp(b0 + b1 * (i + 1));            // probability of adapting
    double uu = randu();
    umat llog = abs(Lambda) < epsilon;
    mat lint = arma::conv_to<arma::mat>::from(llog);
    mat lind = sum(lint, 0) / p;  // proportion of elements in each column less than eps in magnitude
    umat vecs = lind.row(0) >= prop;
    int num = 0;
    for(int h = 0; h < k ; ++h){
      num += vecs(0,h);
    } 
    umat vecs2 = lind < prop;
    
    vec lind2 = arma::conv_to<arma::vec>::from(lind);
    if(uu < prob) {
      if((i > 20) && (num == 0) && all(lind2 < 0.995)) {
        k = k + 1;
        Lambda.col(k-1) = 0;
        eta.col(k-1) = randn(n);
        psijh.col(k-1) = randg(p,distr_param(df/2,df/2));
        theta(k-1) = randg(distr_param(ad2,bd2));
        tauh = cumprod(theta);
        Plam = psijh.each_row() % tauh.t();
      } else {
        if (num > 0) {
          ivec facts = {k - num,1};
          k = max(facts);
          Lambda = Lambda.cols(vecs2);
          psijh = psijh.cols(vecs2);
          eta = eta.cols(vecs2);
          theta = theta.elem(vecs2);
          tauh = cumprod(theta);
          Plam = psijh.each_row() % tauh.t();
        }
      }
    }
    bool thincheck = i - std::floor(i/thin) * thin; // % operator stolen by arma
    if(~thincheck && (i > burn)) {
      mat Omega = (Lambda * Lambda.t() + Sigma) * scaleMat;
      if(covm) COVMEAN += Omega / sp;
      if(cov) OMEGA.slice(ind) = Omega;
      if(fac) LAMBDA[ind] = Lambda;
      if(sig) SIGMA.slice(ind) = Sigma;
      if(nfa) K(ind) = k;
      ind += 1;
    }
  }
  return Rcpp::List::create(Rcpp::Named("covMean") = COVMEAN,
                            Rcpp::Named("covSamps") = OMEGA,
                            Rcpp::Named("factSamps") = LAMBDA,
                            Rcpp::Named("sigSamps") = SIGMA,
                            Rcpp::Named("numFact") = K);
}

