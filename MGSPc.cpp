#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;


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
                    arma::mat scaleMat, Rcpp::StringVector output, 
                    int start, bool verbose){
  // --- initialise output objects --- //
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
  cube OMEGA, SIGMA;
  field<mat> LAMBDA;
  vec K(sp);
  
  if(covm) COVMEAN = zeros<mat>(p, p);
  if(cov) OMEGA = zeros<cube>(p, p, sp);
  if(sig) SIGMA = zeros<cube>(p, p, sp);
  if(fac) LAMBDA = field<mat>(sp);
  if(nfa) vec K = zeros<vec>(sp);
  int ind = 0;
  
  // --- initialise loop objects --- //
  mat Lmsg, Veta1, S, Veta, Meta,
  noise, eta, eta2, Llamt, Llam,Lambda_sq,
  shape, matr, Ytil, bsvec, lint, lind, Omega;
  
  umat llog, vecs;
  
  vec tauh_sub, lind2;
  
  ivec facts;
  
  double ad, bd, prob, uu;
  
  int num;
  
  bool thincheck, printcheck;
  
  // --- loop --- //
  for(int i=0; i<nrun; i++, start++){
    // --- UPDATE ETA --- //
    Lmsg = Lambda.each_col() % ps;
    Veta1 = eye<mat>(k,k) + Lmsg.t() * Lambda;
    S = inv(trimatu(chol(Veta1)));
    Veta = S * S.t();
    Meta = Y * Lmsg * Veta;
    noise = randn(n, k);
    eta = Meta + noise * S.t();
    
    // --- UPDATE LAMBDA --- //
    eta2 = eta.t() * eta;    // prepare eta crossproduct before the loop
    for(int j = 0; j < p; j++) {
      Llamt = trimatu(chol(diagmat(Plam.row(j)) + ps(j)*eta2));
      Llam = trimatl(Llamt.t());
      Lambda.row(j) = (solve(Llamt, randn<vec>(k)) +
        solve(Llamt, solve(Llam, ps(j) * eta.t() * Y.col(j)))).t();
    }
    
    // --- UPDATE psihj --- //
    Lambda_sq = square(Lambda);
    shape = Lambda_sq.each_row() % tauh.t();
    for (int l=0; l < p; l++) {
      for (int j=0; j < k; j++) {
        psijh(l,j) = randg(distr_param(df/2 + 0.5, 1 / (df/2 + shape(l,j))));
        
      }
    }
    
    // --- UPDATE THETA & TAUH --- //
    matr = psijh % square(Lambda);
    ad = ad1 + 0.5 * p * k;
    bd = bd1 + 0.5 / theta(0) * sum(tauh.t() % sum(matr, 0));
    theta(0) = randg(distr_param(ad, 1 / bd));
    tauh = cumprod(theta);
    
    for(int h=1; h<k; h++) {
      ad = ad2 + 0.5 * p *(k-h);
      tauh_sub = tauh.subvec(h,k-1);
      bd = bd2 + 0.5 / theta(h) * sum(tauh_sub.t() % sum(matr.cols(h,k-1), 0));
      theta(h) = randg(distr_param(ad, 1 / bd));
      tauh = cumprod(theta);
    }
    
    // --- UPDATE SIGMA --- //
    Ytil = Y - eta * Lambda.t();
    bsvec =  bs + 0.5 * sum(square(Ytil), 0);
    for(int l = 0; l < p ; l++){
      ps(l) = randg(distr_param(as + 0.5*n, 1 / bsvec(0,l)));
    }
    Sigma = diagmat(1/ps);
    
    // --- UPDATE PLAM --- //
    Plam = psijh.each_row() % tauh.t();
    
    // --- ADAPT K --- //
    prob = 1 / std::exp(b0 + b1 * (start + 1));            // probability of adapting
    uu = randu();
    llog = abs(Lambda) < epsilon;
    lint = arma::conv_to<arma::mat>::from(llog);
    lind = sum(lint, 0) / p;  // proportion of elements in each column less than eps in magnitude
    vecs = lind.row(0) >= prop;
    num = 0;
    for(int h = 0; h < k ; h++){
      num += vecs(0,h);
    }
    lind2 = arma::conv_to<arma::vec>::from(lind);
    if(uu < prob) {
      if((i > 20) && (num == 0) && all(lind2 < 0.995)) {
        k = k + 1;
        Lambda.resize(p, k);
        Lambda.col(k-1) = zeros<vec>(p);
        eta.resize(n, k);
        eta.col(k-1) = randn(n);
        psijh.resize(p, k);
        psijh.col(k-1) = randg(p,distr_param(df/2,1/(df/2)));
        theta.resize(k);
        theta(k-1) = randg(distr_param(ad2,1/bd2));
        tauh = cumprod(theta);
        Plam = psijh.each_row() % tauh.t();
      } else {
        if ((num > 0) && (num < k)) {
          facts = {k - num,1};
          k = max(facts);
          Lambda = Lambda.cols(find(lind < prop));
          psijh = psijh.cols(find(lind < prop));
          eta = eta.cols(find(lind < prop));
          theta = theta.elem(find(lind < prop));
          tauh = cumprod(theta);
          Plam = psijh.each_row() % tauh.t();
        }
      }
    }
    
    thincheck = i - std::floor(i/thin) * thin; // % operator stolen by arma
    if(!thincheck && (i >= burn)) {
      if(covm || cov) Omega = (Lambda * Lambda.t() + Sigma) % scaleMat;
      if(covm) COVMEAN += Omega / std::max(sp+1, 1);
      if(cov) OMEGA.slice(ind) = Omega;
      if(fac) LAMBDA(ind) = Lambda;
      if(sig) SIGMA.slice(ind) = Sigma;
      if(nfa) K(ind) = k;
      ind += 1;
    }
    
    printcheck = (start+1) % 1000;
    if(!printcheck && verbose){
      Rcpp::Rcout << (start+1) << "\n";
    }
  }
  
  Rcpp::List ls = Rcpp::List::create(ps, Sigma, Lambda, meta, 
                                     veta, psijh, theta, tauh, 
                                     Plam, k);
  
  return Rcpp::List::create(Rcpp::Named("covMean") = COVMEAN,
                            Rcpp::Named("covSamps") = OMEGA,
                            Rcpp::Named("factSamps") = LAMBDA,
                            Rcpp::Named("sigSamps") = SIGMA,
                            Rcpp::Named("numFact") = K,
                            Rcpp::Named("lastState") = ls);
}

