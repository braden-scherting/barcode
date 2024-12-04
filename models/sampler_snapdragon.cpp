# include <RcppArmadillo.h>
# include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <fstream>

arma::vec oneMultinomCalt(arma::vec probs, int size, int L){
  Rcpp::IntegerVector ans(L);
  rmultinom(size, probs.begin(), L, ans.begin());
  return(Rcpp::as<arma::vec>(ans));
}

void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

Rcpp::List JpresentList_zinnia(arma::umat Ypresent, int P){
  Rcpp::List YjHere;
  
  for (int j=0; j<P; j++){
    YjHere.insert(j, find(Ypresent.col(j)));
  }
  return YjHere;
}

Rcpp::List IpresentList_zinnia(arma::umat Ypresent, int N){
  Rcpp::List YiHere;
  
  for (int i=0; i<N; i++){
    YiHere.insert(i, find(Ypresent.row(i)));
  }
  return YiHere;
}
  
arma::cube Ypart_snapdragon(arma::mat Y, arma::mat Gamma, arma::mat S, arma::mat Phi, 
                                    arma::mat C, int N, int P, int L, arma::umat Ypresent){
  
  arma::cube Ypart(N,P,L, arma::fill::zeros);
  arma::vec El(L);
  arma::vec parts(L);
  
  for (int i=0; i<N; i++){
    for (int j=0; j<P; j++){
      if (Ypresent(i,j)){
        El = (C.row(i)%Phi.row(i)).t() % (S.row(j)%Gamma.row(j)).t();
        if (max(El)==0) continue;
        parts = oneMultinomCalt(El/accu(El), Y(i,j), L);
        Ypart(arma::span(i), arma::span(j), arma::span::all) = parts; 
      }
    }
  }
  
  return Ypart;
}

double logProbZeroPois(arma::vec y, arma::vec Ey){
  double lp=0;
  arma::uvec EyPlus = find(Ey>0);
  lp = sum(y(EyPlus) % arma::log(Ey(EyPlus)) - Ey(EyPlus));
  y.shed_rows(EyPlus);
  lp += std::log(all(y==0));
  return lp;
}


arma::mat Gamma_snapdragon(arma::cube Ypart, arma::mat S, arma::mat Phi, arma::mat C, int P, int L, 
                       double aGamma, arma::vec bGamma, Rcpp::List YjHere, double tune_rate=20.0){
  
  arma::mat Gamma(P, L);
  
  arma::mat ysum = sum(Ypart, 0);
  arma::vec CPhil;
  double cphisum;
    
  for (int l=0; l<L; l++){
    CPhil = Phi.col(l)%C.col(l);
    for (int j=0; j<P; j++){
      if (S(j,l)==0){
        Gamma(j,l) = R::rgamma(1, 1/tune_rate);
      } else {
        arma::uvec present = YjHere[j];
        cphisum = sum(CPhil.elem(present));
        
        Gamma(j,l) = R::rgamma(aGamma+ysum(j,l), 1/(bGamma(l)+cphisum));
          // R::qgamma(R::runif(
          // R::pgamma(0.1, aGamma+ysum(j,l), 1/(bGamma(l)+cphisum), 1, 0), 1),
          // aGamma+ysum(j,l), 1/(bGamma(l)+cphisum),1,0); 
      }
    }
  }
  return Gamma;
}

arma::mat Phi_snapdragon(arma::cube Ypart, arma::mat Gamma, arma::mat S, arma::mat C, int N, int L,
                         double aPhi, arma::vec bPhi, Rcpp::List YiHere, double tune_rate=20.0){
  arma::mat Phi(N, L);
  
  arma::mat ysum = sum(Ypart, 1);
  arma::vec SGam;
  double sgamsum;
  
  for (int l=0; l<L; l++){
    SGam = S.col(l)%Gamma.col(l);
    for (int i=0; i<N; i++){
      if (C(i,l)==0){
        Phi(i,l) = R::rgamma(1, 1/tune_rate);
      } else{
        arma::uvec present = YiHere[i];
        sgamsum = sum(SGam.elem(present));
        
        Phi(i,l) = R::rgamma(aPhi+ysum(i,l), 1/(bPhi(l)+sgamsum));
        // Phi(i,l) = R::qgamma(R::runif(
        //   R::pgamma(0.1, aPhi+ysum(i,l), 1/(bPhi(l)+sgamsum), 1, 0), 1),
        //   aPhi+ysum(i,l), 1/(bPhi(l)+sgamsum), 1, 0);
      }
    }
  }
  return Phi;
}

arma::vec bPhi_snapdragon(arma::mat C, arma::mat Phi, int L, double aPhi=0.5, double a0=0.5, double b0=0.5){
  arma::vec bPhi(L);
  arma::vec sumPhiC = sum(Phi%C, 0).t();
  arma::vec sumC = sum(C, 0).t();
  
  for (int l=0; l<L; l++){
    bPhi(l) = R::rgamma(a0 + sumC(l)*aPhi, 1/(b0 + sumPhiC(l)));
  }
  return bPhi;
}

arma::mat Zeta_snapdragon(arma::cube Ypart, arma::mat Gamma, arma::mat S, arma::mat C, arma::mat Zeta, int N, int L,
                          double aZeta, Rcpp::List YiHere){
  arma::mat newZeta(N,L);
  
  double ul;
  arma::mat ySimNL = sum(Ypart, 1);
  arma::vec ySumL = sum(ySimNL, 0).t();
  
  for (int l=0; l<L; l++){
    ul = R::rgamma(ySumL(l), 1 / sum(C.col(l) % Zeta.col(l)));
    for (int i=0; i<N; i++){
      if (C(i,l)==0){
        newZeta(i,l) = R::rgamma(aZeta, 1);
      } else {
        newZeta(i,l) = R::rgamma(aZeta + ySimNL(i,l), 1/(1 + ul));
      }
    }
  }
  return newZeta;
}

arma::mat PhiZeta_snapdragon(arma::mat C, arma::mat Zeta, int N, int L){
  arma::mat Phi(N,L);
  for (int l=0; l<L; l++){
    Phi.col(l) = Zeta.col(l) / sum(C.col(l) % Zeta.col(l));
  }
  return Phi;
}

arma::mat Cuni_snapdragon(arma::mat Y, arma::mat Gamma, arma::mat S, arma::mat C, arma::mat Zeta, 
                          arma::mat PrC, int N, int L, Rcpp::List YiHere){
  arma::vec mu0;
  arma::vec mu1;
  
  arma::vec Tl = sum(C%Zeta, 0).t();
  arma::mat Lam = Gamma % S;
  
  for (int i=0; i<N; i++){
    arma::uvec present = YiHere[i];
    arma::vec Yi = Y.row(i).t();
    Yi = Yi(present);
    
    for (int l=0; l<L; l++){
      arma::vec Ci = C.row(i).t();
      
      arma::vec CiZero = Ci; CiZero(l) = 0;
      arma::vec CiOne = Ci; CiOne(l) = 1;

      arma::vec TlZero = Tl; TlZero(l) -= (C(i,l)*Zeta(i,l));
      arma::vec TlOne = Tl; TlOne(l) += ((1-C(i,l))*Zeta(i,l));
      
      mu0 = Lam * ((CiZero % Zeta.row(i).t()) / TlZero);
      mu1 = Lam * ((CiOne % Zeta.row(i).t()) / TlOne);

      double lpc0 = std::log(1-PrC(i,l)) + logProbZeroPois(Yi, mu0(present));
      double lpc1 = std::log(PrC(i,l)) + logProbZeroPois(Yi, mu1(present));
      
      arma::vec u = Rcpp::runif(2, 0., 1.);
      arma::vec z = -arma::log(-arma::log(u));
      z(0) += lpc0;
      z(1) += lpc1;
      C(i,l) = z.index_max();
      Tl(l) = sum(C.col(l) % Zeta.col(l));
    }
  }
  return C;
}

arma::vec bGamma_snapdragon(arma::mat S, arma::mat Gamma, int L, double aGamma=0.5, double a0=0.5, double b0=0.5){
  
  arma::vec bGamma(L);
  arma::vec sumGamS = sum(Gamma%S, 0).t();
  arma::vec sumS = sum(S, 0).t();
  
  for (int l=0; l<L; l++){
    bGamma(l) = R::rgamma(a0 + sumS(l)*aGamma, 1 / (b0 + sumGamS(l)));
  }
  
  return bGamma;
}

double psiS_snapdragon(arma::mat S, int P, int L){
  double psiS = R::rbeta(10 + accu(S), 10 + (P*L) - accu(S));
  return psiS;
}

arma::mat Beta_snapdragon(arma::mat X, arma::mat des01, arma::mat W, arma::mat Xi, int L, int Q, double tauB=4){
  arma::mat B(Q, L);
  arma::vec diagB(Q, arma::fill::value(tauB));
  arma::mat prior_prec=diagmat(diagB);
  arma::mat data_prec=X.t()*X;
  
  arma::mat Sigma = arma::inv_sympd(prior_prec + data_prec);
  arma::mat cSigma  = arma::chol(Sigma, "lower");
  
  arma::vec W_adj_l(Q);
  arma::vec mu_l(Q);
  arma::vec std_norm(Q);
  
  for (int l=0; l<L; l++){
    W_adj_l = W.col(l) - (des01 * Xi.col(l));
    mu_l = Sigma * X.t() * W_adj_l;
    std_norm = Rcpp::rnorm(Q, 0., 1.);
    B.col(l) = cSigma * std_norm + mu_l;
  }
  
  return B;
}


arma::mat XiIndep_snapdragon(arma::mat X, arma::mat des01, arma::mat W, arma::mat Beta,
                        arma::cube condVar, arma::cube condVarChol, int N, int L, int Q, int K){
  arma::mat Xinew(K,L);
  arma::mat DtWadj = des01.t() * (W - (X*Beta));
  
  for (int l=0; l<L; l++){
    arma::vec Z=Rcpp::rnorm(K, 0, 1);
    Xinew.col(l)=(condVarChol.slice(l) * Z + (condVar.slice(l) * DtWadj.col(l)));
  }
  
  return Xinew;
}


Rcpp::NumericVector R_Rcpp_truncnorm(arma::vec a, arma::vec b, arma::vec m, int N){
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("truncnorm");
  Rcpp::Function f = pkg["rtruncnorm"];
  return f( N, Rcpp::Named("a", a), Rcpp::Named("b", b), Rcpp::Named("mean", m));
}


arma::mat W_snapdragon(arma::mat X, arma::mat des01, arma::mat C, arma::mat Beta, 
                      arma::mat Xi, int N, int L){
  
  arma::mat M = X*Beta + des01*Xi;
  arma::mat A(N, L, arma::fill::zeros);
  arma::mat B(N, L, arma::fill::zeros);
  
  A.elem( find(C==0) ).fill(R_NegInf);
  B.elem( find(C==1) ).fill(R_PosInf);
  
  arma::mat out(N,L);
  Rcpp::NumericVector outl(N);
  
  for (int l=0; l<L; l++){
    outl = R_Rcpp_truncnorm(A.col(l), B.col(l), M.col(l), N);
    out.col(l) = Rcpp::as<arma::vec>(outl);
  }
  
  return out;
}


arma::mat ProbitToEtaN_snapdragon(arma::mat X, arma::mat des01, arma::mat Beta, 
                             arma::mat Xi, int N, int L){
  arma::mat EtaN(N, L);
  arma::mat What = X*Beta + des01*Xi;
  for (int i=0; i<N; i++){
    for (int l=0; l<L; l++){
      EtaN(i,l) = R::pnorm(What(i,l), 0, 1, true, false);
    }
  }
  return EtaN;
}

// arma::mat CRandomBlocked_snapdragon(arma::mat Y, arma::mat Gamma, arma::mat S, arma::mat Phi, 
//                           arma::mat C, double aPhi, arma::vec bPhi, arma::mat PrC, int N, int P, int L, 
//                           Rcpp::List YiHere, arma::mat CstarAll, double tune_rate=20.0){
arma::mat CRandomBlocked_snapdragon(arma::mat Y, arma::mat Gamma, arma::mat S, arma::mat Phi, 
                                    arma::mat C, arma::mat Zeta, arma::mat PrC, int N, int P, int L, 
                                    Rcpp::List YiHere, arma::mat CstarAll){
  
  arma::mat Lambda = S % Gamma;
  
  int blockSize = CstarAll.n_cols;
  
  int B = ceil((L-1)/ (1.0*blockSize)); // number of blocks
  arma::vec block_nelem(B, arma::fill::value(blockSize));
  block_nelem(B-1) = L - 1 - (B-1)*blockSize;
  arma::vec block_nbinary = arma::exp2(block_nelem);
  
  for (int i=0; i<N; i++){
    
    arma::uvec present = YiHere[i];
    arma::vec Yi = Y.row(i).t();
    Yi = Yi(present);
    
    arma::rowvec Phi_i = Phi.row(i);
    arma::rowvec Zeta_i = Zeta.row(i);
    arma::uvec cells = arma::regspace<arma::uvec>(1,L-1);
    arma::uvec blk;
    
    for (int b=0; b<B; b++){
      arma::rowvec Ci = C.row(i);
      
      arma::vec logProbs(block_nbinary(b), arma::fill::zeros);
      // arma::mat CstarPhi(block_nbinary(b), L);
      arma::mat Mui(block_nbinary(b),P);
      
      // CstarPhi = arma::repelem(C.row(i)%Phi.row(i), block_nbinary(b), 1);
      
      arma::mat Cstar = CstarAll(arma::span(0, block_nbinary(b)-1), arma::span(blockSize-block_nelem(b), blockSize-1));
      
      arma::uvec tmp = arma::sort(Rcpp::RcppArmadillo::sample(arma::regspace<arma::uvec>(0,cells.n_elem-1), Cstar.n_cols, 0));
      blk = cells(tmp);
      
      arma::rowvec PrCi=PrC.row(i);
      
      // CstarPhi.cols(blk) = Cstar.each_row() % Phi_i(blk).t();
      // Mui = CstarPhi * Lambda.t();
      arma::mat OneMinusCstar = 1 - Cstar;
      
      logProbs += sum(arma::log((Cstar.each_row() % PrCi(blk).t()) +
                      (OneMinusCstar.each_row() % (1-PrCi(blk).t()))), 1);
                      //   +
                      // (Cstar.each_row() % (aPhi*arma::log(bPhi(blk)) + (aPhi-1)*arma::log(Phi_i(blk)) - bPhi(blk)%Phi_i(blk) - std::lgamma(aPhi)).t()) +
                      // (OneMinusCstar.each_row() % (std::log(tune_rate) - tune_rate*Phi_i(blk)).t()), 1);
       
      // 3) likelihood
      
      arma::rowvec tmpTs = sum(C % Zeta, 0);
      tmpTs(blk) -= Ci(blk) % Zeta_i(blk);
      
      for (int bit=0; bit<block_nbinary(b); bit++){
        // arma::vec tmpMu = Mui.row(bit).t();
        arma::rowvec tmpTsbit = tmpTs; 
        tmpTsbit(blk) += Cstar(arma::span(bit), arma::span::all)%Zeta_i(blk).t();
        arma::rowvec XcX = C.row(i) % Zeta_i; 
        XcX(blk) = Cstar(arma::span(bit), arma::span::all)%Zeta_i(blk).t();
        XcX = XcX / tmpTsbit;
        
        arma::vec tmpMu = Lambda * XcX.t();
        
        logProbs(bit) += logProbZeroPois(Yi, tmpMu.elem(present));
      }
      
      arma::vec u = Rcpp::runif(block_nbinary(b), 0., 1.);
      arma::vec z = -arma::log(-arma::log(u));
      z += logProbs;
      arma::rowvec winner = Cstar.row(z.index_max());
      for (int xx=0; xx<Cstar.n_cols; xx++){
        C(i, blk(xx)) = winner(xx);
      }
      cells.shed_rows(tmp);
    }
  }
  return C;
}

arma::mat SRandomBlocked_snapdragon(arma::mat Y, arma::mat Gamma, arma::mat S, arma::mat Phi, arma::mat C, 
                                double aGamma, arma::vec bGamma, double psiS, int N, int P, int L, 
                                Rcpp::List YjHere, arma::mat SstarAll, double tune_rate=20.0){
  
  arma::mat Factors = C % Phi;
  int blockSize = SstarAll.n_cols; // 2
  
  int B = ceil(L/ (1.0*blockSize)); // number of blocks
  arma::vec block_nelem(B, arma::fill::value(blockSize));
  block_nelem(B-1) = L - (B-1)*blockSize;
  arma::vec block_nbinary = arma::exp2(block_nelem);
  
  // block_nbinary(B-1) = pow(2.0, block_nelem(B-1));
  
  for (int j=0; j<P; j++){
    arma::uvec present = YjHere[j];
    arma::vec Yj = Y.col(j);
    Yj = Yj(present);
    
    arma::rowvec Gamma_j = Gamma.row(j);
    arma::uvec cells = arma::regspace<arma::uvec>(0,L-1);
    arma::uvec blk;

    for (int b=0; b<B; b++){
      arma::vec logProbs(block_nbinary(b), arma::fill::zeros);
      arma::mat SstarGamma(block_nbinary(b), L);
      arma::mat Mui(block_nbinary(b),P);
      
      SstarGamma = arma::repelem(S.row(j)%Gamma.row(j), block_nbinary(b), 1);
      
      arma::mat Sstar = SstarAll(arma::span(0, block_nbinary(b)-1), arma::span(blockSize-block_nelem(b), blockSize-1));
      
      arma::uvec tmp = arma::sort(Rcpp::RcppArmadillo::sample(arma::regspace<arma::uvec>(0,cells.n_elem-1), block_nelem(b), 0));
      blk = cells(tmp);
      
      SstarGamma.cols(blk) = Sstar.each_row() % Gamma_j(blk).t();
      
      Mui = SstarGamma * Factors.t();

      arma::mat OneMinusSstar = 1-Sstar;

      logProbs += sum(arma::log((Sstar * psiS) +
        (OneMinusSstar * (1-psiS))) +
        (Sstar.each_row() % (aGamma*arma::log(bGamma(blk)) + (aGamma-1)*arma::log(Gamma_j(blk)) - bGamma(blk)%Gamma_j(blk) - std::lgamma(aGamma)).t()) +
        (OneMinusSstar.each_row() % (std::log(tune_rate) - tune_rate*Gamma_j(blk)).t()), 1);
      
      for (int bit=0; bit<block_nbinary(b); bit++){
        arma::vec tmpMu = Mui.row(bit).t();
        logProbs(bit) += logProbZeroPois(Yj, tmpMu.elem(present));
      }
      
      arma::vec u = Rcpp::runif(block_nbinary(b), 0., 1.);
      arma::vec z = -arma::log(-arma::log(u));
      z += logProbs;
      
      arma::rowvec winner = Sstar.row(z.index_max()); 
      
      for (int xx=0; xx<blk.n_elem; xx++){
        S(j, blk(xx)) = winner(xx); 
      }
      cells.shed_rows(tmp);
    }
  }
  return S;
}


// [[Rcpp::export]]
Rcpp::List sampler_snapdragon(int niter, int nthin, int nburn, Rcpp::List dims, Rcpp::List dat, 
                        Rcpp::List storeInit, unsigned int workerSeed, int x=0, int pred=0) {
  
  set_seed(workerSeed);
  
  std::string file_name = "progressFiles/out" + std::to_string(x) + ".txt";
  
  arma::mat Y = Rcpp::as<arma::mat>(dat["Y"]);
  arma::mat X = Rcpp::as<arma::mat>(dat["X"]);
  arma::mat des01 = Rcpp::as<arma::mat>(dat["des01"]);

  arma::cube condVar=Rcpp::as<arma::cube>(storeInit["condVar"]);
  arma::cube condVarChol=Rcpp::as<arma::cube>(storeInit["condVarChol"]);
  
  arma::mat Cstar = Rcpp::as<arma::mat>(storeInit["Cstar"]);
  arma::mat Sstar = Rcpp::as<arma::mat>(storeInit["Cstar"]);
  
  double a0 = 0.5;
  double b0 = 0.5;
  // double aPhi = 0.5;
  double aGamma = 0.5;
  
  arma::umat Ypresent = Rcpp::as<arma::umat>(storeInit["Ypresent"]);
  
  int N = Rcpp::as<int>(dims["N"]);
  int P = Rcpp::as<int>(dims["P"]);
  int L = Rcpp::as<int>(dims["L"]);
  int Q = Rcpp::as<int>(dims["Q"]);
  int K = Rcpp::as<int>(dims["K"]);
  
  arma::cube C(N, L, niter);
  arma::mat Cit(N, L, arma::fill::ones);
  
  arma::cube S(P, L, niter);
  arma::mat Sit(P, L, arma::fill::ones);
  
  arma::vec psiS(niter);
  double psiSit=1.0/2.0;
  
  arma::cube Gamma(P, L, niter);
  arma::mat Gammait = Rcpp::as<arma::mat>(storeInit["GammaInit"]);
  
  arma::cube Phi(N, L, niter);
  arma::mat Phiit = Rcpp::as<arma::mat>(storeInit["PhiInit"]);
  
  arma::cube Zeta(N, L, niter);
  arma::mat Zetait = Rcpp::as<arma::mat>(storeInit["PhiInit"]);
  // for (int l=0; l<L; l++){
  //   Zetait.col(l) = Phiit.col(l) / sum(Phiit.col(l));
  // }
  
  arma::mat bGamma(L, niter);
  arma::vec bGammait(L, arma::fill::value(0.5));
  
  arma::mat bPhi(L, niter);
  arma::vec bPhiit(L, arma::fill::value(0.5));
   
  arma::cube Eta(N, L, niter);
  arma::mat EtaNit(N,L,arma::fill::value(0.5));
  
  arma::cube Beta(Q, L, niter);
  arma::mat Betait(Q, L, arma::fill::zeros);
   
  arma::cube Xi(K, L, niter);
  arma::mat Xiit(K, L, arma::fill::zeros);

  arma::cube W(N, L, niter);
  arma::mat Wit(N, L, arma::fill::zeros);
  
  arma::cube Ypart(N, P, L);
  
  Rcpp::List YiHere = IpresentList_zinnia(Ypresent, N);
  Rcpp::List YjHere = JpresentList_zinnia(Ypresent, P);
  
  
  
  for (int xx_burn=0; xx_burn<nburn; xx_burn++){
    Sit = SRandomBlocked_snapdragon(Y, Gammait, Sit, Phiit, Cit, aGamma, bGammait, psiSit, N, P, L, YjHere, Sstar);
    Cit = CRandomBlocked_snapdragon(Y, Gammait, Sit, Phiit, Cit, Zetait, EtaNit, N, P, L, YiHere, Cstar);
    // Cit = Cuni_snapdragon(Y, Gammait, Sit, Cit, Zetait, EtaNit, N, L, YiHere);
    Phiit = PhiZeta_snapdragon(Cit, Zetait, N, L);
    
    Ypart = Ypart_snapdragon(Y, Gammait, Sit, Phiit, Cit, N, P, L, Ypresent);
    Gammait = Gamma_snapdragon(Ypart, Sit, Phiit, Cit, P, L, aGamma, bGammait, YjHere);
    // Phiit = Phi_snapdragon(Ypart, Gammait, Sit, Cit, N, L, aPhi, bPhiit, YiHere);
    Zetait = Zeta_snapdragon(Ypart, Gammait, Sit, Cit, Zetait, N, L, 1, YiHere);
    Phiit = PhiZeta_snapdragon(Cit, Zetait, N, L);
    
    bGammait = bGamma_snapdragon(Sit, Gammait, L, aGamma, a0, b0);
    // bPhiit = bPhi_snapdragon(Cit, Phiit, L, aPhi, a0, b0);
    psiSit = psiS_snapdragon(Sit, P, L);
    
    Betait = Beta_snapdragon(X, des01, Wit, Xiit, L, Q);
    Xiit = XiIndep_snapdragon(X, des01, Wit, Betait, condVar, condVarChol, N, L, Q, K);
    Wit = W_snapdragon(X, des01, Cit, Betait, Xiit, N, L);
    EtaNit = ProbitToEtaN_snapdragon(X, des01, Betait, Xiit, N, L);

    if (xx_burn % 100 == 99){
      std::ofstream out(file_name.c_str());
      out << "Burn iteration: " << xx_burn+1 << " / " << nburn << "\n";
      out.close();
    }
  }
  
  for (int it=0; it<niter; it++){
    for (int xit=0; xit<nthin; xit++){
      Sit = SRandomBlocked_snapdragon(Y, Gammait, Sit, Phiit, Cit, aGamma, bGammait, psiSit, N, P, L, YjHere, Sstar);
      Cit = CRandomBlocked_snapdragon(Y, Gammait, Sit, Phiit, Cit, Zetait, EtaNit, N, P, L, YiHere, Cstar);
      // Cit = Cuni_snapdragon(Y, Gammait, Sit, Cit, Zetait, EtaNit, N, L, YiHere);
      Phiit = PhiZeta_snapdragon(Cit, Zetait, N, L);
      
      Ypart = Ypart_snapdragon(Y, Gammait, Sit, Phiit, Cit, N, P, L, Ypresent);
      Gammait = Gamma_snapdragon(Ypart, Sit, Phiit, Cit, P, L, aGamma, bGammait, YjHere);
      // Phiit = Phi_snapdragon(Ypart, Gammait, Sit, Cit, N, L, aPhi, bPhiit, YiHere);
      Zetait = Zeta_snapdragon(Ypart, Gammait, Sit, Cit, Zetait, N, L, 1, YiHere);
      Phiit = PhiZeta_snapdragon(Cit, Zetait, N, L);
      
      bGammait = bGamma_snapdragon(Sit, Gammait, L, aGamma, a0, b0);
      // bPhiit = bPhi_snapdragon(Cit, Phiit, L, aPhi, a0, b0);
      psiSit = psiS_snapdragon(Sit, P, L);
      
      Betait = Beta_snapdragon(X, des01, Wit, Xiit, L, Q);
      Xiit = XiIndep_snapdragon(X, des01, Wit, Betait, condVar, condVarChol, N, L, Q, K);
      Wit = W_snapdragon(X, des01, Cit, Betait, Xiit, N, L);
      EtaNit = ProbitToEtaN_snapdragon(X, des01, Betait, Xiit, N, L);
    }
    
    
    C.slice(it)=Cit;
    S.slice(it)=Sit;
    Gamma.slice(it)=Gammait;
    Phi.slice(it)=Phiit;
    Zeta.slice(it)=Zetait;
    Eta.slice(it)=EtaNit;
    bGamma.col(it)=bGammait;
    bPhi.col(it)=bPhiit;
    psiS(it)=psiSit;
    Beta.slice(it)=Betait;
    Xi.slice(it)=Xiit;
      
    if (it % 100 == 99){
      std::ofstream out(file_name.c_str());
      out << "Thinned iteration: " << it+1 << " / " << niter << "\n";
      out.close();
    }
  }
 
 return Rcpp::List::create(Rcpp::Named("C") = C,
                           Rcpp::Named("S") = S,
                           Rcpp::Named("Gamma") = Gamma,
                           Rcpp::Named("Phi") = Phi,
                           Rcpp::Named("Zeta") = Zeta,
                           Rcpp::Named("Eta") = Eta,
                           Rcpp::Named("bPhi") = bPhi,
                           Rcpp::Named("bGamma") = bGamma,
                           Rcpp::Named("psiS") = psiS,
                           Rcpp::Named("Beta") = Beta,
                           Rcpp::Named("Xi") = Xi);
}

