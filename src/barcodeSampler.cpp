# include <RcppArmadillo.h>
# include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <fstream>

// Utility functions
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

arma::vec oneMultinomCalt(arma::vec probs, int size, int L){
  Rcpp::IntegerVector ans(L);
  rmultinom(size, probs.begin(), L, ans.begin());
  return(Rcpp::as<arma::vec>(ans));
}

double logProbZeroPois(arma::vec y, arma::vec Ey){
  double lp=0;
  arma::uvec EyPlus = find(Ey>0);
  lp = sum(y(EyPlus) % arma::log(Ey(EyPlus)) - Ey(EyPlus));
  y.shed_rows(EyPlus);
  lp += std::log(all(y==0));
  return lp;
}

arma::sp_mat SparseToSparse(Rcpp::S4 RspY) {

  Rcpp::IntegerVector dims = RspY.slot("Dim");
  arma::urowvec i = Rcpp::as<arma::urowvec>(RspY.slot("i"));
  arma::urowvec j = Rcpp::as<arma::urowvec>(RspY.slot("p"));
  arma::vec y = Rcpp::as<arma::vec>(RspY.slot("x"));

  int n = dims[0];
  int p = dims[1];
  arma::sp_mat Ysp(i, j, y, n, p);

  return Ysp;
}

Rcpp::NumericVector R_Rcpp_truncnorm(arma::vec a, arma::vec b, arma::vec m, int N){
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("truncnorm");
  Rcpp::Function f = pkg["rtruncnorm"];
  return f( N, Rcpp::Named("a", a), Rcpp::Named("b", b), Rcpp::Named("mean", m));
}

// Sampler object
struct markov{
  // dimensions
  int N, P, L, Q, K, Peff;

  // data
  arma::mat Y, X, des01;
  arma::sp_mat Ysp;
  arma::cube condVar, condVarChol;

  // parameters
  arma::mat YpartN, YpartP;

  arma::mat C, S; // binary switches
  arma::mat Phi, Gamma, Zeta; // continuous scores
  arma::vec bGamma;
  double psiS;
  arma::mat Beta, Xi, W, Eta; // latent regression parameters

  // hyperparameters
  double a0, b0, aGamma;

  // functional aids
  arma::mat SstarAll;
};

// Sampling functions
void stepYparts(markov* pstate){
  pstate->YpartN.fill(0.0);
  pstate->YpartP.fill(0.0);

  arma::vec El(pstate->L), parts(pstate->L);
  arma::uword i, j;

  arma::sp_mat::const_iterator it     = pstate->Ysp.begin();
  arma::sp_mat::const_iterator it_end = pstate->Ysp.end();

  for(; it != it_end; ++it){
    i = it.row();
    j = it.col();
    El = (pstate->Phi.row(i)%pstate->C.row(i)).t() %
      (pstate->Gamma.row(j)%pstate->S.row(j)).t();
    parts = oneMultinomCalt(El/accu(El), (*it), pstate->L);
    pstate->YpartN.row(i) += parts.t();
    pstate->YpartP.row(j) += parts.t();
  }
}

void stepGamma(markov* pstate){
  for (int l=0; l<pstate->L; ++l){
    // arma::vec CPhil = pstate->Phi.col(l)%pstate->C.col(l);
    double cphisum = arma::accu(pstate->Phi.col(l)%pstate->C.col(l));
    for (int j=0; j<pstate->P; ++j){
      if (pstate->S(j,l)==0){
        pstate->Gamma(j,l) = R::rgamma(1., 1./20.);
      } else{
        pstate->Gamma(j,l) = R::rgamma(pstate->aGamma + pstate->YpartP(j,l),
                      1./(pstate->bGamma(l) + cphisum));
      }
    }
  }
}

void stepZetaPhi(markov* pstate){
  arma::vec ySumL = sum(pstate->YpartN, 0).t();

  for (int l=1; l<pstate->L; ++l){
    double ul = R::rgamma(ySumL(l), 1./ sum(pstate->C.col(l) % pstate->Zeta.col(l)));
    for (int i=0; i<pstate->N; ++i){
      if (pstate->C(i,l)==0){
        pstate->Zeta(i,l) = R::rgamma(1., 1.);
      } else{
        pstate->Zeta(i,l) = R::rgamma(1.+pstate->YpartN(i,l), 1./(1.+ul));
      }
    }
    pstate->Phi.col(l) = pstate->Zeta.col(l) / sum(pstate->C.col(l)%pstate->Zeta.col(l));
  }
}

void stepS(markov* pstate){

  arma::mat Factors = pstate->C % pstate->Phi;
  int blockSize = pstate->SstarAll.n_cols; // 2

  int B = ceil(pstate->L/ (1.0*blockSize)); // number of blocks
  arma::vec block_nelem(B, arma::fill::value(blockSize));
  block_nelem(B-1) = pstate->L - (B-1)*blockSize;
  arma::vec block_nbinary = arma::exp2(block_nelem);

  for (int j=0; j<pstate->P; ++j){
    arma::vec Yj = pstate->Y.col(j);

    arma::rowvec Gamma_j = pstate->Gamma.row(j);
    arma::uvec cells = arma::regspace<arma::uvec>(0,pstate->L-1);
    arma::uvec blk;

    for (int b=0; b<B; ++b){
      arma::vec logProbs(block_nbinary(b), arma::fill::zeros);
      arma::mat SstarGamma(block_nbinary(b), pstate->L);
      arma::mat Mui(block_nbinary(b), pstate->P);

      SstarGamma = arma::repelem(pstate->S.row(j)%pstate->Gamma.row(j), block_nbinary(b), 1);

      arma::mat Sstar = pstate->SstarAll(arma::span(0, block_nbinary(b)-1), arma::span(blockSize-block_nelem(b), blockSize-1));

      arma::uvec tmp = arma::sort(Rcpp::RcppArmadillo::sample(arma::regspace<arma::uvec>(0,cells.n_elem-1), block_nelem(b), 0));
      blk = cells(tmp);

      SstarGamma.cols(blk) = Sstar.each_row() % Gamma_j(blk).t();

      Mui = SstarGamma * Factors.t();

      arma::mat OneMinusSstar = 1-Sstar;

      logProbs += sum(arma::log((Sstar * pstate->psiS) + (OneMinusSstar * (1-pstate->psiS))) +
        (Sstar.each_row() % (pstate->aGamma*arma::log(pstate->bGamma(blk)) +
        (pstate->aGamma-1)*arma::log(Gamma_j(blk)) - pstate->bGamma(blk)%Gamma_j(blk) -
        std::lgamma(pstate->aGamma)).t()) + (OneMinusSstar.each_row() % (std::log(20.) - 20.*Gamma_j(blk)).t()), 1);

      for (int bit=0; bit<block_nbinary(b); ++bit){
        arma::vec tmpMu = Mui.row(bit).t();
        logProbs(bit) += logProbZeroPois(Yj, tmpMu);
      }

      arma::vec u = Rcpp::runif(block_nbinary(b), 0., 1.);
      arma::vec z = -arma::log(-arma::log(u));
      z += logProbs;

      arma::rowvec winner = Sstar.row(z.index_max());

      for (int xx=0; xx<blk.n_elem; xx++){
        pstate->S(j, blk(xx)) = winner(xx);
      }
      cells.shed_rows(tmp);
    }
  }
  return;
}

void stepC(markov* pstate){
  arma::vec Tl = sum(pstate->C%pstate->Zeta, 0).t();
  arma::mat Lam = pstate->Gamma%pstate->S;

  for (int i=0; i<pstate->N; ++i){
    arma::vec Yi = pstate->Y.row(i).t(); 
    
    for (int l=1; l<pstate->L; ++l){
      arma::vec CiZero = pstate->C.row(i).t(); CiZero(l)=0;
      arma::vec CiOne = CiZero; CiOne(l)=1.;

      arma::vec TlZero = Tl; TlZero(l) -= (pstate->C(i,l)*pstate->Zeta(i,l));
      arma::vec TlOne = Tl; TlOne(l) += ((1-pstate->C(i,l))*pstate->Zeta(i,l));

      arma::vec mu0 = Lam * ((CiZero % pstate->Zeta.row(i).t()) / TlZero);
      arma::vec mu1 = Lam * ((CiOne % pstate->Zeta.row(i).t()) / TlOne);

      double lpc0 = std::log(1-pstate->Eta(i,l)) + logProbZeroPois(Yi, mu0);
      double lpc1 = std::log(pstate->Eta(i,l)) + logProbZeroPois(Yi, mu1);

      arma::vec u = Rcpp::runif(2, 0., 1.);
      arma::vec z = -arma::log(-arma::log(u));
      z(0) += lpc0;
      z(1) += lpc1;
      pstate->C(i,l) = z.index_max();
      Tl(l) = sum(pstate->C.col(l) % pstate->Zeta.col(l));
    }
  }
  for (int l=1; l<pstate->L; ++l) pstate->Phi.col(l) = pstate->Zeta.col(l) / sum(pstate->C.col(l)%pstate->Zeta.col(l));
}

void stepRegression(markov* pstate){
  // Data augmentation
  arma::mat M = pstate->X*pstate->Beta + pstate->des01*pstate->Xi;
  arma::mat A(pstate->N, pstate->L, arma::fill::zeros);
  arma::mat B(pstate->N, pstate->L, arma::fill::zeros);

  A.elem( find(pstate->C==0) ).fill(R_NegInf);
  B.elem( find(pstate->C==1) ).fill(R_PosInf);
  Rcpp::NumericVector outl(pstate->N);

  for (int l=0; l<pstate->L; l++){
    outl = R_Rcpp_truncnorm(A.col(l), B.col(l), M.col(l), pstate->N);
    pstate->W.col(l) = Rcpp::as<arma::vec>(outl);
  }

  // Sample Beta
  arma::vec diagB(pstate->Q, arma::fill::value(5.));
  arma::mat prior_prec=diagmat(diagB);
  arma::mat data_prec=pstate->X.t()*pstate->X;

  arma::mat Sigma = arma::inv_sympd(prior_prec + data_prec);
  arma::mat cSigma  = arma::chol(Sigma, "lower");

  arma::vec W_adj_l(pstate->Q);
  arma::vec mu_l(pstate->Q);
  arma::vec std_norm(pstate->Q);

  for (int l=1; l<pstate->L; ++l){
    W_adj_l = pstate->W.col(l) - (pstate->des01 * pstate->Xi.col(l));
    mu_l = Sigma * pstate->X.t() * W_adj_l;
    std_norm = Rcpp::rnorm(pstate->Q, 0., 1.);
    pstate->Beta.col(l) = cSigma * std_norm + mu_l;
  }

  // Sample Xi
  arma::mat DtWadj = pstate->des01.t() * (pstate->W - (pstate->X*pstate->Beta));

  for (int l=1; l<pstate->L; ++l){
    arma::vec Z=Rcpp::rnorm(pstate->K, 0, 1);
    pstate->Xi.col(l)=(pstate->condVarChol.slice(l) * Z + (pstate->condVar.slice(l) * DtWadj.col(l)));
  }

  // Update Pr(c_il = 1)
  arma::mat What = pstate->X*pstate->Beta + pstate->des01*pstate->Xi;
  for (int i=0; i<pstate->N; ++i){
    for (int l=0; l<pstate->L; ++l){
      pstate->Eta(i,l) = R::pnorm(What(i,l), 0, 1, true, false);
    }
  }
}

void stepRegressionNonspatial(markov* pstate){
  // Data augmentation
  arma::mat M = pstate->X*pstate->Beta;
  arma::mat A(pstate->N, pstate->L, arma::fill::zeros);
  arma::mat B(pstate->N, pstate->L, arma::fill::zeros);
  
  A.elem( find(pstate->C==0) ).fill(R_NegInf);
  B.elem( find(pstate->C==1) ).fill(R_PosInf);
  Rcpp::NumericVector outl(pstate->N);
  for (int l=0; l<pstate->L; l++){
    outl = R_Rcpp_truncnorm(A.col(l), B.col(l), M.col(l), pstate->N);
    pstate->W.col(l) = Rcpp::as<arma::vec>(outl);
  }
  
  // Sample Beta
  arma::vec diagB(pstate->Q, arma::fill::value(5.));
  arma::mat prior_prec=diagmat(diagB);
  arma::mat data_prec=pstate->X.t()*pstate->X;
  
  arma::mat Sigma = arma::inv_sympd(prior_prec + data_prec);
  arma::mat cSigma  = arma::chol(Sigma, "lower");
  
  arma::vec W_adj_l(pstate->Q);
  arma::vec mu_l(pstate->Q);
  arma::vec std_norm(pstate->Q);
  for (int l=1; l<pstate->L; ++l){
    mu_l = Sigma * pstate->X.t() * pstate->W.col(l);
    std_norm = Rcpp::rnorm(pstate->Q, 0., 1.);
    pstate->Beta.col(l) = cSigma * std_norm + mu_l;
  }
  
  // Update Pr(c_il = 1)
  arma::mat What = pstate->X*pstate->Beta;
  for (int i=0; i<pstate->N; ++i){
    for (int l=0; l<pstate->L; ++l){
      pstate->Eta(i,l) = R::pnorm(What(i,l), 0, 1, true, false);
    }
  }
}

void stepbGamma(markov* pstate){
  arma::vec sumGamS = sum(pstate->Gamma%pstate->S, 0).t();
  arma::vec sumS = sum(pstate->S, 0).t();

  for (int l=0; l<pstate->L; ++l){
    pstate->bGamma(l) = R::rgamma(pstate->a0 + sumS(l)*pstate->aGamma,
                   1./(pstate->b0 + sumGamS(l)));
  }
}

void steppsiS(markov* pstate){
  pstate->psiS = R::rbeta(10 + accu(pstate->S), 10 + (pstate->P*pstate->L) - accu(pstate->S));
}

void stepAll(markov* pstate){
  stepYparts(pstate);
  stepGamma(pstate);
  stepZetaPhi(pstate);
  stepC(pstate);
  stepS(pstate);
  stepbGamma(pstate);
  steppsiS(pstate);
  stepRegression(pstate);
}

void stepAllOrdination(markov* pstate){
  stepYparts(pstate);
  stepGamma(pstate);
  stepZetaPhi(pstate);
  stepC(pstate);
  stepS(pstate);
  stepbGamma(pstate);
  steppsiS(pstate);
  stepRegressionNonspatial(pstate);
}

// [[Rcpp::export]]
Rcpp::List sampler_camellia (Rcpp::List dat, Rcpp::List dims, Rcpp::List storeInit, unsigned int workerSeed,
                       int niter=10, int nthin=1, int nburn=10, int pred=0, int update=10, int x=0){
  // Bookkeeping
  std::string file_name = "../output/progressFiles/chain" + std::to_string(x) + ".txt";
  set_seed(workerSeed);

  // Instantiate markov chain
  markov chain;
  markov* pchain;
  pchain = &chain;

  // dimensions
  int N = Rcpp::as<int>(dims["N"]); chain.N = N;
  int P = Rcpp::as<int>(dims["P"]); chain.P = P;
  int L = Rcpp::as<int>(dims["L"]); chain.L = L;
  int Q = Rcpp::as<int>(dims["Q"]); chain.Q = Q;
  int K = Rcpp::as<int>(dims["K"]); chain.K = K;
  
  // data
  Rcpp::S4 RspY = Rcpp::as<Rcpp::S4>(dat["RspY"]);
  chain.Ysp = SparseToSparse(RspY);
  arma::mat Ytmp(chain.Ysp); chain.Y = Ytmp;
  chain.X = Rcpp::as<arma::mat>(dat["X"]);
  chain.des01 = Rcpp::as<arma::mat>(dat["des01"]);

  chain.condVar=Rcpp::as<arma::cube>(storeInit["condVar"]);
  chain.condVarChol=Rcpp::as<arma::cube>(storeInit["condVarChol"]);

  chain.SstarAll=Rcpp::as<arma::mat>(storeInit["Cstar"]);

  // parameters
  chain.YpartN = arma::mat(N, L, arma::fill::zeros);
  chain.YpartP = arma::mat(P, L, arma::fill::zeros);

  chain.Gamma = Rcpp::as<arma::mat>(storeInit["Ginit"]);
  chain.S = arma::mat(P, L, arma::fill::ones);
  chain.bGamma = arma::vec(L, arma::fill::value(1.));
  chain.psiS = 1./2.;

  chain.Phi = Rcpp::as<arma::mat>(storeInit["Pinit"]);
  chain.Phi.col(0).fill(1./N);

  chain.Zeta = arma::mat(N, L, arma::fill::randu);
  chain.C = arma::mat(N, L, arma::fill::ones);
  chain.Beta = arma::mat(Q, L, arma::fill::zeros);
  chain.Beta.col(0).fill(10.0);

  chain.Xi = arma::mat(K, L, arma::fill::zeros);
  chain.Xi.col(0).fill(10.0);

  chain.W = arma::mat(N, L, arma::fill::zeros);
  chain.Eta = arma::mat(N, L, arma::fill::value(1./2.));

  // hyperparameters
  chain.a0 = 0.5;
  chain.b0 = 0.5;
  chain.aGamma = 0.5;

  // output
  arma::cube outGamma(P, L, niter);
  arma::cube outS(P, L, niter);
  arma::cube outPhi(N, L, niter);
  arma::cube outC(N, L, niter);

  arma::mat outbGamma(L, niter);
  arma::vec outpsiS(niter);

  arma::cube outBeta(Q, L, niter);
  arma::cube outXi(K, L, niter);

  for (int xx=0; xx<nburn; ++xx){
    stepAll(pchain);

    if (xx % 100 == 99){
      std::ofstream out(file_name.c_str());
      out << "Warmup iteration: " << xx+1 << " / " << nburn << "\n";
      out.close();
    }
  }

  for (int it=0; it<niter; ++it){
    for (int xit=0; xit<nthin; ++xit){
      stepAll(pchain);
    }

    outGamma.slice(it) = chain.Gamma;
    outPhi.slice(it) = chain.Phi;
    outC.slice(it) = chain.C;
    outS.slice(it) = chain.S;
    outBeta.slice(it) = chain.Beta;
    outXi.slice(it) = chain.Xi;
    outpsiS(it) = chain.psiS;
    outbGamma.col(it) = chain.bGamma;

    if (it % 100 == 99){
      std::ofstream out(file_name.c_str());
      out << "Thinned iteration: " << it+1 << " / " << niter << "\n";
      out.close();
    }
  }

  return Rcpp::List::create(Rcpp::Named("YpartN") = chain.YpartN,
                            Rcpp::Named("YpartP") = chain.YpartP,
                            Rcpp::Named("Gamma") = outGamma,
                            Rcpp::Named("Phi") = outPhi,
                            Rcpp::Named("bGamma") = outbGamma,
                            Rcpp::Named("C") = outC,
                            Rcpp::Named("S") = outS,
                            Rcpp::Named("psi") = outpsiS,
                            Rcpp::Named("Beta") = outBeta,
                            Rcpp::Named("Xi") = outXi);
}

// [[Rcpp::export]]
Rcpp::List sampler_camellia_ordination (Rcpp::List dat, Rcpp::List dims, Rcpp::List storeInit, unsigned int workerSeed,
                             int niter=10, int nthin=1, int nburn=10, int pred=0, int update=10, int x=0){
  // Bookkeeping
  std::string file_name = "../output/progressFiles/chain" + std::to_string(x) + ".txt";
  set_seed(workerSeed);
  
  // Instantiate markov chain
  markov chain;
  markov* pchain;
  pchain = &chain;
  
  // dimensions
  int N = Rcpp::as<int>(dims["N"]); chain.N = N;
  int P = Rcpp::as<int>(dims["P"]); chain.P = P;
  int L = Rcpp::as<int>(dims["L"]); chain.L = L;
  int Q = Rcpp::as<int>(dims["Q"]); chain.Q = Q;
  
  // data
  Rcpp::S4 RspY = Rcpp::as<Rcpp::S4>(dat["RspY"]);
  chain.Ysp = SparseToSparse(RspY);
  arma::mat Ytmp(chain.Ysp); chain.Y = Ytmp;
  chain.X = Rcpp::as<arma::mat>(dat["X"]);
  
  chain.SstarAll=Rcpp::as<arma::mat>(storeInit["Cstar"]);
  
  // parameters
  chain.YpartN = arma::mat(N, L, arma::fill::zeros);
  chain.YpartP = arma::mat(P, L, arma::fill::zeros);
  
  chain.Gamma = Rcpp::as<arma::mat>(storeInit["Ginit"]);
  chain.S = arma::mat(P, L, arma::fill::ones);
  chain.bGamma = arma::vec(L, arma::fill::value(1.));
  chain.psiS = 1./2.;
  
  chain.Phi = Rcpp::as<arma::mat>(storeInit["Pinit"]);
  chain.Phi.col(0).fill(1./N);
  
  chain.Zeta = arma::mat(N, L, arma::fill::randu);
  chain.C = arma::mat(N, L, arma::fill::ones);
  chain.Beta = arma::mat(Q, L, arma::fill::zeros);
  chain.Beta.col(0).fill(10.0);
  
  chain.W = arma::mat(N, L, arma::fill::zeros);
  chain.Eta = arma::mat(N, L, arma::fill::value(1./2.));
  
  // hyperparameters
  chain.a0 = 0.5;
  chain.b0 = 0.5;
  chain.aGamma = 0.5;
  
  // output
  arma::cube outGamma(P, L, niter);
  arma::cube outS(P, L, niter);
  arma::cube outPhi(N, L, niter);
  arma::cube outC(N, L, niter);
  
  arma::mat outbGamma(L, niter);
  arma::vec outpsiS(niter);
  
  arma::cube outBeta(Q, L, niter);
  
  for (int xx=0; xx<nburn; ++xx){
    stepAllOrdination(pchain);
    
    if (xx % 100 == 99){
      std::ofstream out(file_name.c_str());
      out << "Warmup iteration: " << xx+1 << " / " << nburn << "\n";
      out.close();
    }
  }
  
  for (int it=0; it<niter; ++it){
    for (int xit=0; xit<nthin; ++xit){
      stepAllOrdination(pchain);
    }
    
    outGamma.slice(it) = chain.Gamma;
    outPhi.slice(it) = chain.Phi;
    outC.slice(it) = chain.C;
    outS.slice(it) = chain.S;
    outBeta.slice(it) = chain.Beta;
    // outXi.slice(it) = chain.Xi;
    outpsiS(it) = chain.psiS;
    outbGamma.col(it) = chain.bGamma;
    
    if (it % 100 == 99){
      std::ofstream out(file_name.c_str());
      out << "Thinned iteration: " << it+1 << " / " << niter << "\n";
      out.close();
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("YpartN") = chain.YpartN,
                            Rcpp::Named("YpartP") = chain.YpartP,
                            Rcpp::Named("Gamma") = outGamma,
                            Rcpp::Named("Phi") = outPhi,
                            Rcpp::Named("bGamma") = outbGamma,
                            Rcpp::Named("C") = outC,
                            Rcpp::Named("S") = outS,
                            Rcpp::Named("psi") = outpsiS,
                            Rcpp::Named("Beta") = outBeta);
}

void stepImpute(markov* pstate, arma::umat Ymissing){
  for (int xx=0; xx<Ymissing.n_rows; ++xx){
    int i = Ymissing(xx,0); 
    int j = Ymissing(xx,1);
    double Eij = arma::accu((pstate->C.row(i)%pstate->Phi.row(i)) % (pstate->S.row(j)%pstate->Gamma.row(j)));
    int yij = R::rpois( Eij );
    pstate->Y(i,j) = yij;
    pstate->Ysp(i,j) = yij;
  }
}

// [[Rcpp::export]]
Rcpp::List sampler_camellia_missing (Rcpp::List dat, Rcpp::List dims, Rcpp::List storeInit, unsigned int workerSeed,
                             int niter=10, int nthin=1, int nburn=10, int pred=0, int update=10, int x=0){
  // Bookkeeping
  std::string file_name = "../output/progressFiles/chain" + std::to_string(x) + ".txt";
  set_seed(workerSeed);
  
  // Instantiate markov chain
  markov chain;
  markov* pchain;
  pchain = &chain;
  
  // dimensions
  int N = Rcpp::as<int>(dims["N"]); chain.N = N;
  int P = Rcpp::as<int>(dims["P"]); chain.P = P;
  int L = Rcpp::as<int>(dims["L"]); chain.L = L;
  int Q = Rcpp::as<int>(dims["Q"]); chain.Q = Q;
  int K = Rcpp::as<int>(dims["K"]); chain.K = K;
  
  // data
  Rcpp::S4 RspY = Rcpp::as<Rcpp::S4>(dat["RspY"]);
  chain.Ysp = SparseToSparse(RspY);
  arma::mat Ytmp(chain.Ysp); chain.Y = Ytmp;
  chain.X = Rcpp::as<arma::mat>(dat["X"]);
  chain.des01 = Rcpp::as<arma::mat>(dat["des01"]);
  
  arma::umat Ymissing = Rcpp::as<arma::umat>(storeInit["Ymissing"]) - 1;

  chain.condVar=Rcpp::as<arma::cube>(storeInit["condVar"]);
  chain.condVarChol=Rcpp::as<arma::cube>(storeInit["condVarChol"]);
  
  chain.SstarAll=Rcpp::as<arma::mat>(storeInit["Cstar"]);
  
  // parameters
  chain.YpartN = arma::mat(N, L, arma::fill::zeros);
  chain.YpartP = arma::mat(P, L, arma::fill::zeros);
  
  chain.Gamma = Rcpp::as<arma::mat>(storeInit["Ginit"]);
  chain.S = arma::mat(P, L, arma::fill::ones);
  chain.bGamma = arma::vec(L, arma::fill::value(1.));
  chain.psiS = 1./2.;
  
  chain.Phi = Rcpp::as<arma::mat>(storeInit["Pinit"]);
  chain.Phi.col(0).fill(1./N);
  
  chain.Zeta = arma::mat(N, L, arma::fill::randu);
  chain.C = arma::mat(N, L, arma::fill::ones);
  chain.Beta = arma::mat(Q, L, arma::fill::zeros);
  chain.Beta.col(0).fill(10.0);
  
  chain.Xi = arma::mat(K, L, arma::fill::zeros);
  chain.Xi.col(0).fill(10.0);
  
  chain.W = arma::mat(N, L, arma::fill::zeros);
  chain.Eta = arma::mat(N, L, arma::fill::value(1./2.));
  
  // hyperparameters
  chain.a0 = 0.5;
  chain.b0 = 0.5;
  chain.aGamma = 0.5;
  
  // output
  arma::cube outGamma(P, L, niter);
  arma::cube outS(P, L, niter);
  arma::cube outPhi(N, L, niter);
  arma::cube outC(N, L, niter);
  
  arma::mat outbGamma(L, niter);
  arma::vec outpsiS(niter);
  
  arma::cube outBeta(Q, L, niter);
  arma::cube outXi(K, L, niter);
  
  arma::mat outE(N, P, arma::fill::zeros);
  
  for (int xx=0; xx<nburn; ++xx){
    stepImpute(pchain, Ymissing);
    stepAll(pchain);
    
    if (xx % 100 == 99){
      std::ofstream out(file_name.c_str());
      out << "Warmup iteration: " << xx+1 << " / " << nburn << "\n";
      out.close();
    }
  }
  
  for (int it=0; it<niter; ++it){
    for (int xit=0; xit<nthin; ++xit){
      stepImpute(pchain, Ymissing);
      stepAll(pchain);
    }
    
    outGamma.slice(it) = chain.Gamma;
    outPhi.slice(it) = chain.Phi;
    outC.slice(it) = chain.C;
    outS.slice(it) = chain.S;
    outBeta.slice(it) = chain.Beta;
    outXi.slice(it) = chain.Xi;
    outpsiS(it) = chain.psiS;
    outbGamma.col(it) = chain.bGamma;
    
    outE += chain.Y;
    
    if (it % 100 == 99){
      std::ofstream out(file_name.c_str());
      out << "Thinned iteration: " << it+1 << " / " << niter << "\n";
      out.close();
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("YpartN") = chain.YpartN,
                            Rcpp::Named("YpartP") = chain.YpartP,
                            Rcpp::Named("Gamma") = outGamma,
                            Rcpp::Named("Phi") = outPhi,
                            Rcpp::Named("bGamma") = outbGamma,
                            Rcpp::Named("C") = outC,
                            Rcpp::Named("S") = outS,
                            Rcpp::Named("psi") = outpsiS,
                            Rcpp::Named("Beta") = outBeta,
                            Rcpp::Named("Xi") = outXi,
                            Rcpp::Named("E") = outE/niter);
}