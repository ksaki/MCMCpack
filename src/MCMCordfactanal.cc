// MCMCordfactanal.cc is C++ code to estimate an ordinal data 
// factor analysis model

#include <iostream>
#include "Scythe_Matrix.h"
#include "Scythe_Simulate.h"
#include "Scythe_Stat.h"
#include "Scythe_LA.h"
#include "Scythe_IDE.h"
#include "Scythe_Math.h"

using namespace SCYTHE;
using namespace std;


extern "C"{

// function called by R to fit model
void
ordfactanalpost (double* sam, const int* samrow, const int* samcol,
		 const int* XX, const int* Xrow, const int* Xcol,
		 const int* burnin, const int* mcmc,  const int* thin,
		 const double* tune, const int* seed, const int* verbose, 
		 const double* Lamstart, const int* Lamstartrow, 
		 const int* Lamstartcol, 
		 const double* gam, const int* gamrow, const int* gamcol,
		 const int* ncat, const int* ncatrow, const int* ncatcol,
		 const double* Lameq, const int* Lameqrow, 
		 const int* Lameqcol,
		 const double* Lamineq, const int* Lamineqrow, 
		 const int* Lamineqcol,
		 const double* Lampmean, const int* Lampmeanrow, 
		 const int* Lampmeancol,
		 const double* Lampprec, const int* Lampprecrow,
		 const int* Lamppreccol, const int* storelambda,
		 const int* storescores,
		 int* accepts, int* outswitch
		 ) {

  // put together matrices
  Matrix<int> X(Xcol[0], Xrow[0], XX);
  X = t(X);
  Matrix<double> Lambda(Lamstartcol[0], Lamstartrow[0], Lamstart);
  Lambda = t(Lambda);
  Matrix<double> gamma(gamcol[0], gamrow[0], gam);
  gamma = t(gamma);
  Matrix<int> ncateg(ncatcol[0], ncatrow[0], ncat);
  ncateg = t(ncateg);
  Matrix<double> Lambda_eq(Lameqcol[0], Lameqrow[0], Lameq);
  Lambda_eq = t(Lambda_eq);
  Matrix<double> Lambda_ineq(Lamineqcol[0], Lamineqrow[0], Lamineq);
  Lambda_ineq = t(Lambda_ineq);
  Matrix<double> Lambda_prior_mean(Lampmeancol[0], Lampmeanrow[0], Lampmean);
  Lambda_prior_mean = t(Lambda_prior_mean);
  Matrix<double> Lambda_prior_prec(Lamppreccol[0], Lampprecrow[0], Lampprec);
  Lambda_prior_prec = t(Lambda_prior_prec);
  

  // initialize seed (mersenne twister / use default seed unless specified)
  if(seed==0) set_mersenne_seed(5489UL);
  else set_mersenne_seed(seed[0]);
      
  // parameters
  int K = X.cols();  // number of manifest variables
  int N = X.rows();  // number of observations
  int D = Lambda.cols();  // number of factors (including constant)
  int tot_iter = burnin[0] + mcmc[0];  

  // constants 
  Matrix<double> I = eye<double>(D-1);
  Matrix<double> Lambda_free_indic = Matrix<double>(K, D);
  for (int i=0; i<(K*D); ++i){
    if (Lambda_eq[i] == -999) Lambda_free_indic[i] = 1.0;
  }

  // starting values for phi, Xstar, and gamma_p
  Matrix<double> phi = Matrix<double>(N,D-1);
  phi = cbind(ones<double>(N,1), phi);
  Matrix<double> Xstar = Matrix<double>(N,K);
  Matrix<double> gamma_p = gamma(_,0);

  // storage matrices (row major order)
  Matrix<double> Lambda_store;
  if (storelambda[0]==1){
    Lambda_store = Matrix<double>(mcmc[0]/thin[0],K*D);
  }
  Matrix<double> gamma_store = Matrix<double>(mcmc[0]/thin[0], 
					      gamrow[0]*gamcol[0]);
  Matrix<double> phi_store;
  if (storescores[0]==1){
    phi_store = Matrix<double>(mcmc[0]/thin[0], N*D);
  }

  int count = 0;
 
  ///////////////////
  // Gibbs Sampler //
  ///////////////////
  
  for (int iter=0; iter < tot_iter; ++iter){


    
    // sample Xstar
    for (int j=0; j<K; ++j){
        Matrix<double> X_mean = phi * t(Lambda(j,_));
      for (int i=0; i<N; ++i){
	if (X(i,j) == -999){ // if missing
	  Xstar(i,j) = rnorm(X_mean[i], 1.0);
	}
	else { // if not missing
	  Xstar(i,j) = rtnorm(X_mean[i], 1.0, 
			      gamma(X(i,j)-1, j), gamma(X(i,j), j));
	}
      }
    }



    // sample phi
    Matrix<double> Lambda_const = Lambda(_,0);
    Matrix<double> Lambda_rest = Lambda(0, 1, K-1, D-1);
    Matrix<double> phi_post_var = invpd(I + t(Lambda_rest) * Lambda_rest);
    Matrix<double> phi_post_C = cholesky(phi_post_var);
    for (int i=0; i<N; ++i){
      Matrix<double> phi_post_mean = phi_post_var * 
	(t(Lambda_rest)  * (t(Xstar(i,_))-Lambda_const));
      Matrix<double> phi_samp = gaxpy(phi_post_C, rnorm(D-1, 1), 
      			      phi_post_mean);
      for (int j=0; j<(D-1); ++j)
	phi(i,j+1) = phi_samp[j];
    }
    

    
    // sample gamma
    for (int j=0; j<K; ++j){ // do the sampling for each manifest var
      Matrix<double> X_mean = phi * t(Lambda(j,_));
      for (int i=2; i<(ncateg[j]); ++i){
	if (i==(ncateg[j]-1)){
	  gamma_p[i] = rtbnorm_combo(gamma(i,j), ::pow(tune[j], 2.0), 
				     gamma_p[i-1]);
	}
	else {
	  gamma_p[i] = rtnorm(gamma(i,j), ::pow(tune[j], 2.0), gamma_p[i-1], 
			      gamma(i+1, j));
	}
      }
      double loglikerat = 0.0;
      double loggendenrat = 0.0;

          
      // loop over observations and construct the acceptance ratio
      for (int i=0; i<N; ++i){
	if (X(i,j) != -999){
	  if (X(i,j) == ncateg[j]){
	    loglikerat = loglikerat 
	      + log(1.0  - 
		    pnorm1(gamma_p[X(i,j)-1] - X_mean[i]) ) 
	      - log(1.0 - 
		    pnorm1(gamma(X(i,j)-1,j) - X_mean[i]) );
	  }
	  else if (X(i,j) == 1){
	    loglikerat = loglikerat 
	      + log(pnorm1(gamma_p[X(i,j)] - X_mean[i])  ) 
	      - log(pnorm1(gamma(X(i,j), j) - X_mean[i]) );
	  }
	  else{
	    loglikerat = loglikerat 
	      + log(pnorm1(gamma_p[X(i,j)] - X_mean[i]) - 
		    pnorm1(gamma_p[X(i,j)-1] - X_mean[i]) ) 
	      - log(pnorm1(gamma(X(i,j), j) - X_mean[i]) - 
		    pnorm1(gamma(X(i,j)-1, j) - X_mean[i]) );
	  }
	}
      }
      for (int k=2; k<(ncateg[j]-1); ++k){
	loggendenrat = loggendenrat 
	  + log(pnorm(gamma(k+1,j), gamma(k,j), tune[j]) - 
		pnorm(gamma(k-1,j), gamma(k,j), tune[j]) )  
	  - log(pnorm(gamma_p[k+1], gamma_p[k], tune[j]) - 
		pnorm(gamma_p[k-1], gamma_p[k], tune[j]) );
      }
      double logacceptrat = loglikerat + loggendenrat;
      if (runif() <= exp(logacceptrat)){
	for (int i=0; i<gamrow[0]; ++i){
	  if (gamma(i,j) == 300) break;
	  gamma(i,j) = gamma_p[i];
	}
	++accepts[0];
      }
    }
  



    // sample Lambda
    for (int i=0; i<K; ++i){
      Matrix<double> free_indic = t(Lambda_free_indic(i,_));
      Matrix<double> not_free_indic = (free_indic-1)*-1;
      if (sumc(free_indic)[0] > 0 && 
	  sumc(not_free_indic)[0] > 0){ // both constrnd & unconstrnd
	Matrix<double> phifree_i =  t(selif(t(phi), free_indic));
	Matrix<double> mulamfree_i = selif(t(Lambda_prior_mean(i,_)), 
					   free_indic); // prior mean
	Matrix<double> hold = selif(t(Lambda_prior_prec(i,_)), 
				    free_indic);
	Matrix<double> sig2lamfree_inv_i = 
	  eye<double>(hold.rows());   // prior prec
	for (int j=0; j<(hold.rows()); ++j)
	  sig2lamfree_inv_i(j,j) = hold[j];
	Matrix<double> Lambdacon_i = selif(t(Lambda(i,_)), not_free_indic);
	Matrix<double> phicon_i  = t(selif(t(phi), not_free_indic));
	
	Matrix<double> newX_i = Xstar(_,i) - phicon_i*Lambdacon_i; 
	Matrix<double> Lam_post_var = invpd(sig2lamfree_inv_i + 
					    t(phifree_i) * phifree_i); 
	Matrix<double> Lam_post_C = cholesky(Lam_post_var);
	Matrix<double> Lam_post_mean = Lam_post_var * 
	  (sig2lamfree_inv_i * mulamfree_i + t(phifree_i) * newX_i);
	
	Matrix<double> Lambdafree_i = 
	  gaxpy(Lam_post_C, rnorm(hold.rows(), 1), Lam_post_mean);
	
	// check to see if inequality constraints hold
	Matrix<double> Lambda_ineq_vec = Lambda_ineq(i,_);
	double ineq_holds = 0;
	int Lam_count = 0;
	for (int j=0; j<D; ++j){
	  if (free_indic[j]==1)
	    ineq_holds = std::min(ineq_holds, 
				  Lambda_ineq_vec[j]*Lambdafree_i[Lam_count]);
	}
	while (ineq_holds < 0){
	  Lambdafree_i = 
	    gaxpy(Lam_post_C, rnorm(hold.rows(), 1), Lam_post_mean);
	  Lam_count = 0;
	  double test = 0;
	  for (int j=0; j<D; ++j){
	    if (free_indic[j]==1){
	      Matrix<double> prodcheck = 
		Lambda_ineq_vec[j]*Lambdafree_i[Lam_count];    
	      test = std::min(test, prodcheck[0]); 	      
	    }
	    ++Lam_count;
	  }
	  ineq_holds = test;
	}
	
	// put draw into Lambda 
	Lam_count = 0;
	for (int j=0; j<D; ++j){
	  if (free_indic[j] == 1){
	    Lambda(i,j) = Lambdafree_i[Lam_count];
	    ++Lam_count;
	  }
	}
      }
      else  if (sumc(free_indic)[0] > 0){ // just unconstrained
	Matrix<double> phifree_i =  t(selif(t(phi), free_indic));
	Matrix<double> mulamfree_i = selif(t(Lambda_prior_mean(i,_)), 
					   free_indic); // prior mean
	Matrix<double> hold = selif(t(Lambda_prior_prec(i,_)), free_indic);
	Matrix<double> sig2lamfree_inv_i = 
	  eye<double>(hold.rows());  // prior prec
	for (int j=0; j<hold.rows(); ++j)
	  sig2lamfree_inv_i(j,j) = hold[j];
	Matrix<double> Lam_post_var = invpd(sig2lamfree_inv_i + 
					    t(phifree_i) * phifree_i); 
	Matrix<double> Lam_post_C = cholesky(Lam_post_var);
	Matrix<double> Lam_post_mean = Lam_post_var * 
	  (sig2lamfree_inv_i * mulamfree_i + t(phifree_i) * Xstar(_,i));
	Matrix<double> Lambdafree_i = 
	  gaxpy(Lam_post_C, rnorm(hold.rows(), 1), Lam_post_mean);
	
	// check to see if inequality constraints hold
	Matrix<double> Lambda_ineq_vec = Lambda_ineq(i,_);
	double ineq_holds = 0;
	int Lam_count = 0;
	for (int j=0; j<D; ++j){
	  ineq_holds = 
	    std::min(ineq_holds, Lambda_ineq_vec[j]*Lambdafree_i[j]);
	}
	while (ineq_holds < 0){
	  Lambdafree_i = 
	    gaxpy(Lam_post_C, rnorm(hold.rows(), 1), Lam_post_mean);
	  Lam_count = 0;
	  double test = 0;
	  for (int j=0; j<D; ++j){
	    double prodcheck = Lambda_ineq_vec[j]*Lambdafree_i[Lam_count];
	    test = std::min(test, prodcheck); 
	    ++Lam_count;
	  }
	  ineq_holds = test;
	}
	
	// put draw into Lambda
	Lam_count = 0;
	for (int j=0; j<D; ++j){
	  Lambda(i,j) = Lambdafree_i[Lam_count];
	  ++Lam_count;
	}
      }	
    }      
    
    
    // print results to screen
    if (verbose[0] == 1 && iter % 500 == 0 && outswitch[0] == 1){
      cout << " MCMCordfactanal Iteration = " << iter << endl;
      cout << " acceptance rate = " << static_cast<double>(accepts[0])/
	    (static_cast<double>(iter+1)*K) << endl << endl;
    }
    if (verbose[0] == 1 && iter % 500 == 0 && outswitch[0] == 2){
      cout << " MCMCirtKd Iteration = " << iter << endl;
    }
    
    // store results
    if ((iter >= burnin[0]) && ((iter % thin[0]==0))) {
      
      // store Lambda
      if (storelambda[0]==1){
	Matrix<double> Lambda_store_vec = reshape(Lambda,1,K*D);
	for (int l=0; l<K*D; ++l)
	  Lambda_store(count, l) = Lambda_store_vec[l];
      }
      
      // store gamma
      Matrix<double> gamma_store_vec = reshape(gamma, 1, gamrow[0]*gamcol[0]);
      for (int l=0; l<gamrow[0]*gamcol[0]; ++l)
	gamma_store(count, l) = gamma_store_vec[l];

      // store phi
      if (storescores[0]==1){
	Matrix<double> phi_store_vec = reshape(phi, 1, N*D);
	for (int l=0; l<N*D; ++l)
	  phi_store(count, l) = phi_store_vec[l];
      }
      count++;
    }
    

  } // end Gibbs loop
  
  
  
  // return output
  Matrix<double> output;
  if (storelambda[0]==1){
    output = cbind(Lambda_store, gamma_store);
  }
  else {
    output = gamma_store;
  }
  if(storescores[0] == 1) {
    output = cbind(output, phi_store);
  }
  int loop = samrow[0] * samcol[0];
  for (int i=0; i<loop; ++i)
    sam[i] = output[i];
  
}

}


