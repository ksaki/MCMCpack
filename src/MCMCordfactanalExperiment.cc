////////////////////////////////////////////////////////////////////////
//  Working in Progress: MCMCordfactanal for Experiment data 
//////////////////////////////////////////////////////////////////////////


#ifndef MCMCORDFACTANALEXPERIMENT_CC
#define MCMCORDFACTANALEXPERIMENT_CC

#include <iostream>

#include "matrix.h"
#include "algorithm.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

typedef Matrix<double,Row,View> rmview;

using namespace std;
using namespace scythe;

template <typename RNGTYPE>
void MCMCordfactanalExperiment_impl(rng<RNGTYPE>& stream,
        const Matrix<int>& X,
        const Matrix<int>& treatment,
        const Matrix<>& cov_phi,
        const Matrix<>& cov_tau,
        Matrix<>& Lambda,
			  Matrix<>& gamma, const Matrix<>& ncateg,
			  const Matrix<>& Lambda_eq,
			  const Matrix<>& Lambda_ineq,
			  const Matrix<>& Lambda_prior_mean,
			  const Matrix<>& Lambda_prior_prec,
        const double* tune,
			  bool storelambda, bool storescores,
			  int outswitch, unsigned int burnin,
			  unsigned int mcmc, unsigned int thin,
			  unsigned int verbose, Matrix<int>& accepts,
			  Matrix<>& output)
{

  // constants 
  const unsigned int K = X.cols();  // number of manifest variables
  const unsigned int N = X.rows();  // number of observations
  const unsigned int D = Lambda.cols();  // # of factors (incl constant)
  const unsigned int H = *max_element(treatment.begin(), treatment.end())+1; 
                    // # of treatment arms (including control)
  const unsigned int tot_iter = burnin + mcmc;  
  const unsigned int nsamp = mcmc / thin;
  const Matrix<> I = eye<double>(D-1);
  const Matrix<bool> Lambda_free_indic(K, D);
  for (unsigned int i = 0; i < (K * D); ++i) 
    if (Lambda_eq(i) == -999) 
      Lambda_free_indic(i) = true;

  const Matrix<> Psi = eye<double>(K);
  const Matrix<> Psi_inv = eye<double>(K);

  //Rprintf("Switches are %i %i %i\n", storelambda, storescores, outswitch);

  // starting values for phi, Xstar, and gamma_p
  Matrix<> phi(N, D-1);
  //Matrix<double> phi = stream->rnorm(N, D-1);
  phi = cbind(ones<double>(N,1), phi);
  Matrix<> Xstar(N, K);
  Matrix<> tau(N, H); // first col of tau should be zero (control group)
  Matrix<> coef_phi(cov_phi.cols(), 1); // coefficient of the mean of phi
  //TODO: add mixture to tau - column size should be more than 1
  Matrix<> coef_tau(cov_tau.cols(), 1); // coefficient of the mean of phi

  // storage matrices (row major order)
  Matrix<> Lambda_store;
  if (storelambda){
    Lambda_store = Matrix<double>(nsamp, K*D);
  }
  Matrix<> gamma_store(nsamp, gamma.size());
  Matrix<> phi_store;
  if (storescores){
    phi_store = Matrix<>(nsamp, N*D);
  }
  Matrix<> tau_store;
  if (storescores){
    tau_store = Matrix<>(nsamp, N*H);
  }
  Matrix<> coef_phi_store(nsamp, coef_phi.rows());
  Matrix<> coef_tau_store(nsamp, coef_tau.rows());

  ///////////////////
  // Gibbs Sampler //
  ///////////////////
  int count = 0;  
  for (unsigned int iter = 0; iter < tot_iter; ++iter) {

    // sample Xstar
    // /////////////////////////////////////////////////////////////////////////
    // Original:
    //for (unsigned int i = 0; i < N; ++i) {
    //  Matrix<> X_mean = Lambda * t(phi(i,_));
    //  for (unsigned int j = 0; j < K; ++j) {
    //    if (X(i,j) == -999) { // if missing
    //      Xstar(i,j) = stream.rnorm(X_mean[j], 1.0);
    //    } else { // if not missing
    //      Xstar(i,j) = stream.rtnorm_combo(X_mean[j], 1.0, 
    //               gamma(X(i,j)-1, j), gamma(X(i,j), j));
    //    }
    //    if (i == N-1 & j == K-1){
    //      cout << "Xstar(i,j): " << Xstar(i,j) << "\n";
    //    }
    //  }
    //}
    // /////////////////////////////////////////////////////////////////////////
    // with treatment
    //XXX: The seeds are not passed? Enabling this chunk causes randomness
    // even if the seed is provided....
    for (unsigned int i = 0; i < N; ++i) {
      Matrix<> X_mean = Lambda * t(phi(i,_));
      for (unsigned int j = 0; j < K; ++j) {
        unsigned int h = treatment(i,j);
        Matrix<> treated_part = Lambda(j,1) * tau(i,h);
        double X_mean_new = X_mean(j) + treated_part(0);
        if (X(i,j) == -999) { // if missing
          Xstar(i,j) = stream.rnorm(X_mean_new, 1.0);
        } else { // if not missing
          Xstar(i,j) = stream.rtnorm_combo(X_mean_new, 1.0, 
                   gamma(X(i,j)-1, j), gamma(X(i,j), j));
        }
        Matrix<> el = Xstar(i,j);
      }
    }
    // /////////////////////////////////////////////////////////////////////////

    // sample phi
    Matrix<> Lambda_const = Lambda(_,0);
    Matrix<> Lambda_rest = Lambda(0, 1, K-1, D-1);
    Matrix<> phi_post_var = invpd(I + crossprod(Lambda_rest) );
    Matrix<> phi_post_C = cholesky(phi_post_var);
    for (unsigned int i = 0; i < N; ++i) {
      // Original
      /////////////////////////////////////////////////////////////////////////
      //Matrix<> phi_post_mean = phi_post_var * (t(Lambda_rest)  
		  // 		       * (t(Xstar(i,_))-Lambda_const));
      /////////////////////////////////////////////////////////////////////////
      // With treatment and covariates for phi
      Matrix<> tau_obs(1, K);
      for (unsigned int j = 0; j < K; ++j) {
        unsigned int h = treatment(i,j);
        tau_obs(1,j) = tau(i,h);
      }
      Matrix<> Lambda_treat = t(Lambda_rest) * t(tau_obs);
      Matrix<> phi_post_mean = phi_post_var * (t(Lambda_rest)  
					       * (t(Xstar(i,_))-Lambda_const-Lambda_treat) + cov_phi(i,_) * coef_phi);

      Matrix<> phi_samp = gaxpy(phi_post_C, stream.rnorm(D-1, 1, 0, 1), 
				phi_post_mean);
      for (unsigned int j = 0; j < (D-1); ++j)
	      phi(i,j+1) = phi_samp(j);
      /////////////////////////////////////////////////////////////////////////
    }
    
    /////////////////////////////////////////////////////////////////////////
    // NOTE: TAU SAMPLER IS UNDER CONSTRUCTION! 
    // sample tau
    // for each treatment arm (H)
    for (unsigned int h = 0; h < H; ++h){
      // for each respondents (N)
      for (unsigned int i = 0; i < N; ++i) {
        // find j such that T_ij = h
        // XXX: What if there is no such j?
        unsigned int j_treat = 0;
        for (unsigned int j = 0; j < K; ++j) {
          if (treatment(i,j) == h){
            j_treat = j;
            break;
          }
        }
        Matrix<> Lambda_const = Lambda(j_treat, 0); // alpha_j
        // beta - submatrix from top-left (j,1) to bottom-right (j,D-1)
        // = the row of j except 0th column
        Matrix<> Lambda_rest = Lambda(j_treat, 1, j_treat, D-1); 
        Matrix<> tau_post_var = invpd(I + crossprod(Lambda_rest) );
        Matrix<> tau_post_C = cholesky(tau_post_var);

        Matrix<> Lambda_phi = t(Lambda_rest) * phi(i,1); // phi's 0th col is constant part
       
        //XXX: one column of coef_tau should be chosen when mixture is added
        Matrix<> tau_post_mean = tau_post_var * (t(Lambda_rest)  
                   * (Xstar(i,j_treat)-Lambda_const-Lambda_phi) + cov_tau(i,_) * coef_tau); 

        //cout << "dim of tau_post_C: " << tau_post_C.rows() << " " << tau_post_C.cols() << "\n"; 
        //cout << "dim of tau_post_mean: " << tau_post_mean.rows() << " " << tau_post_mean.cols() << "\n"; 

        Matrix<> tau_samp = gaxpy(tau_post_C, stream.rnorm(D-1, 1, 0, 1), 
          tau_post_mean);
        for (unsigned int j = 0; j < (D-1); ++j)
          tau(i,j+1) = tau_samp(j);
      }
    }
    /////////////////////////////////////////////////////////////////////////
				
    // sample Lambda
    /////////////////////////////////////////////////////////////////////////
    // Original:
    //NormNormfactanal_Lambda_draw(Lambda, Lambda_free_indic, 
		//		 Lambda_prior_mean, Lambda_prior_prec,
		//		 phi, Xstar, Psi_inv, Lambda_ineq, D, K,
		//		 stream);
    /////////////////////////////////////////////////////////////////////////
    // New scheme:
    // make a copy and fill?
    // This is pretty memory inefficient. Think about how to make it more efficient
    for (unsigned int j = 0; j < K; ++j){
      Matrix<> Lambda_row = Lambda(j,_);
      Matrix<> Lambda_free_indic_row = Lambda_free_indic(j,_);
      Matrix<> Lambda_prior_mean_row = Lambda_prior_mean(j,_);
      Matrix<> Lambda_prior_prec_row = Lambda_prior_prec(j,_);
      Matrix<> phi_obs = phi;
      for (unsigned int i = 0; i < N; ++i){
        unsigned int h = treatment(i,j);
        phi_obs(i,0) += tau(i,h);
      }
      Matrix<> Xstar_row = Xstar(_,j);
      Matrix<> Lambda_ineq_row = Lambda_ineq(j,_);
      NormNormfactanal_Lambda_draw(Lambda_row, Lambda_free_indic_row, 
           Lambda_prior_mean_row, Lambda_prior_prec_row,
           phi_obs, Xstar_row, Psi_inv, Lambda_ineq_row, D, 1, // set K = 1
           stream);
      Lambda(j,_) = Lambda_row;
    }
    /////////////////////////////////////////////////////////////////////////


    /////////////////////////////////////////////////////////////////////////
    // XX: New sections to be added:
    // Sample tau, eta (coef of the mean of tau) etc...
    /////////////////////////////////////////////////////////////////////////

    // sample gamma
    for (unsigned int j = 0; j < K; ++j) { 
      // do the sampling for each manifest var
      Matrix<> gamma_p = gamma(_,j);
      Matrix<> X_mean = phi * t(Lambda(j,_));
      for (unsigned int i = 2; i < (ncateg(j)); ++i) {
        if (i == (ncateg(j)-1)) {
          gamma_p(i) = stream.rtbnorm_combo(gamma(i,j), 
                    std::pow(tune[j], 2.0), gamma_p[i-1]);
        } else {
          gamma_p[i] = stream.rtnorm_combo(gamma(i,j), 
                   std::pow(tune[j], 2.0), gamma_p[i-1], gamma(i+1, j));
        }
      }
      double loglikerat = 0.0;
      double loggendenrat = 0.0;
			
			
      // loop over observations and construct the acceptance ratio
      // ///////////////////////////////////////////////////////////////////////
      // Original
      //for (unsigned int i = 0; i < N; ++i) {
      //  if (X(i,j) != -999) {
      //    if (X(i,j) == ncateg(j)) {
      //      loglikerat = loglikerat + 
      //        log(1.0  - pnorm(gamma_p[X(i,j)-1] - X_mean[i], 0, 1) ) 
      //        - log(1.0 - pnorm(gamma(X(i,j)-1,j) - X_mean[i], 0, 1) );
      //    } else if (X(i,j) == 1) { 
      //      loglikerat = loglikerat + 
      //        log(pnorm(gamma_p[X(i,j)] - X_mean[i], 0, 1)  ) 
      //        - log(pnorm(gamma(X(i,j), j) - X_mean[i], 0, 1) );
      //    } else { 
      //      loglikerat = loglikerat + 
      //        log(pnorm(gamma_p[X(i,j)] - X_mean[i], 0, 1) 
      //      - pnorm(gamma_p[X(i,j)-1] - X_mean[i], 0, 1) ) 
      //        - log(pnorm(gamma(X(i,j), j) - X_mean[i], 0, 1) - 
      //        pnorm(gamma(X(i,j)-1, j) - X_mean[i], 0, 1) ); 
      //    }
      //  }
      //}
      //
      // ///////////////////////////////////////////////////////////////////////
      // To adjust treatment effect, 
      // subtract appropriate beta_j * tau after -X_mean
      for (unsigned int i = 0; i < N; ++i) {
        unsigned int h = treatment(i,j);
        Matrix<> X_mean_treat = Lambda(j,1) * tau(i,h);
        if (X(i,j) != -999) {
          if (X(i,j) == ncateg(j)) {
            loglikerat = loglikerat + 
              log(1.0  - pnorm(gamma_p[X(i,j)-1] - X_mean[i] - X_mean_treat[0], 0, 1) ) 
              - log(1.0 - pnorm(gamma(X(i,j)-1,j) - X_mean[i] - X_mean_treat[0], 0, 1) );
          } else if (X(i,j) == 1) { 
            loglikerat = loglikerat + 
              log(pnorm(gamma_p[X(i,j)] - X_mean[i] - X_mean_treat[0], 0, 1)  ) 
              - log(pnorm(gamma(X(i,j), j) - X_mean[i] - X_mean_treat[0], 0, 1) );
          } else { 
            loglikerat = loglikerat + 
              log(pnorm(gamma_p[X(i,j)] - X_mean[i] - X_mean_treat[0], 0, 1) 
            - pnorm(gamma_p[X(i,j)-1] - X_mean[i] - X_mean_treat[0], 0, 1) ) 
              - log(pnorm(gamma(X(i,j), j) - X_mean[i] - X_mean_treat[0], 0, 1) - 
              pnorm(gamma(X(i,j)-1, j) - X_mean[i] - X_mean_treat[0], 0, 1) ); 
          }
        }
      }
      // ///////////////////////////////////////////////////////////////////////

      for (unsigned int k = 2; k < ncateg(j); ++k) {
        loggendenrat = loggendenrat 
          + log(pnorm(gamma(k+1,j), gamma(k,j), tune[j]) - 
          pnorm(gamma_p[k-1], gamma(k,j), tune[j]) )  - 
          log(pnorm(gamma_p[k+1], gamma_p[k], tune[j]) -
              pnorm(gamma(k-1,j), gamma_p[k], tune[j]) );
      }

      double logacceptrat = loglikerat + loggendenrat;
      if (stream() <= exp(logacceptrat)) { 
        for (unsigned int i = 0; i < gamma.rows(); ++i) {
          if (gamma(i,j) == 300) break;
          gamma(i,j) = gamma_p[i];
        }
        ++accepts(j);
      }
    }
	
		
    // print results to screen
    if (verbose > 0 && iter % verbose == 0 && outswitch == 1) {
      Rprintf("\n\nMCMCordfactanal iteration %i of %i \n", (iter+1),
	      tot_iter);
      Rprintf("Lambda = \n");
      for (unsigned int i = 0; i < K; ++i) {
        for (unsigned int j = 0; j < D; ++j) {
          Rprintf("%10.5f", Lambda(i,j));
	      }
	      Rprintf("\n");
      }
      Rprintf("\nMetropolis-Hastings acceptance rates = \n");
      for (unsigned int j = 0; j < K; ++j) { 
	Rprintf("%6.2f", static_cast<double>(accepts[j]) / 
		static_cast<double>((iter+1))); 
      }
    }
    if (verbose > 0 && iter % verbose == 0 && outswitch == 2) {
      Rprintf("\n\nMCMCirtKd iteration %i of %i \n", (iter+1),
	      tot_iter);
    }
		
    // store results
    if ((iter >= burnin) && ((iter % thin==0))) {      
      // store Lambda
      if (storelambda) { 
        rmview(Lambda_store(count, _)) = Lambda;
      }
			
      // store gamma
      //Matrix<> gamma_store_vec = reshape(gamma, 1, gamma.size());
      //for (unsigned int l = 0; l < gamma.size(); ++l) 
      //	gamma_store(count, l) = gamma_store_vec(l);
      rmview(gamma_store(count, _)) = gamma;

      // store phi
      if (storescores) {
      //Matrix<> phi_store_vec = reshape(phi, 1, N*D);
      //for (unsigned int l = 0; l < N * D; ++l)
      //	phi_store(count, l) = phi_store_vec(l);
        rmview(phi_store(count, _)) = phi;
      }
      
      // store tau
      if (storescores) {
      //Matrix<> phi_store_vec = reshape(phi, 1, N*D);
      //for (unsigned int l = 0; l < N * D; ++l)
      //	phi_store(count, l) = phi_store_vec(l);
        rmview(tau_store(count, _)) = tau;
      }
      count++;
    }

    // allow user interrupts
    R_CheckUserInterrupt();    

  } // end MCMC loop

  if (storelambda) {
    output = cbind(Lambda_store, gamma_store);
  } else {
    output = gamma_store;
  }
  if(storescores) {
    output = cbind(output, phi_store);
    output = cbind(output, tau_store);
  }
}

extern "C"{

  // function called by R to fit model
  void
  ordfactanalpostExperiment (double* sampledata, const int* samplerow, 
		   const int* samplecol,
		   const int* Xdata, const int* Xrow, const int* Xcol,
       const int* treatmentdata,
       const double *cov_phidata, const int* cov_phicol,
       const double *cov_taudata, const int* cov_taucol,
		   const int* burnin, const int* mcmc,  const int* thin,
		   const double* tune, const int *uselecuyer, 
		   const int *seedarray,
		   const int *lecuyerstream, const int* verbose, 
		   const double* Lamstartdata, const int* Lamstartrow, 
		   const int* Lamstartcol, 
		   const double* gamdata, const int* gamrow, const int* gamcol,
		   const int* ncatdata, const int* ncatrow, const int* ncatcol,
		   const double* Lameqdata, const int* Lameqrow, 
		   const int* Lameqcol,
		   const double* Lamineqdata, const int* Lamineqrow, 
		   const int* Lamineqcol,
		   const double* Lampmeandata, const int* Lampmeanrow, 
		   const int* Lampmeancol,
		   const double* Lampprecdata, const int* Lampprecrow,
		   const int* Lamppreccol, const int* storelambda,
		   const int* storescores,
		   int* acceptsdata, const int* acceptsrow, 
		   const int* acceptscol, const int* outswitch) 
  {

    // put together matrices
    const Matrix<int> X(*Xrow, *Xcol, Xdata);
    const Matrix<int> treatment(*Xrow, *Xcol, treatmentdata);
    const Matrix<double> cov_phi(*Xrow, *cov_phicol, cov_phidata);
    const Matrix<double> cov_tau(*Xrow, *cov_taucol, cov_taudata);
    Matrix<> Lambda(*Lamstartrow, *Lamstartcol, Lamstartdata);
    Matrix<> gamma(*gamrow, *gamcol, gamdata);
    const Matrix<> ncateg(*ncatrow, *ncatcol, ncatdata);
    const Matrix<> Lambda_eq(*Lameqrow, *Lameqcol, Lameqdata);
    const Matrix<> Lambda_ineq(*Lamineqrow, *Lamineqcol, Lamineqdata);
    const Matrix<> Lambda_prior_mean(*Lampmeanrow, *Lampmeancol, 
				     Lampmeandata);
    const Matrix<> Lambda_prior_prec(*Lampprecrow, *Lamppreccol,
				     Lampprecdata);  
    Matrix<int> accepts(*acceptsrow, *acceptscol, acceptsdata);
			
			
    // return output
    Matrix<double> output;
    MCMCPACK_PASSRNG2MODEL(MCMCordfactanalExperiment_impl,
         X, treatment, cov_phi, cov_tau, Lambda, gamma,
			   ncateg, Lambda_eq, Lambda_ineq, Lambda_prior_mean,
			   Lambda_prior_prec, tune, *storelambda, 
			   *storescores, *outswitch,
			   *burnin, *mcmc, *thin, *verbose, accepts, output);


    for (unsigned int i = 0; i < output.size(); ++i)
      sampledata[i] = output(i);

    for (unsigned int j = 0; j < X.cols(); ++j)
      acceptsdata[j] = accepts(j);
  }
}

#endif
