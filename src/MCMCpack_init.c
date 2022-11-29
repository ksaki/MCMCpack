#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void cHMMpanelFE(double *deltadraws, double* sigmadraws, double *statedraws, double* betadraws, const int* betarow, const int* betacol, const int* totalstates, const int* nsubj, const int* ntime, const int* nobs, const int* subjectid, const int* m, const int* mmax, const int* mmin, const double* Ydata, const int* Yrow, const int* Ycol, const double* Xdata, const int* Xrow, const int* Xcol, const int* burnin, const int* mcmc, const int* thin, const int* verbose, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const double* betastartdata, const double* sigma2start, const double* deltastartdata, const int* deltastartrow, const double* b0data, const double* B0data, const double* delta0, const double* Delta0, const double* c0, const double* d0, const double* P0data, const int* P0row, const double* Pstartdata, const double* subject_groupinfodata);
extern void cHMMpanelRE(double* betadata, const int* betarow, const int* betacol, double* sigmadata, double* Ddata, double *psout, double *sout, const int* nsubj, const int* ntime, const int* m, const int* nobs, const int* subjectid, const int* timeid, const double* Ydata, const int* Yrow, const int* Ycol, const double* Xdata, const int* Xrow, const int* Xcol, const double* Wdata, const int* Wrow, const int* Wcol, const double* YTdata, const double* XTdata, const double* WTdata, const int* burnin, const int* mcmc, const int* thin, const int* verbose, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const double* betastartdata, const double* sigma2start, const double *Pstart, const double* b0data, const double* B0data, const double* c0, const double* d0, const int* r0, const double* R0data, const double* subject_groupinfodata, const double* time_groupinfodata, double *logmarglikeholder, double *loglikeholder, const int *chib);
extern void cMCMCbinaryChange(double *phiout, double *Pout, double *psout, double *sout, const double *Ydata, const int *Yrow, const int *Ycol, const int *m, const int *burnin, const int *mcmc, const int *thin, const int *verbose, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const double *phistart, const double *Pstart, const double *a, const double *b, const double *c0, const double *d0, const double *A0data, double *logmarglikeholder, const int *chib);
extern void cMCMCdynamicIRT1d_b(double* thetadraws, const int* nrowthetadraws, const int* ncolthetadraws, double* alphadraws, const int* nrowalphadraws, const int* ncolalphadraws, double* betadraws, const int* nrowbetadraws, const int* ncolbetadraws, double* tau2draws, const int* nrowtau2draws, const int* ncoltau2draws, const int* nsubj, const int* nitems, const int* ntime, const int* Ydata, const int* nrowYdata, const int* ncolYdata, const int* ITdata, const int* lengthITdata, const int* burnin, const int* mcmc, const int* thin, const int* uselecuyer, const int* seedarray, const int* lecuyerstream, const int* verbose, const double* thetastartdata, const int* lengththetastart, const int* thetainfodata, const int* nrowthetainfo, const int* ncolthetainfo, double* alphastartdata, const int* lengthalphastart, double* betastartdata, const int* lengthbetastart, double* tau2startdata, const int* lengthtau2start, const double* c0, const int* lengthc0, const double* d0, const int* lengthd0, const double* a0, double* A0, const double* b0, double* B0, const double* e0, const double* E0inv, const double* thetaeqdata, const int* nrowthetaeq, const int* ncolthetaeq, const double* thetaineqdata, const int* nrowthetaineq, const int* ncolthetaineq, const int* storeitem, const int* storeability);
extern void cMCMCdynamicIRT1d(double* thetadraws, const int* nrowthetadraws, const int* ncolthetadraws, double* alphadraws, const int* nrowalphadraws, const int* ncolalphadraws, double* betadraws, const int* nrowbetadraws, const int* ncolbetadraws, double* tau2draws, const int* nrowtau2draws, const int* ncoltau2draws, const int* nsubj, const int* nitems, const int* ntime, const int* Ydata, const int* nrowYdata, const int* ncolYdata, const int* ITdata, const int* lengthITdata, const int* burnin, const int* mcmc, const int* thin, const int* uselecuyer, const int* seedarray, const int* lecuyerstream, const int* verbose, const double* thetastartdata, const int* lengththetastart, const int* thetainfodata, const int* nrowthetainfo, const int* ncolthetainfo, double* alphastartdata, const int* lengthalphastart, double* betastartdata, const int* lengthbetastart, double* tau2startdata, const int* lengthtau2start, const double* c0, const int* lengthc0, const double* d0, const int* lengthd0, const double* a0, const double* A0, const double* b0, const double* B0, const double* e0, const double* E0inv, const double* thetaeqdata, const int* nrowthetaeq, const int* ncolthetaeq, const double* thetaineqdata, const int* nrowthetaineq, const int* ncolthetaineq, const int* storeitem, const int* storeability);
extern void cMCMChlogit(const int *ngibbs, const int *nthin, const int *nburn, const int *nobs, const int *ngroup, const int *np, const int *nq, const int *IdentGroup, const double *Y_vect, const double *X_vect, const double *W_vect, double *beta_vect, double *b_vect, double *Vb_vect, double *V, const double *mubeta_vect, const double *Vbeta_vect, const double *r, const double *R_vect, const double *s1_V, const double *s2_V, double *Deviance, double *theta_pred, const int *seed, const int *verbose, const int *FixOD);
extern void cMCMChpoisson(const int *ngibbs, const int *nthin, const int *nburn, const int *nobs, const int *ngroup, const int *np, const int *nq, const int *IdentGroup, const double *Y_vect, const double *X_vect, const double *W_vect, double *beta_vect, double *b_vect, double *Vb_vect, double *V, const double *mubeta_vect, const double *Vbeta_vect, const double *r, const double *R_vect, const double *s1_V, const double *s2_V, double *Deviance, double *lambda_pred, const int *seed, const int *verbose, const int *FixOD);
extern void cMCMChregress(const int *ngibbs, const int *nthin, const int *nburn, const int *nobs, const int *ngroup, const int *np, const int *nq, const int *IdentGroup, const double *Y_vect, const double *X_vect, const double *W_vect, double *beta_vect, double *b_vect, double *Vb_vect, double *V, const double *mubeta_vect, const double *Vbeta_vect, const double *r, const double *R_vect, const double *s1_V, const double *s2_V, double *Deviance, double *Y_pred, const int *seed, const int *verbose);
extern void cMCMCregress(double *sampledata, const int *samplerow, const int *samplecol, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *burnin, const int *mcmc, const int *thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *betastartdata, const int *betastartrow, const int *betastartcol, const double *b0data, const int *b0row, const int *b0col, const double *B0data, const int *B0row, const int *B0col, const double *c0, const double *d0, double* logmarglikeholder, const int* chib);
extern void cMCMCirt1d(double* sampledata, const int* samplerow, const int* samplecol, const int* Xdata, const int* Xrow, const int* Xcol, const int* burnin, const int* mcmc,  const int* thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int* verbose, const double* thetastartdata, const int* thetastartrow, const int* thetastartcol, const double* astartdata, const int* astartrow, const int* astartcol, const double* bstartdata, const int* bstartrow, const int* bstartcol, const double* t0, const double* T0,	const double* ab0data, const int* ab0row, const int* ab0col, const double* AB0data, const int* AB0row, const int* AB0col, const double* thetaeqdata, const int* thetaeqrow, const int* thetaeqcol, const double* thetaineqdata, const int* thetaineqrow, const int* thetaineqcol, const int* storei, const int* storea);
extern void cMCMCirtHier1d(double* sampledata, const int* samplerow, const int* samplecol, const int* Xdata, const int* Xrow, const int* Xcol, const int* burnin, const int* mcmc,  const int* thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int* verbose, const double* thetastartdata, const int* thetastartrow, const int* thetastartcol, const double* astartdata, const int* astartrow, const int* astartcol, const double* bstartdata, const int* bstartrow, const int* bstartcol, const double* ab0data, const int* ab0row, const int* ab0col, const double* AB0data, const int* AB0row, const int* AB0col, const double* Xjdata, const int* Xjrow, const int* Xjcol, const double* betastartdata, const int* betastartrow, const int* betastartcol, const double* b0data, const int* b0row, const int* b0col, const double* B0data, const int* B0row, const int* B0col, const double* c0, const double* d0, const int* storei, const int* storea, double* logmarglikeholder, const int* chib, const int* px, const double* px_a0, const double* px_b0);
extern void cMCMCoprobitChange(double *betaout, double *betalinearout, double *gammaout, double *Pout, double *psout, double *sout, const double *Ydata, const double *Xdata, const int *Xrow, const int *Xcol, const int *m, const int *ncat, const int *burnin, const int *mcmc, const int *thin, const int *verbose, const double *tunedata, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const double *betastart,  const double *betalinearstart, const double *gammastart, const double *Pstart, const double *sigmastart, const double *a, const double *b, const double *b0data, const double *B0data, const double *A0data, double *logmarglikeholder, double *loglikeholder, const int *chib, const int *gammafixed);
extern void cMCMCpoissonChange(double *betaout, double *Pout, double *psout, double *sout, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *m, const int *burnin, const int *mcmc, const int *thin, const int *verbose, const double *betastart, const double *Pstart, const double *taustart, const double *componentstart, const double *a, const double *b, const double *c0, const double *d0, const int* uselecuyer, const int* seedarray, const int* lecuyerstream, const double *b0data, const double *B0data, const double *A0data, double *logmarglikeholder, double *loglikeholder, const double *wrin, const double *mrin, const double *srin, const int *chib);
extern void cMCMCprobitChange(double *betaout, double *Pout, double *psout, double *sout, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *m, const int *burnin, const int *mcmc, const int *thin, const int *verbose, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const double *betastart,  const double *Pstart, const double *a, const double *b, const double *b0data, const double *B0data, const double *A0data, double *logmarglikeholder, double *loglikeholder, const int *chib);
extern void cMCMCregressChange(double *betaout, double *Sigmaout, double *psout, double *sout,  double *yloglike, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *m, const int *burnin, const int *mcmc, const int *thin, const int *verbose, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const double *betastart, const double *Sigmastart, const double *Pstart, const int *statestart, const double *a, const double *b, const double *b0data, const double *B0data, const double *c0, const double *d0, const double *A0data, double *logmarglikeholder, double *loglikeholder, const int *marginalrun, const int *sos);
extern void cMCMCresidualBreakAnalysis(double *betaout, double *Sigmaout, double *psout,  double *sout, double *yloglike, const double *Ydata, const int *Yrow, const int *Ycol, const int *m, const int *burnin, const int *mcmc, const int *thin, const int *verbose, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const double *betastart, const double *Sigmastart, const double *Pstart, const int *statestart, const double *a, const double *b, const double *b0data, const double *B0data, const double *c0, const double *d0, const double *A0data, double *logmarglikeholder, double *loglikeholder, const int *marginalrun, const int *sos);
extern void cMCMCSVDreg(double *sampledata, const int *samplerow, const int *samplecol, const double *Ydata, const int *Yrow, const int *Ycol, const int *Ymiss, const double *Adata, const int *Arow, const int *Acol, const double *Ddata, const int *Drow, const int *Dcol, const double *Fdata, const int *Frow, const int *Fcol, const int *burnin, const int *mcmc, const int *thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *taustartdata, const int *taustartrow, const int *taustartcol, const double *g0data, const int *g0row, const int *g0col, const double *a0, const double *b0, const double* c0, const double* d0, const double* w0, const int* betasamp);
extern void dynamicEI(double* sample, const int* samrow, const int* samcol, const double* Rr0, const double* Rr1, const double* Rc0, const double* Rc1, const int* Rntables, const int* Rburnin, const int* Rmcmc, const int* Rthin, const double* RW, const double* Rnu0, const double* Rdelta0, const double* Rnu1, const double* Rdelta1, const int* Rverbose, const int *uselecuyer, const int *seedarray, const int *lecuyerstream);
extern void hierBetaBinom(double* sampledata, const int* samplerow, const int* samplecol, const int* y, const int* s, const double* theta_start, const double* alpha_start, const double* beta_start, const double* a, const double* b, const int* ilabels, const int* jlabels, const int* ilabelsunique, const int* jlabelsunique, const int* n, const int* ni, const int* nj, const int* burnin, const int* mcmc,  const int* thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int* verbose, int *accepts, const double* base_sigma);
extern void hierEI(double* sample, const int* samrow, const int* samcol, const double* Rr0, const double* Rr1, const double* Rc0, const double* Rc1, const int* Rntables, const int* Rburnin, const int* Rmcmc, const int* Rthin, const double* Rmu0pm, const double* Rmu0pv, const double* Rmu1pm, const double* Rmu1pv, const double* Rnu0, const double* Rdelta0, const double* Rnu1, const double* Rdelta1, const int* Rverbose, const int *uselecuyer, const int *seedarray, const int *lecuyerstream);
extern void HMMmultivariateGaussian(double* betadata, const int* betarow, const int* betacol, double* sigmadata, double *psout, const int* nsubj, const int* ntime, const int* m, const int* nobs, const int* subjectid, const int* timeid, const double* Ydata, const int* Yrow, const int* Ycol, const double* Xdata, const int* Xrow, const int* Xcol, const double* YTdata, const double* XTdata, const int* burnin, const int* mcmc, const int* thin, const int* verbose, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const double* betastartdata, const double* sigma2start, const double* b0data, const double* B0data, const double* c0, const double* d0, const double* P0data, const int* P0row, const int* P0col, const double* subject_groupinfodata, const double* time_groupinfodata, double *logmarglikeholder, double *loglikeholder, const int *chib);
extern void irtKdHetpost(double *samdata, const int *samrow, const int *samcol, const int *Xdata, const int *Xrow, const int *Xcol, const int *burnin, const int *mcmc, const int *thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *Lamstartdata, const int *Lamstartrow, const int *Lamstartcol, const double *Lameqdata, const int *Lameqrow, const int *Lameqcol, const double *Lamineqdata, const int *Lamineqrow, const int *Lamineqcol, const double *Lampmeandata, const int *Lampmeanrow, const int *Lampmeancol, const double *Lampprecdata, const int *Lampprecrow, const int *Lamppreccol, const int *storelambda, const int *storescores, const int *storesigma, const double *sigmapriorc, const double *sigmapriord);
extern void irtKdRobpost(double* sampledata, const int* samplerow, const int* samplecol, const int* Xdata, const int* Xrow, const int* Xcol, const int* burnin, const int* mcmc,  const int* thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int* verbose, const int* method_step, const double* theta_w, const int* theta_p, const double* lambda_w, const int* lambda_p, const double* delta0_w, const int* delta0_p, const double* delta1_w, const int* delta1_p, const double * delta0start, const double* delta1start, const double* Lamstartdata, const int* Lamstartrow, const int* Lamstartcol, const double* thetstartdata, const int* thetstartrow, const int* thetstartcol, const double* Lameqdata, const int* Lameqrow, const int* Lameqcol, const double* Lamineqdata, const int* Lamineqrow, const int* Lamineqcol, const double* theteqdata, const int* theteqrow, const int* theteqcol, const double* thetineqdata, const int* thetineqrow, const int* thetineqcol, const double* Lampmeandata, const int* Lampmeanrow, const int* Lampmeancol, const double* Lampprecdata, const int* Lampprecrow, const int* Lamppreccol, const double* k0, const double* k1, const double* c0, const double* c1, const double* d0, const double* d1, const int* storeitem, const int* storeability);
extern void mixfactanalpost(double* sampledata, const int* samplerow, const int* samplecol, const double* Xdata, const int* Xrow, const int* Xcol, const int* burnin, const int* mcmc,  const int* thin, const double* tune, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int* verbose, const double* Lamstartdata, const int* Lamstartrow, const int* Lamstartcol, const double* gamdata, const int* gamrow, const int* gamcol, const double* Psistartdata, const int* Psistartrow, const int* Psistartcol, const int* ncatdata, const int* ncatrow, const int* ncatcol, const double* Lameqdata, const int* Lameqrow, const int* Lameqcol, const double* Lamineqdata, const int* Lamineqrow, const int* Lamineqcol, const double* Lampmeandata, const int* Lampmeanrow, const int* Lampmeancol, const double* Lampprecdata, const int* Lampprecrow, const int* Lamppreccol, const double* a0data, const int* a0row, const int* a0col, const double* b0data, const int* b0row, const int* b0col, const int* storelambda, const int* storescores, int* acceptsdata, const int* acceptsrow, const int* acceptscol);
extern void ordfactanalpost(double* sampledata, const int* samplerow, const int* samplecol, const int* Xdata, const int* Xrow, const int* Xcol, const int* burnin, const int* mcmc,  const int* thin, const double* tune, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int* verbose, const double* Lamstartdata, const int* Lamstartrow, const int* Lamstartcol, const double* gamdata, const int* gamrow, const int* gamcol, const int* ncatdata, const int* ncatrow, const int* ncatcol, const double* Lameqdata, const int* Lameqrow, const int* Lameqcol, const double* Lamineqdata, const int* Lamineqrow, const int* Lamineqcol, const double* Lampmeandata, const int* Lampmeanrow, const int* Lampmeancol, const double* Lampprecdata, const int* Lampprecrow, const int* Lamppreccol, const int* storelambda, const int* storescores, int* acceptsdata, const int* acceptsrow, const int* acceptscol, const int* outswitch);
extern void ordfactanalpostExperiment(double* sampledata, const int* samplerow, const int* samplecol, const int* Xdata, const int* Xrow, const int* Xcol, const int* treatment, const int* burnin, const int* mcmc,  const int* thin, const double* tune, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int* verbose, const double* Lamstartdata, const int* Lamstartrow, const int* Lamstartcol, const double* gamdata, const int* gamrow, const int* gamcol, const int* ncatdata, const int* ncatrow, const int* ncatcol, const double* Lameqdata, const int* Lameqrow, const int* Lameqcol, const double* Lamineqdata, const int* Lamineqrow, const int* Lamineqcol, const double* Lampmeandata, const int* Lampmeanrow, const int* Lampmeancol, const double* Lampprecdata, const int* Lampprecrow, const int* Lamppreccol, const int* storelambda, const int* storescores, int* acceptsdata, const int* acceptsrow, const int* acceptscol, const int* outswitch);
extern void cMCMCfactanal(double *sampledata, const int *samplerow, const int *samplecol, const double *Xdata, const int *Xrow, const int *Xcol, const int *burnin, const int *mcmc, const int *thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *Lambdadata, const int *Lambdarow, const int *Lambdacol, const double *Psidata, const int *Psirow, const int *Psicol, const double *Lameqdata, const int *Lameqrow, const int *Lameqcol, const double *Lamineqdata, const int *Lamineqrow, const int *Lamineqcol, const double *Lampmeandata, const int *Lampmeanrow, const int *Lampmeancol, const double *Lampprecdata, const int *Lampprecrow, const int *Lamppreccol, const double *a0data, const int *a0row, const int *a0col, const double *b0data, const int *b0row, const int *b0col, const int *storescores);
extern void cMCMClogit(double *sampledata, const int *samplerow, const int *samplecol, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *burnin, const int *mcmc, const int *thin, const double *tunedata, const int *tunerow, const int *tunecol, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *betastartdata, const int *betastartrow, const int *betastartcol, const double *b0data, const int *b0row, const int *b0col, const double *B0data, const int *B0row, const int *B0col, const double *Vdata, const int *Vrow, const int *Vcol);
extern void MCMCmnlMH(double *sampledata, const int *samplerow, const int *samplecol, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *burnin, const int *mcmc, const int *thin, const double *tunedata, const int *tunerow, const int *tunecol, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *betastartdata, const int *betastartrow, const int *betastartcol, const double *betamodedata, const int *betamoderow, const int *betamodecol, const double *b0data, const int *b0row, const int *b0col, const double *B0data, const int *B0row, const int *B0col, const double *Vdata, const int *Vrow, const int *Vcol, const int* RW, const double* tdf);
extern void MCMCmnlslice(double *sampledata, const int *samplerow, const int *samplecol, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *burnin, const int *mcmc, const int *thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *betastartdata, const int *betastartrow, const int *betastartcol, const double *b0data, const int *b0row, const int *b0col, const double *B0data, const int *B0row, const int *B0col, const double *Vdata, const int *Vrow, const int *Vcol);
extern void cMCMCoprobit(double *sampledata, const int *samplerow, const int *samplecol, const int *Y, const double *nYdata, const int *nYrow, const int *nYcol, const double *Xdata, const int *Xrow, const int *Xcol, const int *burnin, const int *mcmc, const int *thin, const double *tunedata, const int *tunerow, const int *tunecol, const double* tdf, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *betadata, const int *betarow, const int *betacol, const double* gammadata, const int* gammarow, const int* gammacol, const double *b0data, const int *b0row, const int *b0col, const double *B0data, const int *B0row, const int *B0col, const double *a0data, const int *a0row, const int *a0col, const double *A0data, const int *A0row, const int *A0col, const int *cowles);
extern void cMCMCpoisson(double *sampledata, const int *samplerow, const int *samplecol, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *burnin, const int *mcmc, const int *thin, const double *tunedata, const int *tunerow, const int *tunecol, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *betastartdata, const int *betastartrow, const int *betastartcol, const double *b0data, const int *b0row, const int *b0col, const double *B0data, const int *B0row, const int *B0col, const double *Vdata, const int *Vrow, const int *Vcol);
extern void cMCMCprobit(double *sampledata, const int *samplerow, const int *samplecol, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *burnin, const int *mcmc, const int *thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *betastartdata, const int *betastartrow, const int *betastartcol, const double *b0data, const int *b0row, const int *b0col, const double *B0data, const int *B0row, const int *B0col, double *logmarglikeholder, const int *chib);
extern void MCMCprobitres(double *sampledata, const int *samplerow, const int *samplecol, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const double *resvecdata, const int *resvecrow, const int *resveccol, const int *burnin, const int *mcmc, const int *thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *betastartdata, const int *betastartrow, const int *betastartcol, const double *b0data, const int *b0row, const int *b0col, const double *B0data, const int *B0row, const int *B0col, double *logmarglikeholder, const int *chib);
extern void cMCMCquantreg(double *sampledata, const int *samplerow, const int *samplecol, const double *tau, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *burnin, const int *mcmc, const int *thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *betastartdata, const int *betastartrow, const int *betastartcol, const double *b0data, const int *b0row, const int *b0col, const double *B0data, const int *B0row, const int *B0col);
extern void cMCMCtobit(double *sampledata, const int *samplerow, const int *samplecol, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const double *below, const double *above, const int *burnin, const int *mcmc, const int *thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *betastartdata, const int *betastartrow, const int *betastartcol, const double *b0data, const int *b0row, const int *b0col, const double *B0data, const int *B0row, const int *B0col, const double *c0, const double *d0);
extern void cSSVSquantreg(double *sampledata, const int *samplerow, const int *samplecol, const double *tau, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *q, const int *burnin, const int *mcmc, const int *thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int *verbose, const double *pi0a0, const double *pi0b0);
extern void cMCMCnegbin(double *betaout, double *nuout, double *rhoout, double *tau1out, double *tau2out, int *comp1out, int *comp2out, double *sr1out, double *sr2out, double *mr1out, double *mr2out, double *rhosizes, const double *Ydata,  const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *burnin,  const int *mcmc, const int *thin, const int *verbose,  const double *betastart, const double *nustart, const double *rhostart, const double *tau1start, const double *tau2start, const double *component1start, const double *e, const double *f, const double *g, const double *rhostepdata, const int* uselecuyer, const int* seedarray,  const int* lecuyerstream, const double *b0data,  const double *B0data, double *logmarglikeholder, double *loglikeholder, const int *chib);
extern void cMCMCnegbinChange(double *betaout, double *Pout, double *psout, double *sout,  double *nuout, double *rhoout, double *tau1out, double *tau2out, int *comp1out, int *comp2out, double *sr1out, double *sr2out, double *mr1out, double *mr2out, double *rhosizes, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *m, const int *burnin, const int *mcmc, const int *thin, const int *verbose, const double *betastart, const double *Pstart, const double *nustart, const double *rhostart, const double *tau1start, const double *tau2start, const double *component1start, const double *a, const double *b, const double *e, const double *f, const double *g, const double *rhostepdata, const int* uselecuyer, const int* seedarray, const int* lecuyerstream, const double *b0data, const double *B0data, const double *A0data, const int *fixed_m, double *logmarglikeholder, double *loglikeholder, const int *chib);
void cHDPHMMnegbin(double *betaout, double *Pout, double *psout, double *sout, double *nuout, double *rhoout, double *tau1out, double *tau2out, int *comp1out, int *comp2out, double *sr1out, double *sr2out, double *mr1out, double *mr2out, double *gammaout, double *akout, double *thetaout, double *rhosizes, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *K, const int *burnin, const int *mcmc, const int *thin, const int *verbose, const double *betastart, const double *Pstart, const double *nustart, const double *rhostart, const double *tau1start, const double *tau2start, const double *component1start, const double *gammastart, const double *akstart, const double *thetastart, const double *a_alpha, const double *b_alpha, const double *a_gamma, const double *b_gamma, const double *a_theta, const double *b_theta, const double *e, const double *f, const double *g, const double *rhostepdata, const int* uselecuyer, const int* seedarray, const int* lecuyerstream, const double *b0data, const double *B0data);
extern void cHDPHMMpoisson(double *betaout, double *Pout, double *psout, double *sout, double *tau1out, double *tau2out, int *comp1out, int *comp2out, double *sr1out, double *sr2out, double *mr1out, double *mr2out, double *gammaout, double *akout, double *thetaout, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata,const int *Xrow, const int *Xcol, const int *K, const int *burnin, const int *mcmc, const int *thin, const int *verbose, const double *betastart, const double *Pstart, const double *tau1start, const double *tau2start, const double *component1start, const double *gammastart, const double *akstart, const double *thetastart, const double *a_alpha, const double *b_alpha, const double *a_gamma, const double *b_gamma, const double *a_theta, const double *b_theta, const int* uselecuyer, const int* seedarray, const int* lecuyerstream, const double *b0data, const double *B0data);
void cHDPHSMMnegbin(double *betaout, double *Pout, double *omegaout, double *sout,  double *nuout, double *rhoout, double *tau1out, double *tau2out, int *comp1out, int *comp2out, double *sr1out, double *sr2out, double *mr1out, double *mr2out, double *gammaout, double *alphaout, double *rhosizes, const double *Ydata, const int *Yrow, const int *Ycol, const double *Xdata, const int *Xrow, const int *Xcol, const int *K, const int *burnin, const int *mcmc, const int *thin, const int *verbose, const double *betastart, const double *Pstart, const double *nustart, const double *rhostart, const double *tau1start, const double *tau2start, const double *component1start, const double *alphastart, const double *gammastart, const double *omegastart, const double *a_alpha, const double *b_alpha, const double *a_gamma, const double *b_gamma, const double *a_omega, const double *b_omega, const double *e, const double *f, const double *g, const double *r, const double *rhostepdata, const int* uselecuyer, const int* seedarray, const int* lecuyerstream, const double *b0data, const double *B0data);
void cMCMCpaircompare(double* sampledata, const int* samplerow, const int* samplecol, const unsigned int* MDdata, const int* MDrow, const int* MDcol, const int* alphafixed, const int* burnin, const int* mcmc, const int* thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int* verbose, const double* thetastartdata, const int* thetastartrow, const int* thetastartcol, const double* astartdata, const int* astartrow, const int* astartcol, const double* a0, const double* A0, const double* thetaeqdata, const int* thetaeqrow, const int* thetaeqcol, const double* thetaineqdata, const int* thetaineqrow, const int* thetaineqcol, const int* storealpha, const int* storetheta);
void cMCMCpaircompare2d(double* sampledata, const int* samplerow, const int* samplecol, const unsigned int* MDdata, const int* MDrow, const int* MDcol, const int* burnin, const int* mcmc, const int* thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int* verbose, const double* thetastartdata, const int* thetastartrow, const int* thetastartcol, const double* gammastartdata, const int* gammastartrow, const int* gammastartcol, const double* tunevalue, const double* thetaeqdata, const int* thetaeqrow, const int* thetaeqcol, const double* thetaineqdata, const int* thetaineqrow, const int* thetaineqcol, const int* storegamma, const int* storetheta, double* gammaacceptrate);
void cMCMCpaircompare2dDP(double* sampledata, const int* samplerow, const int* samplecol, const unsigned int* MDdata, const int* MDrow, const int* MDcol, const int* burnin, const int* mcmc, const int* clustermcmc, const int* thin, const int *uselecuyer, const int *seedarray, const int *lecuyerstream, const int* verbose, const double* thetastartdata, const int* thetastartrow, const int* thetastartcol, const double* gammastartdata, const int* gammastartrow, const int* gammastartcol, const double* clustergammastartdata, const int* clustergammastartrow, const int* clustergammastartcol, const int* judgeclustermembershipstartdata, const int* judgeclustermembershipstartrow, const int* judgeclustermembershipstartcol, const double* tunevalue, const double* thetaeqdata, const int* thetaeqrow, const int* thetaeqcol, const double* thetaineqdata, const int* thetaineqrow, const int* thetaineqcol, const int* storegamma, const int* storetheta, double* gammaacceptrate, double* alpha, const unsigned int* clustermax, const int* alphafixed, const double* a, const double* b);




/* .Call calls */
extern SEXP MCMClogituserprior_cc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MCMCmetrop1R_cc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"cHMMpanelFE",                (DL_FUNC) &cHMMpanelFE,                41},
  {"cHMMpanelRE",                (DL_FUNC) &cHMMpanelRE,                46},
  {"cMCMCbinaryChange",          (DL_FUNC) &cMCMCbinaryChange,          24},
  {"cMCMCdynamicIRT1d",          (DL_FUNC) &cMCMCdynamicIRT1d,          56},
  {"cMCMCdynamicIRT1d_b",        (DL_FUNC) &cMCMCdynamicIRT1d_b,        56},
  {"cMCMChlogit",                (DL_FUNC) &cMCMChlogit,                26},
  {"cMCMChpoisson",              (DL_FUNC) &cMCMChpoisson,              26},
  {"cMCMChregress",              (DL_FUNC) &cMCMChregress,              25},
  {"cMCMCregress",               (DL_FUNC) &cMCMCregress,               29},
  {"cMCMCirt1d",                 (DL_FUNC) &cMCMCirt1d,                 38},
  {"cMCMCirtHier1d",             (DL_FUNC) &cMCMCirtHier1d,             49},
  {"cMCMCoprobitChange",         (DL_FUNC) &cMCMCoprobitChange,         34},
  {"cMCMCpoissonChange",         (DL_FUNC) &cMCMCpoissonChange,         35},
  {"cMCMCprobitChange",          (DL_FUNC) &cMCMCprobitChange,          28},
  {"cMCMCregressChange",         (DL_FUNC) &cMCMCregressChange,         34},
  {"cMCMCresidualBreakAnalysis", (DL_FUNC) &cMCMCresidualBreakAnalysis, 31},
  {"cMCMCSVDreg",                (DL_FUNC) &cMCMCSVDreg,                35},
  {"dynamicEI",                  (DL_FUNC) &dynamicEI,                  20},
  {"hierBetaBinom",              (DL_FUNC) &hierBetaBinom,              26},
  {"hierEI",                     (DL_FUNC) &hierEI,                     23},
  {"HMMmultivariateGaussian",    (DL_FUNC) &HMMmultivariateGaussian,    40},
  {"irtKdHetpost",               (DL_FUNC) &irtKdHetpost,               34},
  {"irtKdRobpost",               (DL_FUNC) &irtKdRobpost,               56},
  {"mixfactanalpost",            (DL_FUNC) &mixfactanalpost,            49},
  {"ordfactanalpost",            (DL_FUNC) &ordfactanalpost,            41},
  {"ordfactanalpostExperiment",  (DL_FUNC) &ordfactanalpostExperiment,  42},
  {"cMCMCfactanal",              (DL_FUNC) &cMCMCfactanal,              38},
  {"cMCMClogit",                 (DL_FUNC) &cMCMClogit,                 31},
  {"MCMCmnlMH",                  (DL_FUNC) &MCMCmnlMH,                  36},
  {"MCMCmnlslice",               (DL_FUNC) &MCMCmnlslice,               28},
  {"cMCMCoprobit",               (DL_FUNC) &cMCMCoprobit,               40},
  {"cMCMCpoisson",               (DL_FUNC) &cMCMCpoisson,               31},
  {"cMCMCprobit",                (DL_FUNC) &cMCMCprobit,                27},
  {"MCMCprobitres",              (DL_FUNC) &MCMCprobitres,              30},
  {"cMCMCquantreg",              (DL_FUNC) &cMCMCquantreg,              26},
  {"cMCMCtobit",                 (DL_FUNC) &cMCMCtobit,                 29},
  {"cSSVSquantreg",              (DL_FUNC) &cSSVSquantreg,              20},
  {"cMCMCnegbin",                (DL_FUNC) &cMCMCnegbin,                40},
  {"cMCMCnegbinChange",          (DL_FUNC) &cMCMCnegbinChange,          49},
  {"cHDPHMMnegbin",              (DL_FUNC) &cHDPHMMnegbin,              54},
  {"cHDPHMMpoisson",             (DL_FUNC) &cHDPHMMpoisson,             45},
  {"cHDPHSMMnegbin",             (DL_FUNC) &cHDPHSMMnegbin,             54},
  {"cMCMCpaircompare",           (DL_FUNC) &cMCMCpaircompare,           30},
  {"cMCMCpaircompare2d",         (DL_FUNC) &cMCMCpaircompare2d,         29},
  {"cMCMCpaircompare2dDP",       (DL_FUNC) &cMCMCpaircompare2dDP,       41},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"MCMClogituserprior_cc", (DL_FUNC) &MCMClogituserprior_cc, 14},
  {"MCMCmetrop1R_cc",       (DL_FUNC) &MCMCmetrop1R_cc,       12},
  {NULL, NULL, 0}
};

void R_init_MCMCpack(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
