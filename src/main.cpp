#include <RcppArmadillo.h>
#include <truncnorm.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

/* *********** NOTES *********** 
 * X = solve( A, B ); Solve a dense system of linear equations, A*X = B, where X is unknown
 * i.e., solve( A, B ) gives A^{-1}B
 * 
 */

/* *********** GENERAL/HELPER FUNCTIONS *********** 
 * 
 * Used for creating matrices and such needed by the samplers below.
 * Includes a multivariate normal sampler, etc.
 * 
 */

// [[Rcpp::export]]
double rtruncnorm(const double mu, const double sigma, const double a,
                  const double b){
  // See https://cran.r-project.org/web/packages/RcppDist/vignettes/RcppDist.pdf
  return r_truncnorm(mu, sigma, a, b);
}

// [[Rcpp::export]]
arma::mat rmvnrnd(arma::vec mu, arma::mat Sigma, int n){
  return arma::mvnrnd(mu, Sigma, n);
}

// [[Rcpp::export]]
double rnormArma(double mu, double sig_sq) {
  /* A single draw of a univariate normal distribution,
  */ 
  return mu + pow(sig_sq, 0.5) * arma::randn( );
}

// [[Rcpp::export]]
arma::mat get_Omega_i(arma::mat Theta_c, arma::mat xi_i) {
  /* 
  * Return P x K \Omega_i = \Theta_{c[i]} \xi_{c[i]}(x_i) for i in 1,...,N_c.
  * Input Theta_c is P x L coefficient matrix associated with cause c[i] ;
  * and xi_i is L x K loadings matrix \xi_{c[i]}(x_i) associated with death i.
  * 
  */
  return Theta_c * xi_i;
}

/* *********** MAIN TERMS SAMPLER FUNCTIONS *********** 
 * 
 * Used for sampling z_{ij}, \xi_{c,lk}, \psi_{c,k}, 
 * \nu_{i}, \theta_{c,lk}, \delta_{lk}, and \sigma_j^2.
 * 
 */

// [[Rcpp::export]]
double sample_xi_ilk(double sig_sq_xi, arma::vec beta_c, arma::vec X_i,
                     arma::vec Sigma_0_vec, arma::mat Theta_c, arma::mat xi_i, 
                     arma::vec z_i, arma::vec eta_i, int l_samp, int k_samp, int fix_xi) {
  /* 
  * Sample \xi_{c,lk}(x_i) where the result is a double for person i having COD c.
  * 
  * The inputs are comprised of "prior" and "likelihood" components;
  * "Prior" components:
  * sig_sq_xi, the variance term associated with obs of \xi_i ; 
  * X_i, the B x 1 covariate vector (X_{i,*})^T for death i with COD c ; 
  * beta_c is the B x 1 cause-specific B vector of coefficients for the regression of \xi ;
  * "Likelihood" components:
  * Sigma_0_vec is the P-vector diag(\Sigma_0) ;
  * Theta_c is the P x L matrix \Theta specific to COD c person i died of ; 
  * xi_i is the current L x K matrix \Xi person i dying of cod c ; 
  * z_i, the P*1 vector [z_{i1},...,z_{iP}]' of obs for person i ;
  * eta_i is the K-vector [\eta_{i1},...,\eta_{iK}]' for death i ;
  * Other components:
  * l and k specify which entries of \xi you want to sample! Note this is C++ zero-indexing.
  */
  
  if( fix_xi==1 ){
    return (beta_c.t() * X_i).eval()(0);
  } else{ 
    int P = Theta_c.n_rows;
    int L = Theta_c.n_cols;
    int K = eta_i.n_rows;
    
    // Loop through to calculate m_sv and c_sv, shorthand holders for products
    arma::vec m_sv(P); m_sv.fill(0.0);
    arma::vec c_sv(P); c_sv.fill(0.0);
    for( int pp=0; pp<P; pp++ ) {
      for( int ll=0; ll<L; ll++ ) {
        if( ll!=l_samp ){ // Only not-of-interest ll
          for( int kk=0; kk<K; kk++ ) { // ALL kk
            m_sv(pp) += Theta_c(pp,ll)*xi_i(ll,kk)*eta_i(kk); 
          } // for( int kk=0; kk<K; kk++ )
        } else{ // if ll==l_samp, i.e. the of-interest ll
          for( int kk=0; kk<K; kk++ ) { 
            if( kk!=k_samp ) { // Only not-of-interest kk
              m_sv(pp) += Theta_c(pp,ll)*xi_i(ll,kk)*eta_i(kk); 
            }
          } // for( int kk=0; kk<K; kk++ )
        }
      } // for( int ll=0; ll<L; ll++ )
      c_sv(pp) = Theta_c(pp,l_samp) * eta_i(k_samp);
    } // for( int pp=0; pp<P; pp++ ) 
    
    // Loop through to calculate M_star and sigsq_star_inv, for use in final sampling
    double M_star;
    double sigsq_star_inv; 
    M_star = (beta_c.t() * X_i / sig_sq_xi).eval()(0);
    sigsq_star_inv = 1/sig_sq_xi;
    for( int pp=0; pp<P; pp++ ) {
      M_star += c_sv(pp) * (z_i(pp) - m_sv(pp)) / Sigma_0_vec(pp) ;
      sigsq_star_inv += c_sv(pp) * c_sv(pp) / Sigma_0_vec(pp) ; 
    } // for( int pp=0; pp<P; pp++ ) 
    
    // Sample from the final (univariate) normal!
      return rnormArma( M_star/sigsq_star_inv, 1/sigsq_star_inv ); 
  }
}

// [[Rcpp::export]]
Rcpp::List sample_xi_all(int num_causes, int L, int K, Rcpp::List N,
                         arma::mat sigSqXi_all, Rcpp::List beta_c_all,
                         Rcpp::List X_all_sig, arma::vec Sigma_0_vec,
                         Rcpp::List Theta_all, Rcpp::List xi_all,
                         Rcpp::List z_all, Rcpp::List eta_all, int fix_xi) {
  /* FIXME describe
   */ 
  
  for(int c=0; c<num_causes; c++){
    int N_c = N[c];
    arma::cube xi_c = xi_all[c];
    arma::cube beta_c = beta_c_all[c];
    arma::mat X_sig_c = X_all_sig[c];
    arma::mat Theta_c = Theta_all[c];
    arma::mat z_c = z_all[c];
    arma::mat z_c_t = z_c.t();
    arma::mat eta_c = eta_all[c];
    arma::mat eta_c_t = eta_c.t();
    for(int i=0; i<N_c; i++){
      arma::mat xi_c_i = xi_c.slice(i);
      arma::vec X_sig_c_i = X_sig_c.row(i).t();
      arma::vec z_c_t_i = z_c_t.row(i).t();
      arma::vec eta_c_t_i = eta_c_t.row(i).t();
      for(int l=0; l<L; l++){
        for(int k=0; k<K; k++){
          arma::vec beta_c_lk = beta_c.subcube( 0,l,k, beta_c.n_rows-1,l,k );
          double sigSqXi_lk = sigSqXi_all(l,k);
          xi_c_i(l,k) = sample_xi_ilk(sigSqXi_lk, beta_c_lk, X_sig_c_i, Sigma_0_vec, 
                 Theta_c, xi_c_i, z_c_t_i, eta_c_t_i, l, k, fix_xi);
        }
      }
      xi_c.slice(i) = xi_c_i;
    }
    xi_all[c] = xi_c;
  }
  return(xi_all);
}

// [[Rcpp::export]]
double sample_psi_ik(double sig_sq_psi, arma::vec alpha_c, arma::vec X_i,
                     arma::vec Sigma_0_vec, arma::mat Omega_i, arma::vec psi_i, 
                     arma::vec z_i, arma::vec nu_i, int k_samp) {
  /* 
  * Sample \psi_{c,k}(x_i) where the result is a double for person i having COD c.
  * 
  * The inputs are comprised of "prior" and "likelihood" components;
  * "Prior" components:
  * sig_sq_psi, the variance term associated with obs of \psi_i ; 
  * X_i, the B x 1 covariate vector (X_{i,*})^T for death i with COD c ; 
  * alpha_c is the B x 1 cause-specific B vector of coefficients for the regression of \psi ;
  * "Likelihood" components:
  * Sigma_0_vec is the P-vector diag(\Sigma_0) ;
  * Omega_i is the P x K matrix \Omega_i = \Theta_{c[i]} \xi_i specific to person i ; 
  * psi_i is the current K x 1 vector \psi_i person i dying of cod c ; 
  * z_i, the P*1 vector [z_{i1},...,z_{iP}]' of obs for person i ;
  * nu_i is the K-vector [\nu_{i1},...,\nu_{iK}]' for death i ;
  * Other components:
  * k specifies which entry of \psi you want to sample! Note this is C++ zero-indexing.
  */
  
  int P = Omega_i.n_rows;
  int K = Omega_i.n_cols;
  
  // Loop through to calculate m_sv and c_sv, shorthand holders for products
  arma::vec m_sv(P); m_sv.fill(0.0);
  arma::vec Om_nu = Omega_i * nu_i ; // P x 1 vector
  for( int pp=0; pp<P; pp++ ) {
    for( int kk=0; kk<K; kk++ ) { 
      if( kk!=k_samp ) { // Only not-of-interest kk
        m_sv(pp) += Omega_i(pp,kk) * psi_i(kk) ; 
      }
    } // for( int kk=0; kk<K; kk++ )
    m_sv(pp) += Om_nu(pp) ;
  } // for( int pp=0; pp<P; pp++ ) 
  
  // Loop through to calculate M_star and sigsq_star_inv, for use in final sampling
  double M_star;
  double sigsq_star_inv; 
  M_star = (alpha_c.t() * X_i / sig_sq_psi).eval()(0);
  sigsq_star_inv = 1/sig_sq_psi;
  for( int pp=0; pp<P; pp++ ) {
    M_star += Omega_i(pp,k_samp) * (z_i(pp) - m_sv(pp)) / Sigma_0_vec(pp) ;
    sigsq_star_inv += Omega_i(pp,k_samp) * Omega_i(pp,k_samp) / Sigma_0_vec(pp) ; 
  } // for( int pp=0; pp<P; pp++ ) 
  
  // Sample from the final (univariate) normal!
  return rnormArma( M_star/sigsq_star_inv, 1/sigsq_star_inv ); 
}

// [[Rcpp::export]]
Rcpp::List sample_psi_all(Rcpp::List N, int K, int num_causes, arma::vec sig_sq_psi,
                          Rcpp::List alpha_c_all, Rcpp::List X_all_mu, arma::vec Sigma_0_vec, 
                          Rcpp::List Omega_all, Rcpp::List psi_all, Rcpp::List z_all, 
                          Rcpp::List nu_all){
  /* FIXME describe
   */ 
  
  for(int c=0; c<num_causes; c++){
    int N_c = N[c];
    arma::mat alpha_c = alpha_c_all[c];
    arma::mat X_mu_c = X_all_mu[c];
    arma::cube Omega_c = Omega_all[c];
    arma::mat psi_c = psi_all[c];
    arma::mat psi_c_t = psi_c.t();
    arma::mat z_c = z_all[c];
    arma::mat z_c_t = z_c.t();
    arma::mat nu_c = nu_all[c];
    arma::mat nu_c_t = nu_c.t();
    for(int i=0; i<N_c; i++){
      arma::vec X_mu_c_i = X_mu_c.row(i).t();
      arma::mat Omega_c_i = Omega_c.slice(i);
      arma::vec psi_c_i = psi_c_t.row(i).t();
      arma::vec z_c_i = z_c_t.row(i).t();
      arma::vec nu_c_i = nu_c_t.row(i).t();
      for(int k=0; k<K; k++){
        double sigSqpsi_k = sig_sq_psi(k);
        arma::vec alpha_c_k = alpha_c.col(k);
        psi_c_i(k) = sample_psi_ik(sigSqpsi_k, alpha_c_k, X_mu_c_i, Sigma_0_vec, 
                                   Omega_c_i, psi_c_i, z_c_i, nu_c_i, k);
        
      }
      psi_c_t.row(i) = psi_c_i;
    }
    psi_all[c] = psi_c_t.t();
  }
  
  return psi_all;
}

// [[Rcpp::export]]
arma::mat sample_Theta_j(double sig_sq_j, arma::cube xi_c, arma::mat eta_c, arma::vec z_cj, 
                        arma::vec Delta_j, arma::vec phi_j, arma::vec tau) {
  /* 
  * Sample \Theta_{c,j*} ; the result is a length L vector where entry l is \Theta_{c,jl}.
  * 
  * The inputs are sig_sq_j, \sigma_j^2, the "main" error variance ;
  * xi_c is the L x K x N_c cube where xi_c[,,i] is the LxK matrix \xi(x_i) ;
  * eta_c is the K x N_c matrix [\eta_{c,1},...,\eta_{c,N_c}], \eta_{c,i} length-K vec for ith obs; 
  * z_cj is length N_c vector {z_{c, 1 j},...,z_{c, N_c j}}', j is the row index of interest ;
  * Delta_j is the length-L vector {\Delta_{j1},...,\Delta_{jL}}' ; 
  * phi_j is the length-L vector {\phi{j1},...,\phi{jL}}' ; 
  * tau is the length-L vector {\tau_1,...,\tau_L}' ; 
  * 
  */
  int L = xi_c.n_rows;
  int N_c = xi_c.n_slices;
  // Set up necessary matrices
  // tilde_eta_c is L x N_c, \tilde{\eta_c} = {\xi(x_1)\eta_1, ..., \xi(x_{N_c})\eta_{N_c}}'
  arma::mat tilde_eta_c(N_c, L);
  for( int i=0; i<N_c; i++ ) {
    tilde_eta_c.row(i) = (xi_c.slice(i) * eta_c.col(i)).t();
  }
  arma::mat tilde_eta_c_t = tilde_eta_c.t();
  // Get overall mean and cov components
  arma::mat Sigma_star_inv = tilde_eta_c_t*tilde_eta_c/sig_sq_j + diagmat( phi_j%tau ) ;
  arma::vec mu_tmp = tilde_eta_c_t*z_cj/sig_sq_j + Delta_j%phi_j%tau ; // % for element-wise mult
  arma::vec mu_star = arma::solve(Sigma_star_inv,mu_tmp);
  // Sigma_star_inv is LxL, MV normal sample!
  return arma::mvnrnd(mu_star, arma::inv(Sigma_star_inv), 1);
}

// [[Rcpp::export]]
double sample_Delta_jl(double phi_jl_delta, double tau_l_delta,
                       arma::vec Theta_jl, double phi_jl_theta, double tau_l_theta) {
  /* 
  * Sample \Delta_{jl} ; the result is a double \Delta_{jl}.
  * 
  * phi_jl_delta is entry \phi_{\delta,jl} ; 
  * tau_l_delta is the entry \tau_{\delta,l} ; 
  * Theta_jl is C vector {\Theta_{1,jl},...,\Theta_{C,jl}}' for jl entries of Theta across causes; 
  * phi_jl_theta is the entry \phi_{\theta,jl} ; 
  * tau_l_theta is the entry \tau_{\theta,l} ; 
  * 
  */
  int C = Theta_jl.n_rows;
  // Get overall mean and cov components
  double Sigma_star_inv =  phi_jl_delta*tau_l_delta + C*phi_jl_theta*tau_l_theta;
  double mu_tmp = phi_jl_theta*tau_l_theta*arma::sum(Theta_jl) ;
  double mu_star = mu_tmp/Sigma_star_inv ;
  // Sigma_star_inv is LxL, MV normal sample!
  return rnormArma(mu_star, 1/Sigma_star_inv);
}

// [[Rcpp::export]]
arma::mat sample_nu(arma::mat xi_i, arma::mat Theta_c, arma::mat Sigma_0, 
                    arma::mat z_i, arma::vec psi_i){
  /* 
   * Sample \nu_{i} ; the result is a length-k vector.
   * 
   * xi_i is L x K matrix \xi(x_i) for observation i ;
   * Theta_c is the P x L \Theta matrix for COD c ;
   * Sigma_0 is the P x P covariance matrix for the inputs ; 
   * z_i is the length-P vector of latent observations ;
   * psi_i is K vector \psi(x_i) .
   * 
  */
  int K = psi_i.n_elem;
  // Get necessary sub-terms
  // Omega_i is P x K \Omega_i = \Theta \xi_{c[i]}(x_i) ;
  arma::mat Omega_i = Theta_c * xi_i;
  arma::vec tilde_z = z_i - Omega_i * psi_i;
  arma::mat I_k(K,K,arma::fill::eye);
  arma::mat tmp_mult = Omega_i.t() * arma::inv(Sigma_0);
  // Get overall mean and cov components
  arma::mat Sigma_star_inv = I_k + tmp_mult * Omega_i ;
  arma::vec mu_tmp = tmp_mult * tilde_z;
  arma::vec mu_star = arma::solve(Sigma_star_inv,mu_tmp);
  // Sigma_star_inv is kxk, so sample as a MV normal since it's cheap anyway!
  return arma::mvnrnd(mu_star, arma::inv(Sigma_star_inv), 1);
}

/* *********** SHRINKAGE PRIOR SAMPLER FUNCTIONS *********** 
 * 
 * Used for sampling hyper-parameters associated with
 * \theta_{c,lk} and \delta_{lk}.
 * 
 */

// [[Rcpp::export]]
double sample_phi_pl(double nu, double tau_deltaSq, int C) {
  /* 
   * General function for sampling \phi_{pl}, given \phi_{pl} ~ Ga(nu/2,nu/2)
   * where taul_deltaplSq is either \tau_l \Delta_{pl}^2  and C=1
   * (for sampling shrinkage term for the \Delta_{pl} components),
   * or \tau_l \sum_{c=1}^C (\Theta_{c,pl}-\Delta_{pl})^2  and C=num_causes
   * (for sampling common shrinkage term for the \Theta_{c,pl} components).
   * Here p \in {1,...,P} and l \in {1,...,L}.
   * Note that arma::randg uses scale parameterization, i.e. 1/rate. 
   */
  arma::vec samp = arma::randg(1, arma::distr_param((nu+C)/2, 2/(nu+tau_deltaSq)));
  return samp(0);
}

// [[Rcpp::export]]
double sample_delta_theta(double a1, double a2, arma::mat Delta, arma::cube Theta,
                          arma::mat phi_theta, arma::vec delta_theta, int h) {
  /* 
  * Function for sampling \delta_{1}, given \delta_{1} ~ Ga(a1,1) a priori
  * and \delta_{h} ~ Ga(a2,1) h \in {2,...,L} a priori. 
  * Theta is P x L x C cube where slice c is P x L Theta matrix for COD c.
  * Delta and phi_theta are P x L, delta_theta is L x 1 (this is old delta sample).
  * h \in {0,...,L-1} specifies which element of delta is being sampled.
  * Note that arma::randg uses scale parameterization, i.e. 1/rate.
  * 
  */
  int P = Theta.n_rows;
  int L = Theta.n_cols;
  int C = Theta.n_slices;
  double shape;
  double rate;
  
  // Specify shape and initial rate parameters
  if(h==0){
    shape = a1 + C*P*L/2;
    rate = 1.0;
  } else{ 
    shape = a2 + C*P*(L-h)/2;
    rate = 1.0/delta_theta(h-1); // Evan's modification
  }
  
  // Calculate the components of the Gamma posterior
  for( int c=0; c<C; c++ ){
    double tau_minus = 1.0;
    for( int l=0; l<L; l++ ){
      if( l != h ){
        tau_minus = tau_minus * delta_theta(l);
      }
      if( l >= h ){
        double tmp_sum = 0.0;
        for( int j=0; j<P; j++ ){
          double tmp_diff = Theta(j,l,c) - Delta(j,l);
          tmp_sum += phi_theta(j,l) * tmp_diff * tmp_diff;
        }
        rate += 0.5*tau_minus*tmp_sum;
      }
    }
  }
  // Sample a new value for \delta_1
  double delta_samp = arma::randg( 1,arma::distr_param(shape,1/rate) )(0);
  return delta_samp;
}

// [[Rcpp::export]]
double sample_delta_Delta(double a1, double a2, arma::mat Delta, arma::mat phi_delta, 
                          arma::vec delta_delta, int h) {
  /* 
  * Function for sampling \delta_{1}, given \delta_{1} ~ Ga(a1,1) a priori
  * and \delta_{h} ~ Ga(a2,1) h \in {2,...,L} a priori. 
  * Delta and phi_delta are P x L, delta_delta is L x 1 (this is old delta sample).
  * Note that arma::randg uses scale parameterization, i.e. 1/rate.
  * 
  */
  int P = Delta.n_rows;
  int L = Delta.n_cols;
  double shape;
  double rate;
  
  // Specify shape and initial rate parameters
  if(h==0){
    shape = a1 + P*L/2;
    rate = 1.0;
  } else{ 
    shape = a2 + P*(L-h)/2; 
    rate = 1.0/delta_delta(h-1); // Evan's modification
  }
  
  // Calculate the components of the Gamma posterior
  double tau_minus = 1.0;
  for( int l=0; l<L; l++ ){
    if( l != h ){
      tau_minus = tau_minus * delta_delta(l);
    }
    if( l >= h ){
      double tmp_sum = 0.0;
      for( int j=0; j<P; j++ ){
        tmp_sum += phi_delta(j,l) * Delta(j,l) * Delta(j,l);
      }
      rate += 0.5*tau_minus*tmp_sum;
    }
  }
  //return rate; // uncomment and comment return delta_samp; below for debugging
  // Sample a new value for \delta_1
  double delta_samp = arma::randg( 1,arma::distr_param(shape,1/rate) )(0);
  return delta_samp;
}

/* *********** HIERARCHICAL REGRESSION SAMPLER FUNCTIONS *********** 
 * 
 * Used for sampling \alpha_{c,k}, \beta_{c,lk}, 
 * and the population-level means and variances of each.
 * 
 */

// [[Rcpp::export]]
arma::mat sample_beta_c(arma::vec y, arma::vec mu_beta, arma::mat Sigma_beta, 
                        double sig_sq, arma::mat XtX, arma::mat Xt) {
  /* Here y can be either \xi_{c,lk} or \psi_{c,k}, with their respective
  * regression hyperparameters filled in depending on which is input;
  * the output will be \beta_{c,\xi_{lk}} or \alpha_{c,\psi_{k}}, resp,
  * the result of y = X\beta + eps regression with given prior on \beta.
  * The input mu_beta is the B-vector "population" level mean for the beta_c's,
  * and Sigma_beta is the B x B matrix for the "population" level covariance.
  * 
  * \beta_{c,\xi_{lk}} ~ N( \tilde{\Sigma} \tilde{\Mu}, \tilde{\Sigma})
  * \tilde{\Mu} = \Sigma_{\beta_{lk}}^{-1} \mu_{\beta_{lk}} + X^T \xi_{c,lk} / \sigma_{\xi_{lk}}^2
  * \tilde{\Sigma}^{-1} = (\Sigma_{\beta_{lk}}^{-1} + X^T X / \sigma_{\xi_{lk}}^2)
  */
  arma::mat Sigma_tilde_inv = arma::inv(Sigma_beta) + XtX/sig_sq;
  arma::vec mu_tilde = arma::solve(Sigma_beta,mu_beta) + Xt * y / sig_sq;
  return arma::mvnrnd(arma::solve(Sigma_tilde_inv,mu_tilde), arma::inv(Sigma_tilde_inv), 1);
}

// [[Rcpp::export]]
arma::mat sample_beta_mu(arma::vec mu_0, arma::mat Lambda_0, 
                         double C, arma::vec beta_mean, arma::mat Sigma_beta) {
  /* Here we learn the population-level mean of the \beta_{c,lk} draws,
  * i.e. \beta_{lk} from \beta_{c,lk} ~ N(\beta_{lk}, \Sigma_{\beta_{lk}}),
  * given we observe \beta_{c,lk} for c=1,...,C,
  * where the prior \beta_{lk} ~ N(\mu_0, \Lambda_0).
  * 
  * \beta_{lk} ~ N( \tilde{\Sigma} \tilde{\Mu}, \tilde{\Sigma})
  * \tilde{\Mu} = \Lambda_0^{-1} \mu_0 + C \Sigma_{\beta_{lk}}^{-1} \bar{\beta_{c,lk}}
  * \tilde{\Sigma}^{-1} = \Lambda_0^{-1} + C \Sigma_{\beta_{lk}}^{-1}
  */
  arma::mat Lambda_tilde_inv = arma::inv(Lambda_0) + C*arma::inv(Sigma_beta);
  arma::vec mu_tilde = arma::solve(Lambda_0,mu_0) + C*arma::solve(Sigma_beta,beta_mean);
  return arma::mvnrnd(arma::solve(Lambda_tilde_inv,mu_tilde), arma::inv(Lambda_tilde_inv), 1);
}

// [[Rcpp::export]]
arma::mat sample_beta_Sigma(int v0, arma::mat S0, int C,
                            arma::vec beta_mu, arma::mat beta_c) {
  /* Here we learn the population-level variance of the \beta_{c,lk} draws,
   * i.e. \Sigma_{\beta_{lk}} from \beta_{c,lk} ~ N(\beta_{lk}, \Sigma_{\beta_{lk}}),
   * given we observe \beta_{c,lk} for c=1,...,C,
   * where the prior \Sigma_{\beta_{lk}}) ~ inv. Wishart(v_0, S_0).
   * 
   * Enter in \beta_c as [\beta_{1,lk}, \beta_{2,lk}, ..., \beta_{C,lk}],
   * that is make each group-level \beta a column of \beta_c (\beta_{c,lk} is B-vector).
   * 
   * The posterior is \Sigma_{\beta_{lk}}) ~ inv. Wishart(v_0+C, (S_0+S_C)^{-1})
   * where S_C = \sum_{c=1}^C (\beta_{c,lk}-\beta_{lk})(\beta_{c,lk}-\beta_{lk})^T
   */
  arma::mat SC(beta_mu.n_rows,beta_mu.n_rows);
  SC.fill(0.0);
  for( int c=0; c<C; c++ ) {
    SC = SC + (beta_c.col(c)-beta_mu) * ((beta_c.col(c)-beta_mu).t());
  }
  // https://en.wikipedia.org/wiki/Inverse-Wishart_distribution#Conjugate_distribution
  return arma::iwishrnd( S0 + SC, v0 + C);
}

// [[Rcpp::export]]
double sample_sigsq(double a, double b, int n, double RSS) {
  /* 
   * General function for sampling \sigma^2, given 1/\sigma^2 ~ Ga(a/2,b/2)
   * and you observe n samples with residual sum of squares RSS.
   * Note that arma::randg uses scale parameterization, i.e. 1/rate. 
   */
  arma::vec samp = arma::randg(1, arma::distr_param((a+n)/2,2/(b+RSS)));
  return 1/samp(0);
}


/* *********** NEW/TEST DATA SAMPLER FUNCTIONS *********** 
* 
* Used for monte carlo samples of \xi_{c,i*} and \eta_{c,i*} 
* for decedent i* having unknown COD (so done for each c in main sampler).
* 
*/

// [[Rcpp::export]]
arma::mat sample_xi_istar(arma::mat sig_sq_xi, arma::cube beta_c, arma::vec X_istar, int fix) {
  /* 
  * Sample \xi_{c,lk}(x_i*) where the result is a double for person i* having COD c.
  * 
  * The sampler are drawn from the "prior" on \xi;
  * sig_sq_xi is L x K, the variance terms associated with each entry of \xi ; 
  * X_i, the B x 1 covariate vector (X_{i*,:})^T for death i* ; 
  * beta_c is the B x L x K cube of cause-specific vectors of coefficients for the regression of \xi_{c,lk} ;
  */
  int B = beta_c.n_rows;
  int L = beta_c.n_cols;
  int K = beta_c.n_slices;
  arma::mat xi_istar(L,K);
  arma::vec beta_clk;
  for( int ll=0; ll<L; ll++ ){
    for( int kk=0; kk<K; kk++ ){
      beta_clk = beta_c.subcube( 0, ll, kk, B-1, ll, kk);
      if( fix==0 ){
        xi_istar(ll,kk) = rnormArma( (beta_clk.t() * X_istar).eval()(0), sig_sq_xi(ll,kk) );
      } else{ xi_istar(ll,kk) = (beta_clk.t() * X_istar).eval()(0); }
    }
  }
  // Return the L x K matrix xi_istar!
  return xi_istar;  
}

// [[Rcpp::export]]
arma::mat sample_psi_istar(arma::vec sig_sq_psi, arma::mat alpha_c, arma::vec X_istar, int fix_psi) {
  /* 
  * Sample \xi_{c,lk}(x_i*) where the result is a double for person i* having COD c.
  * 
  * The sampler are drawn from the "prior" on \xi;
  * sig_sq_xi is L x K, the variance terms associated with each entry of \xi ; 
  * X_i, the B x 1 covariate vector (X_{i*,:})^T for death i* ; 
  * beta_c is the B x L x K cube of cause-specific vectors of coefficients for the regression of \xi_{c,lk} ;
  */
  int B = alpha_c.n_rows;
  int K = alpha_c.n_cols;
  arma::vec psi_istar(K);
  arma::vec alpha_ck;
  for( int kk=0; kk<K; kk++ ){
    alpha_ck = alpha_c.col(kk);
    if( fix_psi== 0 ){
      psi_istar(kk) = rnormArma( (alpha_ck.t() * X_istar).eval()(0), sig_sq_psi(kk) );
    } else{
      psi_istar(kk) = (alpha_ck.t() * X_istar).eval()(0);
    }
  }
  // Return the L x K matrix xi_istar!
  return psi_istar;  
}

/* *********** EXPANDED SAMPLER FUNCTIONS (SPEEDING UP SLOW R LOOPS) *********** 
 * 
 * FIXME DESCRIBE!
 * 
 */

// [[Rcpp::export]]
Rcpp::List sample_z_mean_cov_all(Rcpp::List S_mat, Rcpp::List Omega_all, Rcpp::List eta_all, arma::vec N, 
                                 arma::vec Sig0vec, Rcpp::List z_all_nominus, Rcpp::List z_all, int P,
                                 arma::vec is_binary, Rcpp::List mu_all, bool mu_collapse, Rcpp::List psi_all){
  /* 
   * FIXME DESCRIBE!
   */
  int num_causes = S_mat.size();
  double inf = std::numeric_limits<double>::infinity();
  // Initialize lists to hold results
  Rcpp::List mean_all_list(num_causes);
  Rcpp::List cov_all_list(num_causes);
  Rcpp::List z_all_list(num_causes);
  Rcpp::List z_all_nominus_list(num_causes);
  for(int c=0; c<num_causes; c++){
    // Initialize empty mean_all and cov_all arrays
    arma::mat mean_all_t(N(c),P);
    arma::cube cov_all(P, P, N(c)); // FIXME better order?
    // Save list element as matrix/vector for quick access later
    arma::mat S_mat_c = S_mat[c];
    arma::mat eta_all_c = eta_all[c];
    arma::cube Omega_all_c = Omega_all[c]; // FIXME better order?
    arma::mat psi_all_c = psi_all[c];
    arma::mat psi_all_c_t = psi_all_c.t();
    // Update/sample things (provided to sampler, but also filled in)
    arma::mat z_all_nominus_c = z_all_nominus[c];
    arma::mat z_all_nominus_c_t = z_all_nominus_c.t();
    arma::mat z_all_c = z_all[c];
    arma::mat z_all_c_t = z_all_c.t();
    arma::mat mu_all_c = mu_all[c];
    arma::mat mu_all_c_t = mu_all_c.t();
    for(int i=0; i<N(c); i++){
      arma::mat Si_vec = S_mat_c.row(i);
      arma::mat Omega_all_c_i = Omega_all_c.slice(i);
      for(int j=0; j<P; j++){
        if( NumericVector::is_na( Si_vec(j) ) ){ // Impute z_{i,j} for missing s_{i,j}.
          double mu_tmp = (Omega_all_c_i.row(j) * eta_all_c.col(i)).eval()(0);
          z_all_nominus_c_t(i,j) = rnormArma(mu_tmp, Sig0vec(j));
        } else{ // Sample z_{i,j} for binary j.
          if( is_binary(j) ){
            // Omega_all[[c]][j,,i] is 1 x K row vector
            // eta_all[[c]][,i] is K x 1 col vector
            double mu_tmp = (Omega_all_c_i.row(j) * eta_all_c.col(i)).eval()(0);
            double min_norm, max_norm;
            if( Si_vec(j) == 1 ){
              min_norm = 0; 
              max_norm = inf;
            } else {
              min_norm = -inf; 
              max_norm = 0;
            }
            z_all_nominus_c_t(i,j) = rtruncnorm(mu_tmp, 1, min_norm, max_norm);
          }
          // Subtract off the mean from z_all_nominus to get z_all (no matter where z_all_nominus came from).
          if( mu_collapse ){
            z_all_c_t(i,j) = z_all_nominus_c_t(i,j);
          } else{
            z_all_c_t(i,j) = z_all_nominus_c_t(i,j) - mu_all_c_t(i,j); // Equal to z_all_nominus when mu is 0.
          }
        }
      } // FIXME mu_all_c and psi_all_c slow?
      mean_all_t.row(i) = mu_all_c_t.row(i) + (Omega_all_c_i * (psi_all_c_t.row(i)).t()).t(); // Equal to \lambda\eta when mu is 0
      cov_all.slice(i) = Omega_all_c_i * (Omega_all_c_i).t();
      for(int j=0; j<P; j++){
        cov_all.slice(i)(j,j) += Sig0vec(j);
      }
    } // for(int i=0; i<N(c); i++)
    // Save mean, cov, and z in list to return.
    mean_all_list(c) = mean_all_t.t();
    cov_all_list(c) = cov_all;
    z_all_list(c) = z_all_c_t.t();
    z_all_nominus_list(c) = z_all_nominus_c_t.t();
  } // for(int c=0; c<num_causes; c++)
  Rcpp::List L = List::create(Named("mean_all") = mean_all_list, 
                              _["cov_all"] = cov_all_list, 
                              _["z_all"] = z_all_list, 
                              _["z_all_nominus"] = z_all_nominus_list);
  return(L);
}

// [[Rcpp::export]]
arma::mat get_piSgivenY(int N_test, int num_causes, int P, int mc_tot, bool cov_incl, int fix_xi, 
                        arma::mat X_test_mu, arma::mat X_test_sig, arma::mat S_test,
                        arma::mat sigSqXi_all, Rcpp::List beta_c_all,
                        Rcpp::List Theta_all, bool mu_collapse, Rcpp::List gamma_c_all,
                        Rcpp::List alpha_c_all, 
                        arma::mat sigSqpsi_all, arma::mat Sigma_0, 
                        arma::vec is_binary){
  /* FIXME describe
  */ 
  arma::mat pi_SgivenY(N_test, num_causes);
  pi_SgivenY.fill(0.0);
  arma::mat xi_star;
  arma::mat Omega_star;
  arma::mat mu_star(P,1); 
  arma::mat psi_star;
  arma::mat mean_star(P,1);
  arma::mat cov_star(P,P);
  for( int c=0; c<num_causes; c++ ){
    arma::mat alpha_c = alpha_c_all[c];
    arma::mat Theta_c = Theta_all[c];
    arma::cube beta_c = beta_c_all[c];
    arma::mat gamma_c = gamma_c_all[c];
    if( !cov_incl ){ // no covariate info here, so only need to do once per cause
      // "Sample" xi (but since it's fixed, just get beta X for person i)
      xi_star = sample_xi_istar(sigSqXi_all, beta_c, X_test_sig.row(1).t(), fix_xi);
      // Calculate Omega_star
      Omega_star = get_Omega_i(Theta_c, xi_star);
      // Get mean and covariance for person i assuming their COD is c
      mu_star.fill(0);
      if( !mu_collapse ){ // If we have z_i = \mu_{c[i]} + \Lambda\eta_i, \eta_i~N(0, I), set mu_all[[c]].
        for( int p=0; p<P; p++ ){ mu_star(p) = (X_test_mu.row(1) * gamma_c.col(p)).eval()(0); }
      } // if( !mu_collapse )
      psi_star = sample_psi_istar(sigSqpsi_all, alpha_c, X_test_mu.row(1).t(), 1); // Fixed bc 1
      mean_star = mu_star + Omega_star * psi_star; // NOTE using eta = psi here, i.e. eta~N(alpha x, 1).
      cov_star = Omega_star * Omega_star.t() + Sigma_0;
    } // if( !cov_incl )
    for( int i=0; i<N_test; i++ ){
      if( cov_incl ){
        // "Sample" xi (but since it's fixed, just get beta X for person i)
        xi_star = sample_xi_istar(sigSqXi_all, beta_c, X_test_sig.row(i).t(), fix_xi);
        // Calculate Omega_star
        Omega_star = get_Omega_i(Theta_c, xi_star);
      } // if( cov_incl )
      // Calculate via Monte Carlo Approximation
      // Define variables that don't vary at each iter of MC approximation
      arma::mat S_test_i = S_test.row(i);
      arma::mat X_test_i = X_test_mu.row(i);
      arma::mat X_test_i_t = X_test_i.t();
      std::vector<bool> obs_ind(P);
      for( int p=0; p<P; p++ ){
        obs_ind[p] = !(NumericVector::is_na( S_test_i(p) ));
      }
      cov_star = Omega_star * Omega_star.t() + Sigma_0;
      mu_star.fill(0);
      if( !mu_collapse ){ // If we have z_i = \mu_{c[i]} + \Lambda\eta_i, \eta_i~N(0, I), set mu_all[[c]].
        for( int p=0; p<P; p++ ){ mu_star(p) = (X_test_i * gamma_c.col(p)).eval()(0); }
      }
      // Get monte carlo approximation of \pi(s_i | y_i)
      arma::vec logp_tmp(P);
      for( int mc=0; mc<mc_tot; mc++ ){
        // Sample eta_star, a K-vector of \eta for person i* with covs X_test[i*,];
        // alpha_c_all[[c]] is B x K matrix, sigSqpsi_all is K-vector
        psi_star = sample_psi_istar(sigSqpsi_all, alpha_c, X_test_i_t, 0);
        // Get mean and covariance for person i assuming their COD is c
        mean_star = mu_star + Omega_star * psi_star; // NOTE using eta = psi here, i.e. eta~N(alpha x, 1).
        // Calculate the multivariate normal probability of observing the data given the cause c.
        // For s_{i*} binary, this is in the proper direction (positive when s is pos, neg when neg). 
        // For s_{i*} continuous, you can directly compute each pr(s_{ij}).
        logp_tmp.fill(0.0);
        for(int j=0; j<P; j++ ){
          if( obs_ind[j] ){
            double mu_tmp = mean_star(j);
            if( is_binary(j) ){ // ...then Sigma_0_vec[j] is fixed to 1
              if( S_test_i(j) == 1 ){ 
                logp_tmp(j) = log(1 - arma::normcdf(0.0, mu_tmp, 1.0)); // log probability that z_ij>0
              } else{ logp_tmp(j) = log(arma::normcdf(0.0, mu_tmp, 1.0)); } // log probability that z_ij<0 
            } else{ // ...then Sigma_0_vec[j] is free, and we actually OBSERVE z_{ij}
              logp_tmp(j) = log(arma::normpdf(S_test_i(j), mu_tmp, 1.0)); // std deviation parameterization
            } // if( is_binary(j) )
          } // if( obs_ind(j) )
        } // for(int j=0; j<P; j++ )
        // Pr(x1,...,xn) = Pr(x1)...Pr(xn) for x independent, so log(Pr(x1,...,xn)) = log(Pr(x1)...Pr(xn))
        // = log(Pr(x1)) + ... + log(Pr(xn)). So below, I can do exp(sum(logp_tmp)) to get Pr(x1,...,xn).
        pi_SgivenY(i,c) = pi_SgivenY(i,c) + exp(arma::sum(logp_tmp));
      } // for( int mc=0; mc<mc_tot; mc++ )
    } // for( int i=0; i<N_test; i++ )
  } // for( int c=0; c<num_causes; c++ )
  return(pi_SgivenY);
}

/* *********** R CODE BLOCKS
 * 
 * You can include R code blocks in C++ files processed with sourceCpp 
 * (useful for testing and development). The R code will be automatically
 * run after the compilation.
 * 
 */


/*** R

# S <- 1000; samp <- rep(NA,S); a<-1; b<-4
# for(i in 1:S){
#   samp[i] <- sample_sigsq(a, b, 0, 0)
# }
# hist(samp)
# hist(rgamma(S, shape=(a/2), rate = (b/2)))

# # make sure sample_xi works
# n_c <- 10;
# p <- 3;
# B <- 2;
# z_c <- matrix(rnorm(n_c*p))
# X_c <- matrix(runif(n_c* B), nrow=n_c, ncol=B)
# beta_c <- matrix(c(1,3))
# eta_c <- matrix(rnorm(n_c))
# theta_col_c <- matrix(rnorm(p))
# Sigma_0_vec <- matrix(rnorm(p))
# mu_eps <- matrix(rnorm(n_c*p))
# sig_sq_xi <- 0.01;
# sample_xi(sig_sq_xi, z_c, X_c, beta_c, eta_c, theta_col_c, Sigma_0_vec, mu_eps)
# 


# # Make sure sample_delta works
# # Off by 1 from how Isaac has it bc a1/a2 show up in shape AND rate for him
# p <- 5; k <- 4; l <- 3
# a1 = 1;
# a2 = 2;
# Delta = matrix(rnorm(p*l),nrow=p,ncol=l)
# phi = matrix(runif(p*l),nrow=p,ncol=l)
# delta = matrix(runif(l))
# # Isaac's code (I'm not sampling, just saving rates)
# delta_isaac <- delta
# tauh = cumprod(delta_isaac)
# Lambdasq_phi_sum = apply(Delta^2 * phi, 2, sum)
# rate_saves <- matrix(NA,nrow=l,ncol=2)
# rate = a1 + (1/2) * (1/delta_isaac[1]) * sum(tauh*Lambdasq_phi_sum)
# rate_saves[1,] <- c(rate, sample_delta_Delta(a1, a2, Delta, phi, delta, 0))
# #delta_isaac[1] = rgamma(1, shape = a1 + p*l/2, rate = rate)
# tauh = cumprod(delta_isaac)
# for(h in 2:l){
#   # This is the version from the original paper
#   rate=a2 + (1/2)*(1/delta_isaac[h]) * sum(tauh[h:l] * Lambdasq_phi_sum[h:l])
#   rate_saves[h,] <- c(rate, sample_delta_Delta(a1, a2, Delta, phi, delta, h-1))
#   #delta_isaac[h] = rgamma(1, a2 + p*(l - h + 1)/2, rate=rate)
#   tauh = cumprod(delta_isaac)
#   # This is modified based on the Durante note
#   #rate=a2 * (1/delta_isaac[h-1]) + (1/2)*(1/delta_isaac[h]) * sum(tauh[h:l] * Lambdasq_phi_sum[h:l])
#   #delta_isaac[h] = rgamma(1, a2 + p*(l - h + 1)/2, rate=rate)
#   #tauh = cumprod(delta_isaac)
# }
# rate_saves


#sourceCpp("/Users/Kelly/indite_all/papers_presentations/FARVA/code/main.cpp")

*/

