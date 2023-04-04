run_sampler_init <- function(S_mat, L, K, P, N, num_causes, B_sig, B_mu, 
                             a1_del_del, a2_del_del, a1_del_the, a2_del_the,
                             nu_delta, nu_theta, a_sigsq, b_sigsq,
                             mu_collapse, is_binary){

  suppressMessages( library(MCMCpack) )
  
  # Initialize everything based on its prior!
  # Unless you're in debug mode, then initialize by just random sampling N(0,1)
  
  #####################################################################
  
  # Build up \Delta and \Theta_c,
  # including extra shrinkage to keep sampler well-behaved,
  # starting with priors on the variance terms of each.
  # Sample \delta_l for l in 1...L for Theta and Delta.
  extra_shrink = 100
  delta_delta <- delta_theta <- rep(NA,L)
  delta_delta[1] <- rgamma(n=1,shape=a1_del_del,rate=1)
  delta_theta[1] <- rgamma(n=1,shape=a1_del_the,rate=1)
  for(l in 2:L){
    delta_delta[l] <- rgamma(n=1,shape=a2_del_del,rate=1)
    delta_theta[l] <- rgamma(n=1,shape=a2_del_the,rate=1)
  }
  tau_delta <- cumprod(delta_delta)
  tau_theta <- cumprod(delta_theta)
  # Sample \phi_{pl} for p in 1...P and l in 1...L for Theta and Delta,
  # and at the same time, sample \Delta_{pl} for p in 1...P and l in 1...L
  phi_delta <- phi_theta <- Delta <- matrix(NA,nrow=P,ncol=L)
  for(p in 1:P){
    for(l in 1:L){
      phi_delta[p,l] <- rgamma(1,shape=nu_delta/2,rate=nu_delta/2)
      phi_theta[p,l] <- rgamma(1,shape=nu_theta/2,rate=nu_theta/2)
      Delta[p,l] <- rnorm(1, mean=0, sd=(phi_delta[p,l]*tau_delta[l])^(-0.5)/extra_shrink)
    }
  }
  # Sample \Theta_{c,pl} for p in 1...P and l in 1...L and causes c in 1...num_causes
  Theta_all <- list()
  for(c in 1:num_causes){
    Theta_all[[c]] <- matrix(NA,nrow=P,ncol=L)
    for(p in 1:P){
      for(l in 1:L){
        Theta_all[[c]][p,l] <- rnorm(1, mean=Delta[p,l], sd=(phi_theta[p,l]*tau_theta[l])^(-0.5)/extra_shrink)
      }
    }
  }
  
  #####################################################################
  
  # Build up \xi, starting with priors on the hierarchical regression parameters for each.
  # First, the population-level coefficient mean and covariance for each entry of \Xi.
  # For each l, k, \beta_{lk} ~ N(\mu_{\beta,0}, \Lambda_{\beta,0})
  # and \Sigma_{\beta_{lk}} ~ inv. Wishart(v0_beta, S0_beta)
  beta_all <- array(0,c(B_sig,L,K))
  beta_Sigma <- array(NA,c(B_sig,B_sig,L,K))
  for(l in 1:L){
    for(k in 1:K){
      beta_Sigma[,,l,k] <- diag(B_sig)
    }
  }
  # Next, the cause-specific beta coefficients.
  # For each l, k, \beta_{c,lk} ~ N(\beta_{lk}, \Sigma_{\beta_{lk}})
  beta_c_all <- list()
  for(c in 1:num_causes){
    beta_c_all[[c]] <- array(0,c(B_sig,L,K))
  }
  # Then, fix \sigma_{\xi}^2 
  sigSqXi_all <- matrix(1,nrow=L,ncol=K)
  # Finally, sample each entry \xi_c ~ N(X_c * \beta_c, \sigma_{\xi}^2)
  xi_all <- list()
  for(c in 1:num_causes){
    xi_all[[c]] <- array(rnorm(n=prod(L,K,N[[c]]), mean=0, sd=0.01) ,c(L,K,N[[c]]))
  }
  
  Omega_all <- list()
  for(c in 1:num_causes){
    Omega_all[[c]] <- array(0,c(P,K,N[[c]])) # set up empty array to house obs
    #for(i in 1:(N[[c]]) ){ Omega_all[[c]][,,i] <- get_Omega_i(Theta_all[[c]], xi_all[[c]][,,i]) }
  }
  
  #####################################################################
  
  # Build up \psi, starting with priors on the hierarchical regression parameters for each.
  # First, the population-level coefficient mean and covariance for each entry of \psi.
  # For each k, \alpha_{k} ~ N(\mu_{\alpha,0}, \Lambda_{\alpha,0})
  # and \Sigma_{\alpha_{k}} ~ inv. Wishart(v0_alpha, S0_alpha)
  alpha_all <- matrix(0,nrow=B_mu,ncol=K)
  alpha_Sigma <- array(NA,c(B_mu,B_mu,K))
  for(k in 1:K){
    alpha_Sigma[,,k] <- diag(B_mu)
  }
  # Next, the cause-specific alpha coefficients.
  # For each k, \alpha_{c,k} ~ N(\alpha_{k}, \Sigma_{\alpha_{k}})
  alpha_c_all <- list()
  for(c in 1:num_causes){
    alpha_c_all[[c]] <- matrix(0,nrow=B_mu,ncol=K)
  }
  # Then, fix \sigma_{\psi}^2 to 1
  sigSqpsi_all <- rep(1,K)
  # Finally, sample each entry \psi_c ~ N(X_c * \alpha_c, \sigma_{\psi}^2)
  psi_all <- list()
  for(c in 1:num_causes){
    psi_all[[c]] <- matrix(rnorm(n=K*N[[c]], mean=0, sd=0.01),nrow=K,ncol=N[[c]])
  }
  
  #####################################################################
  
  # Get \eta_i for each i in 1...N[[c]] and and c in 1...num_causes.
  # Simply, \eta_i ~ N(\psi(x_i), I_k).
  eta_all <- nu_all <- list()
  for(c in 1:num_causes){
    nu_all[[c]] <- matrix(0,nrow=K,ncol=N[[c]])
    eta_all[[c]] <- psi_all[[c]] + nu_all[[c]]
  }
  
  # Finally, sample \sigma_p^2, i.e. elements of Sigma_0
  Sigma_0 <- matrix(0, nrow=P, ncol=P)
  for(p in 1:P){
    if( is_binary[p] ){
      Sigma_0[p,p] <- 1
    } else{
      Sigma_0[p,p] <- 1/rgamma(1,shape=a_sigsq/2,rate=b_sigsq/2)
    }
  }
  # Make Sigma_0_vec from Sigma_0
  Sigma_0_vec <- matrix(diag(Sigma_0))
  
  #####################################################################
  
  # Need to sample z_ij for s_ij binary!
  z_all_nominus <- list()
  for(c in 1:num_causes){
    z_all_nominus[[c]] <- matrix(NA,nrow=P,ncol=N[[c]])
    for(i in 1:N[[c]]){
      for(p in 1:P){
        if(is.na(S_mat[[c]][i,p])){
          z_all_nominus[[c]][p,i] <- rnorm(1, mean=0, sd=sqrt(0.01))
        } else{
          if(is_binary[p]){ # If binary, then we need to sample z!
            # Note I was having trouble with Inf and -Inf from rtruncnorm, so I modified.
            if(S_mat[[c]][i,p]==1){ # truncate z to live in (0,Inf)
              z_all_nominus[[c]][p,i] <- rtruncnorm(0, 0.01, 0, Inf)
            } else{ # truncate z to live in (-Inf,0]
              z_all_nominus[[c]][p,i] <- rtruncnorm(0, 0.01, -Inf, 0)
            }
          } else{ # If not binary, then z is just s!
            z_all_nominus[[c]][p,i] <- S_mat[[c]][i,p] # no need to sample if not binary!
          }
        }
      }
    }
  }
  
  
  #####################################################################
  # Initialize things associated with the mean mu, may just be fixed to 0 if mu_collapse==T
  
  mu_all <- z_all <- list()
  for(c in 1:num_causes){ 
    mu_all[[c]] <- matrix(0,nrow=P,ncol=N[[c]]) 
    z_all[[c]] <- z_all_nominus[[c]]
  }
  if( !(mu_collapse==T) ){ # If we have z_i = \mu_{c[i]} + \Lambda\eta_i, \eta_i~N(0, I), sample mu_all[[c]].
    # Sample \gamma_all_{p} (population-level mean) for each entry p\in{1,..,P};
    # the result is a length B_mu vector with entry b being the coefficient for covariate b
    gamma_all <- matrix(0,nrow=B_mu,ncol=P)
    gamma_Sigma <- array(NA,c(B_mu,B_mu,P))
    for(p in 1:P){
      gamma_Sigma[,,p] <- diag(B_mu)
    }
    
    # Next, the cause-specific gamma coefficients.
    gamma_c_all <- list()
    for(c in 1:num_causes){
      gamma_c_all[[c]] <- matrix(0,nrow=B_mu,ncol=P)
    }
    
    init_list=list('delta_delta'=delta_delta, 'delta_theta'=delta_theta, 
                   'tau_delta'=tau_delta, 'tau_theta'=tau_theta, 'phi_delta'=phi_delta, 
                   'phi_theta'=phi_theta, 'Delta'=Delta, 'Theta_all'=Theta_all, 
                   'beta_all'=beta_all, 'beta_Sigma'=beta_Sigma, 'beta_c_all'=beta_c_all, 
                   'sigSqXi_all'=sigSqXi_all, 'xi_all'=xi_all, 'Omega_all'=Omega_all, 
                   'alpha_all'=alpha_all, 'alpha_Sigma'=alpha_Sigma, 'alpha_c_all'=alpha_c_all, 
                   'sigSqpsi_all'=sigSqpsi_all, 'psi_all'=psi_all, 'eta_all'=eta_all, 
                   'nu_all'=nu_all, 'Sigma_0'=Sigma_0, 'Sigma_0_vec'=Sigma_0_vec, 
                   'z_all_nominus'=z_all_nominus, 'mu_all'=mu_all, 'z_all'=z_all,
                   'gamma_all'=gamma_all, 'gamma_Sigma'=gamma_Sigma, 'gamma_c_all'=gamma_c_all) 
    
  } else{
    init_list=list('delta_delta'=delta_delta, 'delta_theta'=delta_theta, 
                   'tau_delta'=tau_delta, 'tau_theta'=tau_theta, 'phi_delta'=phi_delta, 
                   'phi_theta'=phi_theta, 'Delta'=Delta, 'Theta_all'=Theta_all, 
                   'beta_all'=beta_all, 'beta_Sigma'=beta_Sigma, 'beta_c_all'=beta_c_all, 
                   'sigSqXi_all'=sigSqXi_all, 'xi_all'=xi_all, 'Omega_all'=Omega_all, 
                   'alpha_all'=alpha_all, 'alpha_Sigma'=alpha_Sigma, 'alpha_c_all'=alpha_c_all, 
                   'sigSqpsi_all'=sigSqpsi_all, 'psi_all'=psi_all, 'eta_all'=eta_all, 
                   'nu_all'=nu_all, 'Sigma_0'=Sigma_0, 'Sigma_0_vec'=Sigma_0_vec, 
                   'z_all_nominus'=z_all_nominus, 'mu_all'=mu_all, 'z_all'=z_all) 
  }
  
  return(init_list)

}
