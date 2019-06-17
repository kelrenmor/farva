sampleDist = function(n, num_causes, probs) { 
  sample(x = c(1:num_causes), n, replace = T, prob = probs) 
}

get_piSgivenY_analytic = function(sigSqXi_all, sigSqpsi_all, beta_c_all, alpha_c_all, 
                                  X_test_sig, X_test_mu, fix_xi, Theta_all,
                                  P, num_causes, N_test, mu_collapse, gamma_c_all, Sigma_0, 
                                  S_test, is_binary) {
  pi_SgivenY <- matrix(0, nrow=N_test, ncol=num_causes)
  for( c in 1:num_causes ){
    if( !cov_incl ){ # no covariate info here, so only need to do once per cause
      # "Sample" xi (but since it's fixed, just get beta X for person i)
      xi_star <- sample_xi_istar(sigSqXi_all, beta_c_all[[c]], matrix(X_test_sig[1,]), fix_xi)
      # Calculate Omega_star
      Omega_star <- get_Omega_i(Theta_all[[c]], xi_star)
      # Get mean and covariance for person i assuming their COD is c
      mu_star <- matrix(0,nrow=P,ncol=1)
      if( !(mu_collapse==T) ){ # If we have z_i = \mu_{c[i]} + \Lambda\eta_i, \eta_i~N(0, I), set mu_all[[c]].
        for(p in 1:P){mu_star[p] <- matrix(X_test_mu[1,], nrow=1) %*% matrix(gamma_c_all[[c]][,p])}
      }
      psi_star <- sample_psi_istar(matrix(sigSqpsi_all), alpha_c_all[[c]], matrix(X_test_mu[1,]), 1) # Fixed bc 1
      mean_star <- mu_star + Omega_star %*% psi_star # NOTE using eta = psi here, i.e. eta~N(alpha x, 1).
      cov_star <- Omega_star %*% t(Omega_star) + Sigma_0
    }
    for( i in 1:N_test ){
      if( cov_incl ){
        # "Sample" xi (but since it's fixed, just get beta X for person i)
        xi_star <- sample_xi_istar(sigSqXi_all, beta_c_all[[c]], matrix(X_test_sig[i,]), fix_xi)
        # Calculate Omega_star
        Omega_star <- get_Omega_i(Theta_all[[c]], xi_star)
        # Get mean and covariance for person i assuming their COD is c
        mu_star <- matrix(0,nrow=P,ncol=1)
        if( !(mu_collapse==T) ){ # If we have z_i = \mu_{c[i]} + \Lambda\eta_i, \eta_i~N(0, I), set mu_all[[c]].
          for(p in 1:P){mu_star[p] <- matrix(X_test_mu[i,], nrow=1) %*% matrix(gamma_c_all[[c]][,p])}
        }
        psi_star <- sample_psi_istar(matrix(sigSqpsi_all), alpha_c_all[[c]], matrix(X_test_mu[i,]), 1) # Fixed bc 1
        mean_star <- mu_star + Omega_star %*% psi_star # NOTE using eta = psi here, i.e. eta~N(alpha x, 1).
        cov_star <- Omega_star %*% t(Omega_star) + Sigma_0
        # plot(mean_star, dat$mu_test[[i]]);abline(0,1)
        # plot(cov_star, dat$sigma_test[[i]]); abline(0,1)
      }
      # Only use observed values when computing the probability (MAR assumption)
      obs_ind <- !is.na(S_test[i, ])
      S_test_obs <- S_test[i, obs_ind]
      mean_obs <- mean_star[obs_ind]
      cov_obs <- cov_star[obs_ind,obs_ind]
      if( sum(is_binary)==P ){ # All binary
        # Set lower and upper limits depending on whether S is 0/1; get binary prob.
        low_lims <- sapply(S_test_obs, function(x) ifelse(x==0,-Inf,0))
        upp_lims <- sapply(S_test_obs, function(x) ifelse(x==0,0,Inf))
        # https://stackoverflow.com/questions/51290014/rcpp-implementation-of-mvtnormpmvnorm-slower-than-original-r-function
        pi_SgivenY[i,c] <- as.numeric(pmvnorm(lower=low_lims, upper=upp_lims, mean=c(mean_obs), sigma=cov_obs))
      } else{ # all non-binary
        pi_SgivenY[i,c] <- as.numeric(dmvnorm(S_test_obs, mean=c(mean_obs), sigma=cov_obs))
      }
    } # for( i in 1:N_test )
  } # for( c in 1:num_causes )
  
  return(pi_SgivenY)
}
