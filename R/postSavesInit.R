
postSavesInit <- function(num_causes,P,K,L,save_inds_sig,save_inds_mu){
  
  # Set up 0 matrices/arrays/vectors to hold saved parameter values
  Omega_all_post <- eta_all_post <- mean_all_post <- cov_all_post <- Theta_all_post <- list()
  for(c in 1:num_causes){
    l_sis <- length(save_inds_sig[[c]])
    l_sim <- length(save_inds_mu[[c]])
    Omega_all_post[[c]] <- array(0, c(P,K,z))
    eta_all_post[[c]] <- matrix(0, nrow=K, ncol=l_sim)
    mean_all_post[[c]] <- matrix(0, nrow=P, ncol=l_sim)
    cov_all_post[[c]] <- array(0, c(P,P,l_sis))
    Theta_all_post[[c]] <- matrix(0, nrow=P, ncol=L)
  }
  Delta_post <- matrix(0, nrow=P, ncol=L)
  tau_delta_post <- rep(0, L)
  tau_theta_post <- rep(0, L)
  sigsq_RSS_post <- Sigma_0_vec_post <- matrix(0,nrow=P,ncol=1)
  
  post_saves_list <- list('Omega_all_post'=Omega_all_post, 'eta_all_post'=eta_all_post, 
                          'mean_all_post'=mean_all_post, 'cov_all_post'=cov_all_post, 
                          'Theta_all_post'=Theta_all_post, 'Delta_post'=Delta_post, 
                          'tau_delta_post'=tau_delta_post, 'tau_theta_post'=tau_theta_post, 
                          'sigsq_RSS_post'=sigsq_RSS_post, 'Sigma_0_vec_post'=Sigma_0_vec_post)

  return(post_saves_list)
}
