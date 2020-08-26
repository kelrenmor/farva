postSavesInitInference <- function(num_causes,P,K,L,save_num_tot,save_inds_sig,save_inds_mu){
  
  # Set up 0 matrices/arrays/vectors to hold saved parameter values
  z_indiv_inf <- mean_indiv_inf <- cov_indiv_inf <- mu_inf <- list()
  for(c in 1:num_causes){
    l_sis <- length(save_inds_sig[[c]])
    l_sim <- length(save_inds_mu[[c]])
    z_indiv_inf[[c]] <- array(NA, c(P, l_sim, save_num_tot) )
    mu_inf[[c]] <- array(NA, c(P, l_sim, save_num_tot) )
    mean_indiv_inf[[c]] <- array(NA, c(P, l_sim, save_num_tot) )
    cov_indiv_inf[[c]] <- array(NA, c(P, P, l_sis, save_num_tot) )
  }
  tau_delta_inf <- tau_theta_inf <- matrix(NA, nrow=L, ncol=save_num_tot)
  Sigma_0_vec_inf <- matrix(NA, nrow=P, ncol=save_num_tot)
  
  post_saves_inf_list = list('z_indiv_inf'=z_indiv_inf, 'mu_inf'=mu_inf,
                             'mean_indiv_inf'=mean_indiv_inf, 'cov_indiv_inf'=cov_indiv_inf, 
                             'tau_delta_inf'=tau_delta_inf, 'tau_theta_inf'=tau_theta_inf, 
                             'Sigma_0_vec_inf'=Sigma_0_vec_inf)
  
  return(post_saves_inf_list)
}