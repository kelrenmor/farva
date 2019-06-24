postSavesInitInference <- function(num_causes,P,L,save_num_tot,save_inds_sig,save_inds_mu){

  # Set up 0 matrices/arrays/vectors to hold saved parameter values
  mean_all_inf <- cov_all_inf <- list()
  for(c in 1:num_causes){
    l_sis <- length(save_inds_sig[[c]])
    l_sim <- length(save_inds_mu[[c]])
    mean_all_inf[[c]] <- array(NA, c(P, l_sim, save_num_tot) )
    cov_all_inf[[c]] <- array(NA, c(P, P, l_sis, save_num_tot) )
  }
  tau_delta_inf <- tau_theta_inf <- matrix(NA, nrow=L, ncol=save_num_tot)
  Sigma_0_vec_inf <- matrix(NA, nrow=P, ncol=save_num_tot)

  post_saves_inf_list = list('mean_all_inf'=mean_all_inf, 'cov_all_inf'=cov_all_inf, 
                             'tau_delta_inf'=tau_delta_inf, 'tau_theta_inf'=tau_theta_inf, 
                             'Sigma_0_vec_inf'=Sigma_0_vec_inf)
  
  return(post_saves_inf_list)
}