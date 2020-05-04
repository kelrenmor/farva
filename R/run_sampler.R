farva_run <- function(S_mat, X_all_mu=NULL, X_all_sig=NULL, 
                      S_test=NULL, X_test_mu=NULL, X_test_sig=NULL, impossible_cause=NULL,
                      save_num_tot=500, burnin=1000, thin=10, L=5, K=5, inference=F,
                      mu_collapse=F, mc_tot=200, fast_test_samp=T, prec_equal_deltaTheta=T,
                      print_prog=T, save_inds_mu="unique", save_inds_sig="unique", return_data=T){
  # - S_mat (required) should be a num_causes list of matrices with TRAINING data with indexing
  #   from 1 to the number of training causes, i.e. { S_mat[[1]], ..., S_mat[[C]] }. 
  #   Alternatively, S_mat can me a matrix with the first column being a unique COD identifier
  #   and each row corresponding to an individual COD.
  #   In this case, X_all_mu and/or X_all_sig (if not null) must also be matrices with 
  #   rows corresponding to the equivalent row of S_mat.
  # - X_all_mu and X_all_sig should be length num_causes list of matrices (see above for exception)
  #   and will be auto-filled to just matrices of 1s if left as NULL.
  # - S_test should be a matrix with TEST data. If no test data given, only S_mat modeled.
  # - X_test_mu and X_test_sig should be matrices with covariates to be used for mean/cov
  #   and will be auto-filled to just matrices of 1s if left as NULL.
  # - impossible_cause (if specified) is a N_test x num_causes boolean matrix with a TRUE in the i x c cell 
  #   if test person i CANNOT have died of cause c (e.g., person i is female, cause c is prostate cancer).
  # - save_num_tot is the total number of saved sample desired...
  # - burnin is the number of samples to take before beginning saving...
  # - thin says how many samples to take for which 1 is saved...
  #   thus nsamps is the ( total number of samples taken = burnin + thin * save_num_tot ).
  # - L and K are the dimension reduction parameters (see paper).
  # - mu_collapse controls whether to use \z_i = \mu + \lam\eta (F) or \z_i = \lam\eta (T, default).
  # - mc_tot is the number of Monte carlo samples to be used for probability calc of test data (default 200).
  # - fast_test_samp (default T) says whether to use MC samples (faster) for test data even when data are not mixed.
  # - prec_equal_deltaTheta (default T) says whether to assume phi_delta = phi_theta and delta_delta = delta_theta
  # - burnin is number of samples to be discarded as burnin (default 1/5 of nsamps).
  # - thin is the total number of samples per one saved sample (default 10).
  # - save_inds_mu and save_inds_sig specify which of cov_all, mean_all, omega_all, eta_all should be 
  #   saved, "unique" (default) for unique cov combos, "all" for all (slowest), "first" for first
  
  ############ SETUP (PACKAGES) ############
  library(Rcpp)
  library(RcppArmadillo)
  library(RcppDist)
  library(mvtnorm)
  
  ############ ONE-TIME CALCULATIONS/TRANSFORMS ############
  
  if( is.null(X_all_sig) ){ 
    L = K 
    xi_identity = TRUE 
  } else if( dim(X_all_sig[[1]])[2] == 1 ){
    L = K
    xi_identity = TRUE 
  } else{ xi_identity = FALSE }
  
  # If S_mat is entered in matrix form, convert it to list form used by code.
  if( !is.list(S_mat) ){
    nullmu = is.null(X_all_mu)
    nullsig = is.null(X_all_sig)
    N_sum = nrow(S_mat)
    if('COD' %in% colnames(S_mat)){
      ind_rm = which(colnames(S_mat)=='COD')
    } else{
      ind_rm = 1
      print("No column named 'COD', so assuming first column of S_mat is COD.")
    }
    if( !nullmu ){if( !(N_sum==nrow(X_all_mu)) ){stop("Need nrow(X_all_mu) to equal nrow(S_mat)")}}
    if( !nullsig ){if( !(N_sum==nrow(X_all_sig)) ){stop("Need nrow(X_all_sig) to equal nrow(S_mat)")}}
    cods = S_mat[,ind_rm]
    un_cods = sort( unique(cods) )
    num_causes = length(un_cods)
    S_mat_tmp = list()
    if( !nullmu ){X_all_mu_tmp = list()}
    if( !nullsig ){X_all_sig_tmp = list()}
    for(c in 1:num_causes){
      ind_cod = which(cods == un_cods[c])
      S_mat_tmp[[c]] = S_mat[ind_cod,-ind_rm, drop=F]
      if( !nullmu ){X_all_mu_tmp[[c]] = X_all_mu[ind_cod,,drop=F]}
      if( !nullsig ){X_all_sig_tmp[[c]] = X_all_sig[ind_cod,,drop=F]}
    }
    S_mat = S_mat_tmp # replace S_mat with list version!
    if( !nullmu ){X_all_mu = X_all_mu_tmp} # replace X_all_mu with list version!
    if( !nullsig ){X_all_sig = X_all_sig_tmp} # replace X_all_sig with list version!
  } else{
    un_cods = NULL
  }
  
  # Get dimensions of things
  N <- lapply(S_mat,nrow)
  num_causes <- length(N)
  N_sum <- Reduce("+", N) 
  P <- ncol(S_mat[[1]])
  if( is.null(X_all_mu) ){ 
    X_all_mu = list()
    for(c in 1:num_causes){ X_all_mu[[c]] <- matrix(1,nrow=N[[c]],ncol=1) }
  }
  if( is.null(X_all_sig) ){ 
    X_all_sig = list()
    for(c in 1:num_causes){ X_all_sig[[c]] <- matrix(1,nrow=N[[c]],ncol=1) }
  }
  B_mu <- dim(X_all_mu[[1]])[2]; cov_incl = (B_mu>1) 
  B_sig <- dim(X_all_sig[[1]])[2]
  
  # Data structure (figure out which variables are BINARY)
  is_binary <- rep(T,P)
  tmp_mat <- do.call(rbind, S_mat);
  tmp_mat[is.na(tmp_mat)] <- 0; # for cleanliness, make NAs 0
  for(p in 1:P){ if( length(unique(tmp_mat[,p]))>2 ){ is_binary[p] <- F } }
  if( !is.matrix(S_mat[[1]]) ){ S_mat = lapply(S_mat, as.matrix) }
  if( !is.null(S_test) ){ if( !is.matrix(S_test) ){ S_mat = as.matrix(S_test) } }
  # Note for non-binary data, the init code sets z_{ij} = s_{ij}
  
  ############ DEFINE HYPER-HYPER PARAMETERS ############
  
  # Lam0_mult adjusts the prior on \alpha_{\mu}, \beta_{\mu} ~ N(0, Lam0_mult * I).
  # S0_mult and nu_add adjust the prior on \alpha_{\Sigma}, \beta_{\Sigma} ~ IW(B + 2 + nu_add, S0_mult * I).
  # Fix these here, but users can manually change in code if they need to.
  Lam0_mult=1; S0_mult=1; nu_add=0
  
  # Assign priors for the population-level mean and covariance for the 
  # coefficient vectors associated with entries of the matrices \xi_c 
  # and \psi_c (c in 1...C)
  mu_0_beta <- matrix(0,B_sig) # B_sig x 1 vector
  Lambda_0_beta <- Lam0_mult * diag(1,B_sig) # B_sig x B_sig matrix
  v0_beta <- B_sig + 2 + nu_add # integer
  S0_beta <- (v0_beta - B_mu - 1) * S0_mult * diag(1,B_sig) # B_sig x B_sig matrix
  mu_0_alpha <- mu_0_gamma <- matrix(0,B_mu) # B_mu x 1 vector
  Lambda_0_alpha <- Lambda_0_gamma <- Lam0_mult * diag(1,B_mu) # B_mu x B_mu matrix
  v0_alpha <- v0_gamma <- B_mu + 2 + nu_add # integer
  S0_alpha <- S0_gamma <- (v0_alpha - B_mu - 1) * S0_mult * diag(1,B_mu) # B_mu x B_mu matrix
  a_xi <- a_psi <- 1 # positive real
  b_xi <- b_psi <- 1 # positive real
  # Assign priors for the population-level mean and shrinkage prior for the 
  # \Theta_c and \Delta parameters.
  nu_theta <- nu_delta <- 3 # positive real (what Evan used in fixedfact.R)
  a1_del_del <- a1_del_the <- 2.1 # positive real (what Evan used in fixedfact.R)
  a2_del_del <- a2_del_the <- 3.1 # positive real (what Evan used in fixedfact.R)
  # Assign prior for the data covariance
  a_sigsq <- 1; b_sigsq <- 1
  
  ############ INITIALIZE THINGS ############
  
  Xt_all_mu <- Xt_all_sig <- list()
  XtX_all_mu <- XtX_all_sig <- XXt_all_sig <- list()
  for(c in 1:num_causes){
    # Pre-generate X matrices for mean components
    Xt_all_mu[[c]] <- t(X_all_mu[[c]])
    XtX_all_mu[[c]] <- as.matrix(Xt_all_mu[[c]] %*% X_all_mu[[c]])
    # Pre-generate X matrices for covariance components
    Xt_all_sig[[c]] <- t(X_all_sig[[c]])
    XtX_all_sig[[c]] <- as.matrix(Xt_all_sig[[c]] %*% X_all_sig[[c]])
    XXt_all_sig[[c]] <- array(NA, dim=c(B_sig, B_sig, N[[c]]))
    for(i in 1:N[[c]]){ XXt_all_sig[[c]][,,i] = Xt_all_sig[[c]][,i,drop=FALSE] %*% X_all_sig[[c]][i,,drop=FALSE] }
  }
  
  # Make list of save_inds_mu (if not provided) specifying individual-specific indices to save
  if( (save_inds_mu=="all")[1] ){
    save_inds_mu=list(); for(c in 1:num_causes){ save_inds_mu[[c]]=1:N[[c]] }
  } else if( (save_inds_mu=="unique")[1] ){
    save_inds_mu=list(); for(c in 1:num_causes){ save_inds_mu[[c]]=which(!duplicated(X_all_mu[[c]])) }
  } else if( (save_inds_mu=="first")[1] ){
    save_inds_mu=list(); for(c in 1:num_causes){ save_inds_mu[[c]]=1 }
  } else if( is.list(save_inds_mu) ){
    if( !(length(save_inds_mu)==num_causes) ){
      print("Need save_inds_mu to be 'all', 'unique', 'first', or num_causes list of indices.")
    }
  } else{
    print("Need save_inds_mu to be 'all', 'unique', 'first', or num_causes list of indices.")
  }
  # Make list of save_inds_sig (if not provided) specifying individual-specific indices to save
  if( (save_inds_sig=="all")[1] ){
    save_inds_sig=list(); for(c in 1:num_causes){ save_inds_sig[[c]]=1:N[[c]] }
  } else if( (save_inds_sig=="unique")[1] ){
    save_inds_sig=list(); for(c in 1:num_causes){ save_inds_sig[[c]]=which(!duplicated(X_all_sig[[c]])) }
  } else if( (save_inds_sig=="first")[1] ){
    save_inds_sig=list(); for(c in 1:num_causes){ save_inds_sig[[c]]=1 }
  } else if( is.list(save_inds_sig) ){
    if( !(length(save_inds_sig)==num_causes) ){
      print("Need save_inds_sig to be 'all', 'unique', 'first', or num_causes list of indices.")
    }
  } else{
    print("Need save_inds_sig to be 'all', 'unique', 'first', or num_causes list of indices.")
  }
  
  # Initialize parameters using set hyperparameter values
  init_list = run_sampler_init(S_mat, L, K, P, N, num_causes, B_sig, B_mu, 
                               a1_del_del, a2_del_del, a1_del_the, a2_del_the,
                               nu_delta, nu_theta, a_sigsq, b_sigsq,
                               mu_collapse, is_binary)
  list2env(init_list, environment()) # puts list elements in environment
  # For binary data (s_{ij} \in {0,1}), z_{ij} is also initialized
  
  ############ SETUP POST MEAN SAVES AND TEST COD/CSMF SAVES ############
  
  # Determine how many samples are to be saved through sampler
  # nsamps is the total number of samples (including burnin, thinned, etc) from which saved samples are pulled.
  nsamps <- burnin+save_num_tot*thin
  print(paste("Total number of samples will be:", nsamps))
  print(paste("Total saved samples will be:", save_num_tot))
  
  # Set up values to be used in the sampler loop for saving things
  save_num <- 1 # initialize save number
  if( !is.null(S_test) ){
    N_test <- nrow(S_test)
    if( is.null(X_test_mu) ){ X_test_mu <- matrix(1, nrow=N_test, ncol=1) }
    if( is.null(X_test_sig) ){ X_test_sig <- matrix(1, nrow=N_test, ncol=1) }
    cod_test_save <- matrix(NA, nrow=N_test, ncol=save_num_tot)
    csmf_test_save <- matrix(NA, nrow=num_causes, ncol=save_num_tot)
    indiv_prob <- matrix(0, nrow=N_test, ncol=num_causes)
  }
  
  # Set up posterior mean save variables
  postSavesnum <- 0
  post_saves_list <- postSavesInit(num_causes,P,K,L,save_inds_sig,save_inds_mu)
  list2env(post_saves_list, environment()) # puts list elements in environment
  if(inference){ 
    post_saves_inf_list <- postSavesInitInference(num_causes,P,L,save_num_tot,
                                                  save_inds_sig,save_inds_mu)
    list2env(post_saves_inf_list, environment()) # puts list elements in environment
  }
  
  ############ RUN SAMPLER ############
  ten_perc <- round(seq(1,nsamps,length=10))
  for(ss in 1:nsamps){
    
    if( ss%in%ten_perc & print_prog ){ print(paste("sample",ss,"of",nsamps)) }
    print(paste("sample",ss,"of",nsamps))
    
    ############ UPDATE \mu_{c,j} ############
    ############ and hierarchical regression params for each
    
    if( mu_collapse ){ # If we have z_i=\Lambda\eta_i, \eta_i~N(\psi x, I), fix mu_all[[c]] to be 0-vec.
      
      if( ss==1 ){ # Only need to set once, then it just stays the same!
        for(c in 1:num_causes){
          mu_all[[c]] <- matrix(0,nrow=P,ncol=N[[c]])
        }
      }
      
    } else{ # If we have z_i = \mu_{c[i]} + \Lambda\eta_i, \eta_i~N(0, I), sample mu_all[[c]].
      
      # Sample \mu vector for each cause c
      for(c in 1:num_causes){
        for(p in 1:P){
          for(i in 1:N[[c]]){
            mu_all[[c]][p,i] <- matrix(X_all_mu[[c]][i,], nrow=1) %*% matrix(gamma_c_all[[c]][,p])
          }
        }
      }
      
      # Sample regression parameters gamma_c_all for each cause c
      for(c in 1:num_causes){
        # Get temporary matrix to hold observed deviations of z_i from \lambda \eta_i
        dev_tmp = matrix(NA,nrow=N[[c]],ncol=P)
        for(i in 1:N[[c]]){
          dev_tmp[i,] <- c(z_all_nominus[[c]][,i] - Omega_all[[c]][,,i] %*% matrix(psi_all[[c]][,i]))
        }
        for(p in 1:P){
          gamma_c_all[[c]][,p] <- sample_beta_c(matrix(dev_tmp[,p]), matrix(gamma_all[,p]), as.matrix(gamma_Sigma[,,p]), 
                                                Sigma_0_vec[p], XtX_all_mu[[c]], Xt_all_mu[[c]])
        }
      }
      
      # Sample \gamma_all_{p} (population-level mean) for each entry p\in{1,..,P};
      # the result is a length B_mu vector with entry b being the coefficient for covariate b
      gamma_means <- Reduce("+", gamma_c_all) / length(gamma_c_all) # BxP array of avg of cause-spec betas
      for(p in 1:P){
        gamma_all[,p] <- sample_beta_mu(mu_0_gamma, Lambda_0_gamma, num_causes, 
                                        matrix(gamma_means[,p]), as.matrix(gamma_Sigma[,,p]))
      }
      
      # Sample \gamma_Sigma{p} (population-level covariance) for each entry p\in{1,..,P};
      # the result is a B_mu x B_mu covariance matrix with entry b1, b2 being the cov between coef b1 and b2
      for(p in 1:P){
        gamma_Sigma[,,p] <- sample_beta_Sigma(v0_gamma, S0_gamma, num_causes, matrix(gamma_all[,p]), 
                                              do.call(cbind, 
                                                      rapply(gamma_c_all, classes='matrix', how='list', 
                                                             f=function(x) x[, p, drop=FALSE])) )
      }
      
    }
    
    ############ UPDATE \xi_{c,lk} ############
    ############ and hierarchical regression params for each
    
    if( xi_identity ){
      if( ss==1 ){
        # Fix xi to the identity matrix (only need to be done once)
        I_k = diag(K)
        for(c in 1:num_causes){
          for(i in 1:N[[c]]){
            xi_all[[c]][,,i] = I_k
          }
        }
      }
      # Update Omega_all ( \Omega_i = \Theta \xi_{c[i]}(x_i) for i in 1,...,N_c )
      for(c in 1:num_causes){
        for(i in 1:(N[[c]]) ){
          Omega_all[[c]][,,i] <- Theta_all[[c]]
        }
      }
    } else{
      # Sample \beta_{c,lk} for each cause c and entry l\in{1,..,L} and k\in{1,...,K};
      # the result is a length B_sig vector with entry b being the coefficient for covariate b for \xi_{c,lk}
      for(c in 1:num_causes){
        for(l in 1:L){
          for(k in 1:K){
            beta_c_all[[c]][,l,k] <- sample_betaxi_c(z_c=z_all[[c]], Sig0vec=Sigma_0_vec, 
                                                     eta_c=eta_all[[c]], Theta_c=Theta_all[[c]], beta_c=beta_c_all[[c]],
                                                     mu_beta=matrix(beta_all[,l,k]), Sigma_beta=as.matrix(beta_Sigma[,,l,k]), 
                                                     XXt_c=XXt_all_sig[[c]],  X_c=Xt_all_sig[[c]],
                                                     k_get=k-1, l_get=l-1)
            xi_all[[c]][l,k,] = (X_all_sig[[c]])%*%matrix(beta_c_all[[c]][,l,k])
          }
        }
      }
      
      # Rescale Xi (i.e., beta_c_all) and Theta/Delta, along with variance/precision
      # parameters for each, so as to fix multiplicative identifiability issue.
      if( ss==1 ){ # only need to calculate this once
        I_k = diag(K); nI_k = norm(I_k, type='F')
        nI = sum(unlist(N)) * nI_k
      }
      nXi = 0
      for(c in 1:num_causes){ 
        for(i in 1:N[[c]]){ 
          nXi = nXi + norm(xi_all[[c]][,,i], type='F') 
        } 
      }
      mult = nI / nXi
      for(c in 1:num_causes){
        for(l in 1:L){
          for(k in 1:K){
            beta_c_all[[c]][,l,k] <- mult * beta_c_all[[c]][,l,k]
            xi_all[[c]][l,k,] = mult * xi_all[[c]][l,k,]
          }
        }
        Theta_all[[c]] = Theta_all[[c]] / mult
      }
      phi_theta = mult * mult * phi_theta # precision term, so x~N(0,phi^-1) --> x/m~N(0,(m^2phi)^-1)
      phi_delta = mult * mult * phi_delta
      Delta = Delta / mult
      for(l in 1:L){
        for(k in 1:K){
          beta_all[,l,k] <- mult * beta_all[,l,k]
          beta_Sigma[,,l,k] <- mult * mult * beta_Sigma[,,l,k]
        }
      }
      
      # Sample \beta_{lk} (population-level mean) for each entry l\in{1,..,L} and k\in{1,...,K};
      # the result is a length B_sig vector with entry b being the coefficient for covariate b
      beta_means <- Reduce("+", beta_c_all) / length(beta_c_all) # BxLxK array of avg of cause-spec betas
      for(l in 1:L){
        for(k in 1:K){
          beta_all[,l,k] <- sample_beta_mu(mu_0_beta, Lambda_0_beta, num_causes, 
                                           beta_means[,l,k,drop=F], as.matrix(beta_Sigma[,,l,k]))
        }
      }
      
      # Sample \beta_Sigma{lk} (population-level covariance) for each entry l\in{1,..,L} and k\in{1,...,K};
      # the result is a B_sig x B_sig covariance matrix with entry b1, b2 being the cov between coef b1 and b2
      for(l in 1:L){
        for(k in 1:K){
          beta_Sigma[,,l,k] <- sample_beta_Sigma(v0_beta, S0_beta, num_causes, beta_all[,l,k,drop=F], 
                                                 do.call(cbind, 
                                                         rapply(beta_c_all, classes='array', how='list', 
                                                                f=function(x) x[, l, k, drop=FALSE])) )
        }
      }
      
      # Update Omega_all ( \Omega_i = \Theta \xi_{c[i]}(x_i) for i in 1,...,N_c )
      for(c in 1:num_causes){
        for(i in 1:(N[[c]]) ){
          Omega_all[[c]][,,i] <- get_Omega_i(Theta_all[[c]], xi_all[[c]][,,i])
        }
      }
      
    }
    
    ############ UPDATE \psi_{c,k} ############
    ############ and hierarchical regression params for each
    
    # Fix sigSqpsi_all to 1 for identifiabilitys
    if( ss==1 ){ for(k in 1:K){ sigSqpsi_all[k] <- 1 } }
    
    # Sample \psi_{c,i,k} for entries k\in{1,...,K} for each death i with cause c ;
    # each iter samples a single numer, the draw for component k of death i of cause c
    for(c in 1:num_causes){
      for(i in 1:N[[c]]){
        psi_all[[c]][,i] <- sample_psi_i(sigSqpsi_all, alpha_c_all[[c]], matrix(X_all_mu[[c]][i,]),
                                         Sigma_0_vec, Omega_all[[c]][,,i], matrix(z_all[[c]][,i]))
      }
    }
    
    if( mu_collapse ){ # If we have z_i=\Lambda\eta_i, \eta_i~N(\psi x, I) 
      
      # Sample \alpha_{c,k} for each cause c and entry k\in{1,...,K};
      # \alpha_{c,k} is length B_mu vector with entry b the coefficient for covariate b for \psi_{c,k}
      psi_RSS <- matrix(0,nrow=K)
      for(c in 1:num_causes){
        for(k in 1:K){
          alpha_c_all[[c]][,k] <- sample_beta_c(matrix(psi_all[[c]][k,]), matrix(alpha_all[,k]), 
                                                as.matrix(alpha_Sigma[,,k]), sigSqpsi_all[k], 
                                                XtX_all_mu[[c]], Xt_all_mu[[c]])
        }
      }
      
      # Sample \alpha_{k} (population-level mean) for each entry k\in{1,...,K};
      # the result is a length B_mu vector with entry b being the coefficient for covariate b
      alpha_means <- Reduce("+", alpha_c_all) / length(alpha_c_all) # BxK matrix of avg of cause-spec betas
      for(k in 1:K){
        alpha_all[,k] <- sample_beta_mu(mu_0_alpha, Lambda_0_alpha, num_causes, 
                                        matrix(alpha_means[,k]), as.matrix(alpha_Sigma[,,k]))
      }
      
      # Sample \alpha_Sigma{lk} (population-level covariance) for each entry k\in{1,...,K};
      # the result is a B_mu x B_mu covariance matrix with entry b1, b2 being the cov between coef b1 and b2
      for(k in 1:K){
        alpha_Sigma[,,k] <- sample_beta_Sigma(v0_alpha, S0_alpha, num_causes, alpha_all[,k], 
                                              do.call(cbind, 
                                                      rapply(alpha_c_all, classes='matrix', how='list', 
                                                             f=function(x) x[, k, drop=FALSE])) )
      }
      
    } else{ # If we have z_i = \mu + \Lambda\eta_i, \eta_i~N(0, I) 
      
      # Fix alpha_c_all to 0, and alpha_Sigma to the identity, then you end up with \eta_i~N(0, I)
      if(ss==1){ 
        for(c in 1:num_causes){
          for(k in 1:K){
            alpha_c_all[[c]][,k] <- matrix(0,nrow=nrow(XtX_all_mu[[c]]),ncol=1) # allows for mean-0 spec of eta
          }
        }
        for(k in 1:K){ alpha_Sigma[,,k] <- diag(1,B_sig) } 
      }
      
    }
    
    ############ UPDATE \eta_{i} ############
    ############ (via sampling \nu_{i})
    
    # Repsent latent factor process \eta_{i} equivalently as
    # \eta_{i} = \psi(x_i) + \nu_i, where \nu_i~N(0,I_k).
    # BUT I'm going to fix \nu_i to 0, and model \eta_{i,k} ~ (\beta_k x_i, 1) for k \in 1,...,K.
    for(c in 1:num_causes){
      for(i in 1:N[[c]]){
        # Fix nu to 0, based on my adjusted parameterization
        nu_all[[c]][,i] <- 0
        # Now \eta_{i} is simply the sum of psi and nu!
        eta_all[[c]][,i] <- psi_all[[c]][,i] + nu_all[[c]][,i]
      }
    }
    
    ############ UPDATE \sigma^2_{p} ############
    ############ (for non-binary symptoms)
    
    # Get RSS for each 
    sigsq_RSS <- matrix(0,nrow=P)
    for(c in 1:num_causes){
      for(i in 1:N[[c]]){
        tmp_rs <- (Theta_all[[c]]) %*% (xi_all[[c]][,,i]) %*% matrix(eta_all[[c]][,i]) - 
          matrix(z_all[[c]][,i])
        sigsq_RSS = sigsq_RSS + tmp_rs^2 
      }
    }
    
    for(p in 1:P){
      if( !is_binary[p] ){
        Sigma_0[p,p] <- sample_sigsq(a_sigsq, b_sigsq, N_sum, sigsq_RSS[p])
      }
    }
    # Make Sigma_0_vec from updated Sigma_0
    Sigma_0_vec <- matrix(diag(Sigma_0))
    
    ############ UPDATE \Theta_{c,p*} ############
    ############ and shrinkage params for each
    
    # Sample each row of PxL matrix \Theta_c for each c
    for(c in 1:num_causes){
      for(p in 1:P){
        Theta_all[[c]][p,] <- sample_Theta_j(Sigma_0_vec[p], xi_all[[c]], eta_all[[c]], 
                                             matrix(z_all[[c]][p,]), matrix(Delta[p,]),
                                             matrix(phi_theta[p,]), tau_theta)
      }
    }
    
    if( !prec_equal_deltaTheta ){    
      # Sample each element of \phi_{\theta}
      for(p in 1:P){
        for(l in 1:L){
          sum_tmp <- 0
          for(c in 1:num_causes){
            sum_tmp <- sum_tmp + (Theta_all[[c]][p,l]-Delta[p,l])^2
          }
          phi_theta[p,l] <- sample_phi_pl(nu_theta, tau_theta[l] * sum_tmp, num_causes)
        }
      }
  
      # Sample each element of L-vec \delta_{\theta} and calculate tau_theta (cum prod) on the way
      Theta_cube <- array(unlist(Theta_all),
                          dim = c(nrow(Theta_all[[1]]), ncol(Theta_all[[1]]), length(Theta_all)))
      for(l in 1:L){
        delta_theta[l] <- sample_delta_theta(a1_del_the, a2_del_the, Delta, Theta_cube,
                                             phi_theta, delta_theta, l-1)
        if(l==1){
          # initialize tau_theta first value
          tau_theta[l] <- delta_theta[l]
        } else{
          # update tau to be product of deltas as you go
          tau_theta[l] <- tau_theta[l-1]*delta_theta[l]
        }
      }
    }
    
    ############ UPDATE \Delta_{c,pl} ############
    ############ and shrinkage params for each
    
    # Sample each row of PxL matrix \Delta
    for(p in 1:P){ 
      for(l in 1:L){
        Delta[p,l] <- sample_Delta_jl(phi_delta[p,l], tau_delta[l],
                                      matrix(unlist( rapply(Theta_all, classes='matrix', how='list', 
                                                            f=function(x) x[p,l, drop=FALSE]) )), 
                                      phi_theta[p,l], tau_theta[l])
      }
    }
    
    if( prec_equal_deltaTheta ){ # If Delta shares its precision parameters with Theta
      
      # Sample each element of \phi_{\theta} (use this for \phi_{\delta} too)
      for(p in 1:P){
        for(l in 1:L){
          sum_tmp <- tau_theta[l]*(Delta[p,l])^2
          for(c in 1:num_causes){
            sum_tmp <- sum_tmp + tau_theta[l]*(Theta_all[[c]][p,l]-Delta[p,l])^2
          }
          sample_phi_pl(nu_theta, tau_theta[l] * sum_tmp, num_causes+1)
        }
        phi_delta[p,] <- phi_theta[p,]
      }  
      
      # Sample each element of L-vec \delta_{\theta} and calculate tau_theta (cum prod) on the way
      Theta_cube <- array(c(unlist(Theta_all),rep(0,P*L)), 
                          dim = c(nrow(Theta_all[[1]]), ncol(Theta_all[[1]]), length(Theta_all) + 1))
      for(l in 1:L){
        delta_theta[l] <- sample_delta_theta(a1_del_the, a2_del_the, Delta, Theta_cube,
                                             phi_theta, delta_theta, l-1)
        delta_delta[l] <- delta_theta[l]
        if(l==1){
          # initialize tau_theta first value
          tau_theta[l] <- delta_theta[l]
          tau_delta[l] <- tau_theta[l]
        } else{ 
          # update tau to be product of deltas as you go
          tau_theta[l] <- tau_theta[l-1]*delta_theta[l]
          tau_delta[l] <- tau_theta[l]
        }
      }
      
    } else{ # If Delta has its own unique precision parameters
      
      # Sample each element of PxL matrix \phi_{\Delta}
      for(p in 1:P){
        for(l in 1:L){
          phi_delta[p,l] <- sample_phi_pl(nu_delta, tau_delta[l] * (Delta[p,l]^2), 1)
        }
      }

      # Sample each element of L vector \delta_{\Delta} and calculate tau_delta (cum prod) on the way
      for(l in 1:L){
        delta_delta[l] <- sample_delta_Delta(a1_del_del, a2_del_del, Delta, phi_delta, delta_delta, l-1)
        if(l==1){
          tau_delta[l] <- delta_delta[l]
        } else{ # update to be product as you go
          tau_delta[l] <- tau_delta[l-1]*delta_delta[l]
        }
      }
      
    }
    
    ############ UPDATE z_{ij} ############
    ############ (for binary symptoms)
    
    # For observations s_{ij} with symptom j binary,
    # we need to sample to get z_{ij} at each iteration.
    # (For non-binary symptom j, z_{ij} = s_{ij} always.)
    
    # Update Omega_all ( \Omega_i = \Theta \xi_{c[i]}(x_i) for i in 1,...,N_c )
    for(c in 1:num_causes){
      for(i in 1:(N[[c]]) ){
        Omega_all[[c]][,,i] <- get_Omega_i(Theta_all[[c]], xi_all[[c]][,,i])
      }
    }
    
    # Sample each element of z_all and calculate mean and covariance on the way!
    samps <- sample_z_mean_cov_all(S_mat, Omega_all, eta_all, unlist(N), 
                                   Sigma_0_vec, z_all_nominus, z_all, P,
                                   is_binary, mu_all, mu_collapse, psi_all)
    mean_all = samps[["mean_all"]]
    cov_all = samps[["cov_all"]]
    z_all = samps[["z_all"]]
    z_all_nominus = samps[["z_all_nominus"]]
    
    # Want to keep Inf/-Inf values from breaking my sampler, but print warnings as this is used
    bad_z = 0; z_all_sum = sum(unlist(z_all))
    if( is.infinite(z_all_sum) | is.na(z_all_sum) ){
      for(c in 1:num_causes){
        for(i in 1:(N[[c]]) ){
          for(j in 1:P){
            if(is_binary[j]){
              bad_cond = is.infinite(z_all[[c]][j,i]); bad_replace = 50
              if(bad_cond){
                if( S_mat[[c]][i,j] == 1 ){ z_all[[c]][j,i] <- bad_replace; bad_z = bad_z+1 }
                if( S_mat[[c]][i,j] == 0 ){ z_all[[c]][j,i] <- -bad_replace; bad_z = bad_z+1 }
              }
            }
          }
        }
      }
      print(paste("bad_z=",bad_z,sep=""))
    }
    
    ############ CLEANUP AND SAVE, AND SAMPLE FOR Y_TEST ############
    ############ 
    
    # Save and add things to running mean total depending on conditions
    if( (ss%%thin==0) & (ss>burnin) & (save_num<=save_num_tot) ){
      
      # Save the current values of the higher level parameters to get posterior mean
      postSavesnum <- postSavesnum + 1
      # Collect ongoing sum of select parameters to get posterior mean
      Omega_all_post <- Map("+", Omega_all_post, lapply(1:num_causes, function(x) Omega_all[[x]][, , save_inds_sig[[x]], drop=FALSE]))
      eta_all_post <- Map("+", eta_all_post, lapply(1:num_causes, function(x) eta_all[[x]][, save_inds_mu[[x]], drop=FALSE]))
      sigsq_RSS_post <- sigsq_RSS_post + sigsq_RSS
      Sigma_0_vec_post <- Sigma_0_vec_post + Sigma_0_vec
      mean_all_post <- Map("+", mean_all_post, lapply(1:num_causes, function(x) mean_all[[x]][, save_inds_mu[[x]], drop=FALSE]))
      cov_all_post <- Map("+", cov_all_post, lapply(1:num_causes, function(x) cov_all[[x]][, , save_inds_sig[[x]], drop=FALSE]))
      Theta_all_post <- Map("+", Theta_all_post, Theta_all)
      Delta_post <- Delta_post + Delta
      tau_delta_post <- tau_delta_post + tau_delta
      tau_theta_post <- tau_theta_post + tau_theta
      if(inference){ 
        # Save current parameter values for inference-parameters
        for(c in 1:num_causes){
          mean_all_inf[[c]][,,save_num] <- mean_all[[c]][, save_inds_mu[[c]]] # , drop=FALSE
          cov_all_inf[[c]][,,,save_num] <- cov_all[[c]][, , save_inds_sig[[c]]] # , drop=FALSE
        }
        tau_delta_inf[,save_num] <- tau_delta
        tau_theta_inf[,save_num] <- tau_theta
        Sigma_0_vec_inf[,save_num] <- Sigma_0_vec
      }
      
      ############ SAMPLE y_{i*} FOR TEST DATA ############
      ############ (if you have people with unknown COD)
      
      if( !is.null(S_test) ){
        if( mu_collapse ){ gamma_c_all=list(); for(c in 1:num_causes){gamma_c_all[[c]]=matrix(0,nrow=1,ncol=P)} }
        if( fast_test_samp | !(sum(is_binary)==P | sum(!is_binary)==P) ){
          # If you have mixed data, or want faster sampling, then sample via MC approximation.
          pi_SgivenY <- get_piSgivenY(N_test, num_causes, P, mc_tot, cov_incl, X_test_mu, X_test_sig, S_test,
                                      beta_c_all, Theta_all, mu_collapse, gamma_c_all,
                                      alpha_c_all, matrix(sigSqpsi_all), Sigma_0, is_binary)
        } else{
          pi_SgivenY <- get_piSgivenY_analytic(sigSqpsi_all, beta_c_all, alpha_c_all, 
                                               X_test_sig, X_test_mu, Theta_all,
                                               P, num_causes, N_test, mu_collapse, gamma_c_all, 
                                               Sigma_0, S_test, is_binary)
        }
        if( !is.null(impossible_cause) ){ # Set impossible-cause entries of pi_SgivenY to 0.
          pi_SgivenY[impossible_cause] = 0
        }
        pi_SgivenY <- sweep(pi_SgivenY,1,rowSums(pi_SgivenY),`/`) # divide each column entry by the sum of its rows
        # sum(pi_SgivenY[1,]); round(pi_SgivenY,3)
        # cbind(cs_test, round(pi_SgivenY,3))
        
        # Now get individual CODs and CSMF
        cod_test_tmp <- rep(NA,N_test)
        csmf_tmp <- rep(NA,num_causes)
        for(i in 1:N_test){ cod_test_tmp[i] <- sampleDist(n=1, num_causes=num_causes, probs=pi_SgivenY[i,]) }
        for(c in 1:num_causes){ csmf_tmp[c] <- mean(cod_test_tmp==c) }
        indiv_prob <- indiv_prob + pi_SgivenY # replace(pi_SgivenY, is.na(pi_SgivenY), 0)
        cod_test_save[,save_num] <- cod_test_tmp
        csmf_test_save[,save_num] <- csmf_tmp
      }
      
      save_num = save_num+1
    }
    
  } # for(ss in 1:nsamps)
  # Get posterior mean from ongoing sum of select parameters
  Omega_all_post <- lapply(Omega_all_post, function(x) x/postSavesnum)
  eta_all_post <- lapply(eta_all_post, function(x) x/postSavesnum)
  sigsq_RSS_post <- sigsq_RSS_post/postSavesnum
  Sigma_0_vec_post <- Sigma_0_vec_post/postSavesnum
  mean_all_post <- lapply(mean_all_post, function(x) x/postSavesnum)
  cov_all_post <- lapply(cov_all_post, function(x) x/postSavesnum)
  if( !is.null(S_test) ){indiv_prob <- sweep(indiv_prob,1,rowSums(indiv_prob),`/`)}
  Theta_all_post <- lapply(Theta_all_post, function(x) x/postSavesnum)
  Delta_post <- Delta_post/postSavesnum
  tau_delta_post <- tau_delta_post/postSavesnum
  tau_theta_post <- tau_theta_post/postSavesnum
  
  if( is.null(S_test) ){ 
    csmf_test_save = cod_test_save = indiv_prob = NULL 
  } else{ 
    rownames(csmf_test_save) = un_cods 
    colnames(indiv_prob) <- un_cods
  }
  if( !inference ){ mean_all_inf = cov_all_inf = tau_delta_inf = tau_theta_post = Sigma_0_vec_inf = NULL }
  
  farva_res = list("csmf_test_save"=csmf_test_save, "cod_test_save"=cod_test_save,
                   "indiv_prob"=indiv_prob, 
                   "save_inds_mu"=save_inds_mu, "save_inds_sig"=save_inds_sig, 
                   "mean_all_post"=mean_all_post, "cov_all_post"=cov_all_post,
                   "Theta_all_post"=Theta_all_post, "tau_theta_post"=tau_theta_post,
                   "Delta_post"=Delta_post, "tau_delta_post"=tau_delta_post,
                   "K"=K,"L"=L,
                   "un_cods"=un_cods,
                   "impossible_cause"=impossible_cause, "burnin"=burnin, "thin"=thin,
                   "mu_collapse"=mu_collapse, "mc_tot"=mc_tot, "fast_test_samp"=fast_test_samp)
  
  if( return_data ){
    dat_list = list("S_mat"=S_mat, "X_all_mu"=X_all_mu, "X_all_sig"=X_all_sig, 
                    "S_test"=S_test, "X_test_mu"=X_test_mu, "X_test_sig"=X_test_sig)
  } else{
    dat_list = list()
  }
  
  if( inference ){ 
    inf_list = list("mean_all_inf"=mean_all_inf,"cov_all_inf"=cov_all_inf,
                    "tau_delta_inf"=tau_delta_inf, "tau_theta_inf"=tau_theta_inf,
                    "Sigma_0_vec_inf"=Sigma_0_vec_inf)
  } else{
    inf_list = list()
  }
  
  return( c(farva_res, dat_list, inf_list) )
}
