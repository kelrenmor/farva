# Set up functions (inspiration and first two from Isaac)
simulate_lambda <- function(k, p, sigmasq){
  sapply(1:k, function(j) simulate_lambda_column(p, j))
}

simulate_lambda_column <- function(p, j){ # p is number of symptoms, j is column number being simulated
  library(tidyverse)
  value = runif(n = p, min = .5, max = 1) * sample(c(-1, 1), size = p, replace=TRUE)
  nonzero = rbernoulli(n = p, p = .4 + .2/j) # higher column number means more chance of 0s
  value[!nonzero] = 0 # replace the not-nonzero entries with 0
  return(value)
}

generate_mean <- function(p){
  library(truncnorm)
  # Try to mimic the mean I observe in real data, with lots of mostly-0s, 
  # some semi-common, and some super frequent 1s
  # hist(rtruncnorm(1000, a=-0.8, b=0.8, mean=-0.25, sd=0.3), breaks=30)
  # mean(rtruncnorm(1000, a=-0.8, b=0.8, mean=-0.25, sd=0.3)>0)
  c(truncnorm::rtruncnorm(p, a=-0.8, b=0.8, mean=-0.25, sd=0.3))
}

simulate_x_simABCD <- function(num_causes, N, p, k, sigmasq, sim, lambda_generator=simulate_lambda, 
                               N_test=NULL, bin=T, perturb_scale=1, cause_wt=0.5){
  library(mvtnorm)
  # num_causes is the number of causes to simulate
  # N should be a length num_causes list of the number of simulated observations desired for each cause
  # p is number of symptoms, k is underlying dimension of lambda, sigmaaq is noise variance
  # lambda_generator is the function used to make lambda 
  # N_test is a length num_causes list of the number of simulated "test" observations desired for each cause
  # bin is a boolean denoting whether resulting symptom observations should be binary
  # sim ("A","B","C",or "D") is a character letter denoting which simulation to run
  if( !(length(N)==num_causes) ){ stop('N must be length num_causes') }
  if( !(p>=k) ){ stop('Need p>=k') }
  if( !(sigmasq>0) ){ stop('Need sigmasq>0') }
  if( !(sim%in%c("A","B","C","D")) ){ stop('Need sim to be one of "A","B","C","D"') }
  lambda = lambda_generator(k, p)
  if( sim%in%c("C","D") ){
    perturb = map(1:num_causes, function(x) perturb_scale * lambda_generator(k, p))
  }
  sigma = mu = list()
  S_mat = X_all = list()
  if( !is.null(N_test) ){
    if( !(length(N_test)==num_causes) ){ stop('N_test must be length num_causes') }
    N_test_sum = Reduce("+", N_test)
    S_test = matrix(NA, nrow=N_test_sum, ncol=p) # Will fill in below.
    X_test = matrix(1, nrow=N_test_sum, ncol=1)  # sims A-D have no covariates, so just make 'intercept' only covariate.
    test_ind = 1
    cods_test_true = c()
  }
  cods_train_true = c()
  if( sim=="A" ) {
    # Common independent covariance
    common_cov = sigmasq*diag(p)
  } else if(sim=="B"){
    # Common dependent covariance
    common_cov = sigmasq*diag(p) + tcrossprod(lambda)
  } else if( sim=="C"|sim=="D" ){
    # Common mean
    common_mu = generate_mean(p)
  }
  for( cs in 1:num_causes ){
    cods_train_true = c(cods_train_true, rep(cs, N[[cs]]))
    if( sim=="A" ){
      # Cause-specific mean and common independent cov
      mu[[cs]] = generate_mean(p)
      sigma[[cs]] = common_cov
    } else if( sim=="B" ){
      # Cause-specific mean and common dependent cov
      mu[[cs]] = generate_mean(p)
      sigma[[cs]] = common_cov
    } else if( sim=="C"){
      # Common mean and cause-specific dependent cov
      mu[[cs]] = common_mu
      sigma[[cs]] = sigmasq*diag(p) + (1-cause_wt)*tcrossprod(lambda) + cause_wt*tcrossprod(perturb[[cs]])
    } else if( sim=="D" ){
      # Cause-specific mean and cause-specific dependent cov
      mu[[cs]] = (1-cause_wt)*common_mu + cause_wt*generate_mean(p) # add some "common" symptom mean levels
      sigma[[cs]] = sigmasq*diag(p) + (1-cause_wt)*tcrossprod(lambda) + cause_wt*tcrossprod(perturb[[cs]])
    }
    if( bin ){
      sigma[[cs]] = cov2cor(sigma[[cs]])
    }
    S_mat[[cs]] = rmvnorm(n=N[[cs]], mean = mu[[cs]], sigma = sigma[[cs]])
    if( bin ){
      S_mat[[cs]] = 1*(S_mat[[cs]]>0)
    }
    X_all[[cs]] = matrix(1, nrow=N[[cs]], ncol=1)
    if( !is.null(N_test) ){
      if( !(N_test[[cs]]==0) ){
        S_test[test_ind:(test_ind+N_test[[cs]]-1),] = rmvnorm(n=N_test[[cs]], mean = mu[[cs]], sigma = sigma[[cs]])
      }
      test_ind = test_ind+N_test[[cs]]
      cods_test_true = c(cods_test_true, rep(cs, N_test[[cs]]))
    }
  } # for( cs in 1:num_causes )
  if( !is.null(N_test) ){
    if( bin ){
      S_test = 1*(S_test>0)
    }
    csmf_test_true = unlist(lapply(N_test,function(l) l/N_test_sum))
    res = list("S_mat" = S_mat, "X_all" = X_all, "mu" = mu, "sigma" = sigma, "S_test" = S_test, "X_test" = X_test,
               "cods_train_true" = cods_train_true, "cods_test_true" = cods_test_true, "csmf_test_true" = csmf_test_true)
  } else{
    res = list("S_mat" = S_mat, "X_all" = X_all, "mu" = mu, "sigma" = sigma)
  }
  return(res)
}

#num_causes=3; N=list(4,5,2); p=20; k=2; sigmasq=1; sim="F"; lambda_generator=simulate_lambda; 
#N_test=list(0,3,1); bin=T; perturb_scale=1.2; bin_pert_wt=0.8; binary_common=T; common_keep_prop=0.5
simulate_x_simEFG <- function(num_causes, N, p, k, sigmasq, sim, lambda_generator=simulate_lambda, 
                              N_test=NULL, bin_prop=0.8, perturb_scale=1, bin_pert_wt=0.5,
                              binary_common=F, common_keep_prop=0.5, old_meth_bool=T){
  # COVARIATE DEPENDENT MEAN AND/OR COV WITH JUST BINARY INDICATOR!
  # num_causes is the number of causes to simulate
  # N should be a length num_causes list of the number of simulated observations desired for each cause
  # p is number of symptoms, k is underlying dimension of lambda, sigmaaq is noise variance
  # lambda_generator is the function used to make lambda 
  # N_test is a length num_causes list of the number of simulated "test" observations desired for each cause
  # bin is a boolean denoting whether resulting symptom observations should be binary
  # sim ("E","F","G") is a character letter denoting which simulation to run (note H, I, J can use G)
  # bin_pert_wt measures how much weight should go to 'common' mu/cov and how much to binary-specific mu/cov
  # binary_common is a boolean for whether the binary additive term is common (T) or cause-specific (F, default)
  # common_keep_prop the proportion of the mean vector that is common vs. cause-specific, on average
  # old_meth_bool (T/F) specifies whether common mu/sig should be partly weighted by cause as before 3/24 or as new
  library(mvtnorm)
  if( !(length(N)==num_causes) ){ stop('N must be length num_causes') }
  if( !(p>=k) ){ stop('Need p>=k') }
  if( !(sigmasq>0) ){ stop('Need sigmasq>0') }
  if( !(sim%in%c("E","F","G")) ){ stop('Need sim to be one of "E","F","G"') }
  lambda = lambda_generator(k, p)
  perturb = map(1:num_causes, function(x) perturb_scale * lambda_generator(k, p))
  perturb_bin = map(1:num_causes, function(x) perturb_scale * lambda_generator(k, p))
  sigma = mu = list()
  S_mat = S_mat_allbin = S_mat_allcont = X_all = list()
  # Make random set of observations binary
  if( sim=="G" ){
    is_bin = rbernoulli(p, p=bin_prop)
  } else{ is_bin = rep(T, p) } # Make all binary, unless special case.
  if( !is.null(N_test) ){
    if( !(length(N_test)==num_causes) ){ stop('N_test must be length num_causes') }
    N_test_sum = Reduce("+", N_test)
    S_test = S_test_allbin = S_test_allcont = matrix(NA, nrow=N_test_sum, ncol=p) # Will fill in below.
    X_test = matrix(1, nrow=N_test_sum, ncol=2)  # First column here will be 'intercept'.
    X_test[,2] = 1*rbernoulli(N_test_sum,p=0.5) # Make one binary 'indicator' covariate.
    X_test[,1] = X_test[,1] - X_test[,2] # Make one binary 'indicator' covariate.
    test_ind = 1
    cods_test_true = c()
  }
  cods_train_true = c()
  # Make adjustment terms and tuning knobs for covariance adjustment
  common_mu = generate_mean(p);
  common_keep = 1*rbernoulli(p, p=common_keep_prop); # Have ~half of mean elements be cause-dependent, half common
  if( binary_common ){ # Common covariate dependent mean and covariate dependent cov elements
    bin0_mu = generate_mean(p);
    bin1_mu = generate_mean(p);
    lam_baseline = lambda_generator(k, p);
    lam_binary = lambda_generator(k, p);
  }
  mu_baseline = mu_baseline_bin = sigma_baseline = sigma_baseline_bin = mu_bin = list()
  # Create empty list to save true test mean/covariance for each individual
  mu_test = sigma_test = list()
  for( cs in 1:num_causes ){
    mu[[cs]] = list()
    sigma[[cs]] = list()
    S_mat[[cs]] = S_mat_allbin[[cs]] = S_mat_allcont[[cs]] = matrix(NA, nrow=N[[cs]], ncol=p)
    cods_train_true = c(cods_train_true, rep(cs, N[[cs]]))
    # Modify some mean elements to be cause-specific and some to be common
    cause_mu = common_keep*common_mu + (1-common_keep)*generate_mean(p);
    if( !binary_common ){ # Cause-specific covariate dependent mean
      bin0_mu = generate_mean(p);
      bin1_mu = generate_mean(p);
    }
    if(sim=="F"){ # Then common mean for everyone!
      mu_baseline[[cs]] = mu_bin[[cs]] = common_mu; 
    } else{ # Baseline comes partly from overall common mean, and at cause level split by binary
      if( old_meth_bool ){
        mu_baseline[[cs]] = (1-bin_pert_wt)*cause_mu + bin_pert_wt*bin0_mu
        mu_bin[[cs]] = (1-bin_pert_wt)*cause_mu - bin_pert_wt*bin1_mu
      } else{
        mu_baseline[[cs]] = (1-bin_pert_wt)*common_mu + bin_pert_wt*bin0_mu
        mu_bin[[cs]] = (1-bin_pert_wt)*common_mu - bin_pert_wt*bin1_mu
      }
    }
    if( sim=="E" ){ # Keep the same baseline when binary cov == 1
      sigma_baseline_bin[[cs]] = sigma_baseline[[cs]] = diag(rep(1,p)) # for E, shared independent
    } else{ # Make a different "baseline" for when binary cov == 1
      if( (!binary_common) | old_meth_bool ){
        sigma_baseline[[cs]] = cov2cor(sigmasq*diag(p) + (1-bin_pert_wt)*tcrossprod(lambda) + 
                                       bin_pert_wt*tcrossprod(perturb[[cs]]))
        sigma_baseline_bin[[cs]] = cov2cor(sigmasq*diag(p) + (1-bin_pert_wt)*tcrossprod(lambda) + 
                                           bin_pert_wt*tcrossprod(perturb_bin[[cs]]))
      } else{
        cause_sig = common_keep_prop*tcrossprod(lambda) + (1-common_keep_prop)*tcrossprod(perturb[[cs]])
        sigma_baseline[[cs]] = cov2cor(sigmasq*diag(p) + (1-bin_pert_wt)*cause_sig + 
                                         bin_pert_wt*tcrossprod(lam_baseline))
        sigma_baseline_bin[[cs]] = cov2cor(sigmasq*diag(p) + (1-bin_pert_wt)*cause_sig + 
                                             bin_pert_wt*tcrossprod(lam_binary))
      }
    }
    # Create covariate matrix
    X_all[[cs]] = matrix(1, nrow=N[[cs]], ncol=2)
    X_all[[cs]][,2] = 1*rbernoulli(N[[cs]],p=0.5) # Make one binary 'indicator' covariate.
    X_all[[cs]][,1] = X_all[[cs]][,1] - X_all[[cs]][,2]
    # Sample from mvnormal for each person
    for( i in 1:N[[cs]] ){
      bin_cov_tmp = X_all[[cs]][i,2] # 1 when person i has binary ind, 0 else
      if(bin_cov_tmp==1){mu[[cs]][[i]] = mu_bin[[cs]]}
      if(bin_cov_tmp==0){mu[[cs]][[i]] = mu_baseline[[cs]]}
      if( !(sim=="E") ){
        sigma[[cs]][[i]] = (1-bin_cov_tmp)*sigma_baseline[[cs]] + bin_cov_tmp*sigma_baseline_bin[[cs]]
      } else{
        sigma[[cs]][[i]] = sigma_baseline[[cs]]
      }
      S_mat_allcont[[cs]][i,] = rmvnorm(n=1, mean = mu[[cs]][[i]], sigma = sigma[[cs]][[i]])
      if(F){
        # For plotting (visualization):
        image(sigma_baseline[[cs]])
        image(sigma_baseline_bin[[cs]])
        image(sigma[[cs]][[i]])
      }
    }
    for( j in 1:p ){
      if( is_bin[j] ){
        S_mat[[cs]][,j] = 1*(S_mat_allcont[[cs]][,j]>0)
      } else{
        S_mat[[cs]][,j] = S_mat_allcont[[cs]][,j]
      }
      S_mat_allbin[[cs]][,j] = 1*(S_mat_allcont[[cs]][,j]>0)
    }
    
    if( !is.null(N_test) ){
      if( !(N_test[[cs]]==0) ){
        for(i in test_ind:(test_ind+N_test[[cs]]-1)){
          bin_cov_tmp = X_test[i,2] # 1 when person i has binary ind, 0 else
          if(bin_cov_tmp==1){mu_test[[i]] = mu_bin[[cs]]}
          if(bin_cov_tmp==0){mu_test[[i]] = mu_baseline[[cs]]}
          if( !(sim=="E") ){
            sigma_test[[i]] = (1-bin_cov_tmp)*sigma_baseline[[cs]] + bin_cov_tmp*sigma_baseline_bin[[cs]]
          } else{
            sigma_test[[i]] = sigma_baseline[[cs]]
          }
          S_test_allcont[i,] = rmvnorm(n=1, mean = mu_test[[i]], sigma = sigma_test[[i]])
        }
        
      }
      test_ind = test_ind+N_test[[cs]]
      cods_test_true = c(cods_test_true, rep(cs, N_test[[cs]]))
      
    }
  } # for( cs in 1:num_causes )
  if( !is.null(N_test) ){
    for( j in 1:p ){
      if( is_bin[j] ){
        S_test[,j] = 1*(S_test_allcont[,j]>0)
      } else{
        S_test[,j] = S_test_allcont[,j]
      }
      S_test_allbin[,j] = 1*(S_test_allcont[,j]>0)
    }
    csmf_test_true = unlist(lapply(N_test,function(l) l/N_test_sum))
    res = list("S_mat" = S_mat, "S_mat_allbin" = S_mat_allbin, "S_mat_allcont" = S_mat_allcont,
               "X_all" = X_all, "mu_baseline" = mu_baseline, "mu_bin" = mu_bin, 
               "sigma_baseline" = sigma_baseline, "sigma_baseline_bin" = sigma_baseline_bin,
               "S_test" = S_test, "S_test_allbin" = S_test_allbin, "S_test_allcont" = S_test_allcont, "X_test" = X_test,
               "cods_train_true" = cods_train_true, "cods_test_true" = cods_test_true, "csmf_test_true" = csmf_test_true,
               "sigma" = sigma, "sigma_test" = sigma_test)
  } else{
    res = list("S_mat" = S_mat, "S_mat_allbin" = S_mat_allbin, 
               "X_all" = X_all, "mu_baseline" = mu_baseline, "mu_bin" = mu_bin, 
               "sigma_baseline" = sigma_baseline, "sigma_baseline_bin" = sigma_baseline_bin,
               "cods_train_true" = cods_train_true)
  }
  return(res)
}
# plot(dat$sigma_baseline[[1]], dat$sigma_baseline_bin[[1]]); abline(0,1)
# plot(dat$sigma_baseline[[1]], dat$sigma_baseline[[2]]); abline(0,1)

