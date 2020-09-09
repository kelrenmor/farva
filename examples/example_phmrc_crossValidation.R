### Demonstration of how to perform cross validation for model selection
### using the farva package on the PHMRC data set.
### For illustrative purposes, deaths from Mexico having observed 
### age data are extracted and the model is run:
### (1) without age as a covariage, and
### (2) with age impacting the symptom mean.

################# Load required libraries and data #################

rm(list=ls())
set.seed(1204)

# Install and load farva library, load other required libraries
remotes::install_github('kelrenmor/farva', quiet=T, dependencies=T)
library(farva)
library(openVA)

# Read in PHMRC data using openVA::getPHMRC_url() function
PHMRC_adult <- read.csv(openVA::getPHMRC_url("adult"))
phmrc_convert <- ConvertData.phmrc(PHMRC_adult, phmrc.type = "adult", cutoff = "default", cause = "va34")[[1]]

################# Extract certain columns as vectors and make a test/train split of the data #################

# Get site, cod, age as vectors
site <- PHMRC_adult$site
cods <- PHMRC_adult$gs_text34; num_causes = length(unique(cods))
age <- as.numeric(as.character(PHMRC_adult$g1_07a))

# Make test/train indices for datasets, setting random number seed for reproducibility
# (Note here I'm only keeping deaths from Mexico and discarding any obs where age is NA)
inds_keep <- (site=="Mexico") & (!is.na(age))
phmrc_red <- phmrc_convert[inds_keep,-1] # Keep only obs from Mexico with age observed
age_red <- age[inds_keep]

# Set k_fold, for k-fold cross validation
k_fold <- 5
# Get numbers of training/test data for k-fold cross validation
train_num <- round((k_fold-1)*nrow(phmrc_red)/k_fold)
test_num <- nrow(phmrc_red) - train_num
# Get fold id (random sorted) for each person in data set
fold_id <- c()
for(kk in 1:k_fold){
  if(kk<k_fold){
    fold_id <- c( fold_id, rep(kk,test_num) )
  } else{
    fold_id <- c( fold_id, rep(kk, nrow(phmrc_red) - test_num*(k_fold-1)) )
  }
}
fold_id <- sample(fold_id)

################# Run the model for each fold, save results #################

acc_1_all = acc_csmf_all = rep(NA, k_fold)
cacc_1_age_all = acc_csmf_age_all = rep(NA, k_fold)
for(kk in 1:k_fold){
  
  # Sample training indices, use rest as test data
  train_ind <- which(fold_id==kk)
  test_ind <- setdiff(1:nrow(phmrc_red), train_ind)
  
  # Split data into testing and training data set and convert to numeric (0/1 vs Y/N)
  phmrc_train <- as.matrix( farva::conv_to_numeric( phmrc_red[train_ind,] ) )
  phmrc_test <- as.matrix( farva::conv_to_numeric( phmrc_red[test_ind,] ) )
  age_train <- age_red[train_ind]; age_test <- age_red[test_ind]
  
  #### Create the covariate matrices and run the model ####
  
  # Get X matrix for train/test data (i.e., covariates for each data set)
  N_train <- nrow(phmrc_train); N_test <- nrow(phmrc_test)
  X_train <- matrix(1,nrow=N_train,ncol=2); X_test <- matrix(1,nrow=N_test,ncol=2)
  X_train[,1] <- 1*(age_train>=65); X_train[,2] <- 1-X_train[,1]
  X_test[,1] <- 1*(age_test>=65); X_test[,2] <- 1-X_test[,1]
  
  # Run FARVA model without age!
  # Note that because X_.._mu and X_.._sig are undefined, the symptom-level covariance will not depend on age.
  farva_res = farva_run(S_mat=phmrc_train, S_test=phmrc_test, 
                        save_num_tot=500, thin=10, burnin=5000, K=10, verbose=F)
  
  # Run FARVA model with age!
  # Note that because X_.._sig is set to NULL, the symptom-level covariance will not depend on age.
  farva_res_age = farva_run(S_mat=phmrc_train, X_all_mu=X_train, X_all_sig=NULL, 
                            S_test=phmrc_test, X_test_mu=X_test, X_test_sig=NULL, 
                            save_num_tot=500, thin=10, burnin=5000, K=10, verbose=F)
  
  #### Get out matrix of individual probabilities and CSMFs ####
  
  for(ii in 1:2){
    # Define data depending on with/without age
    if(ii==1){ # without age
      cod_probs_test = farva_res$indiv_prob # N_test by number of training causes dimensional
      csmf_test = apply(farva_res$csmf_test_save,1,mean) # number of training causes length
    } else{ # with age
      cod_probs_test = farva_res_age$indiv_prob # N_test by number of training causes dimensional
      csmf_test = apply(farva_res_age$csmf_test_save,1,mean) # number of training causes length
    }
    
    # Get matrix of probabilities for each individual in the test set to have each training COD
    # 0-pad to get probs of test individual to have each POSSIBLE COD from original set of all CODs
    # (i.e., just put 0 probability for CODs that are unseen in the training data)
    cod_probs = matrix(0, nrow=N_test, ncol=num_causes) # N_test by num_causes matrix
    colnames(cod_probs) = 1:num_causes
    inds_match = match(colnames(cod_probs_test), colnames(cod_probs))
    cod_probs[,inds_match] = cod_probs_test
    # Get top cause accuracy
    acc1_tmp = acc_num(indiv_prob=cod_probs, true_causes=phmrc_test[,1],top_num=1)
    
    # Get model-predicted CSMF for the TEST data
    # As above, can 0-pad to get the CSMF where you assume no test death types that weren't in training set
    csmf_test_all = rep(0, num_causes) # num_causes length
    csmf_test_all[inds_match] = csmf_test
    # Get CSMF accuracy
    tmp_causes = factor(phmrc_test[,'Cause'], levels=1:num_causes)
    acc_csmf_tmp = csmf_accuracy(csmf_true=table(tmp_causes)/nrow(phmrc_test), csmf_compare=csmf_test_all)
    
    # Save results
    if(ii==1){ # without age
      acc_1_all[kk] = acc1_tmp
      acc_csmf_all[kk] = acc_csmf_tmp
    } else{ # with age
      acc_1_age_all[kk] = acc1_tmp
      acc_csmf_age_all[kk] = acc_csmf_tmp
    }
  }

}

################# (View performance of models) #################

acc_1_all
acc_csmf_all

acc_1_age_all
acc_csmf_age_all
