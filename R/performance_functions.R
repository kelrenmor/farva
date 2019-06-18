whichpart <- function(x, n) {
  nx <- length(x)
  p <- nx-n
  xp <- as.numeric(sort(x, partial=p)[p])
  which(x > xp)
}

acc_num = function(indiv_prob, true_causes, top_num){
  # indiv_prob should be N x num_causes matrix.
  # true_causes should be N-length vector.
  # top_num is how many of top probs to include (=1 for ACC_1, =3 for ACC_3)
  # result will be ACC_{top_num}
  top_match = rep(NA,nrow(indiv_prob))
  for(ii in 1:nrow(indiv_prob)){
    top_match[ii] = (true_causes[ii] %in% whichpart(indiv_prob[ii,], n=top_num))
  }
  return(mean(top_match))
}

ch_conc = function(indiv_prob, true_causes, wts="equal"){
  # indiv_prob should be N x num_causes matrix.
  # true_causes should be N-length vector.
  # wts should be a num_causes length vector of weights to be used in calc
  # result will be the overall CCC (from Murray 2011 Robust Metrics paper)
  top_match = apply(indiv_prob,1,which.max)
  causes = sort(unique(true_causes))
  num_causes = length(causes)
  ccc_j = rep(NA, num_causes)
  for(cc in 1:num_causes){
    cause = causes[[cc]]
    tmp_vec = top_match[true_causes==cause]
    ccc_j[cc] = (sum(tmp_vec==cause)/(length(tmp_vec)) - 1/num_causes)/(1-1/num_causes)
  }
  if(wts=="equal"){
    return(mean(ccc_j))
  } else(
    return(mean(wts*ccc_j))
  )
}

csmf_accuracy <- function(csmf_true, csmf_compare){
  num_causes <- length(csmf_true)
  return(1 - sum(abs(csmf_true-csmf_compare))/2*(1-min(csmf_true)))
}