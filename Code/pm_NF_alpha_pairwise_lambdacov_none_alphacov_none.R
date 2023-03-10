#' Ricker model with pairwise alphas and no covariate effects
#'
#' @param par 1d vector of initial parameters: 'lambda', 'alpha_intra' (optional), 'alpha_inter', and 'sigma'
#' @param fitness 1d vector of fitness observations, in log scale
#' @param neigh_intra_matrix optional matrix of one column, number of intraspecific neighbours for each observation
#' @param neigh_inter_matrix matrix of arbitrary columns, number of interspecific neighbours for each observation
#' @param covariates included for compatibility, not used in this model
#' @param fixed_parameters optional list specifying values of fixed parameters, 
#' with components "lambda","alpha_intra","alpha_inter".
#'
#' @return log-likelihood value
#' @export
pm_NF_alpha_pairwise_lambdacov_none_alphacov_none <- function(par,
                                                              fitness,
                                                              neigh_intra_matrix = NULL,
                                                              neigh_inter_matrix,
                                                              covariates,
                                                              fixed_parameters){
  
  
  # quickest way to include germination rates is to hard-code them here
  germ <- c(BEMA = 0.38, CETE = 0.53, CHFU = 0.8, CHMI = 0.76, HOMA = 0.94,
            LEMA = 0.89, LYTR = 0.62, MEEL = 0.59, MESU = 0.77, PAIN = 0.33,
            PLCO = 0.84, POMA = 0.85, POMO = 0.91, PUPA = 0.84, SASO = 0.52,
            SCLA = 0.69, SOAS = 0.6, SPRU = 0.44, SUSP = 0.61)
  
  # retrieve parameters -----------------------------------------------------
  # parameters to fit are all in the "par" vector,
  # so we need to retrieve them one by one
  # order is {lambda,lambda_cov,alpha_intra,alpha_inter,alpha_cov,sigma}
  
  # comment or uncomment sections for the different parameters
  # depending on whether your model includes them
  pos <- 1
  
  # if a parameter is passed within the "par" vector,
  # it should be NULL in the "fixed_parameters" list
  if(is.null(fixed_parameters[["lambda"]])){
    lambda <- par[pos]
    pos <- pos + 1
  }else{
    lambda <- fixed_parameters[["lambda"]]
  }
  
  # if(is.null(fixed_parameters$lambda_cov)){
  #   lambda_cov <- par[pos:(pos+ncol(covariates)-1)]
  #   pos <- pos + ncol(covariates)
  # }else{
  #   lambda_cov <- fixed_parameters[["lambda_cov"]]
  # }
  
  if(!is.null(neigh_intra_matrix)){
    # intra
    if(is.null(fixed_parameters[["alpha_intra"]])){
      alpha_intra <- par[pos]
      pos <- pos + 1
    }else{
      alpha_intra <- fixed_parameters[["alpha_intra"]]
    }
  }else{
    alpha_intra <- NULL
  }
  
  # inter
  if(is.null(fixed_parameters[["alpha_inter"]])){
    alpha_inter <- par[pos:(pos+ncol(neigh_inter_matrix)-1)]
    pos <- pos + ncol(neigh_inter_matrix) -1
  }else{
    alpha_inter <- fixed_parameters[["alpha_inter"]]
  }
  
  # if(is.null(fixed_parameters$alpha_cov)){
  #   alpha.cov <- par[pos:(pos+(ncol(covariates)*ncol(neigh_matrix))-1)]
  #   pos <- pos + (ncol(covariates)*ncol(neigh_matrix))
  # }else{
  #   alpha.cov <- fixed_parameters[["alpha.cov"]]
  # }
  
  sigma <- par[length(par)]
  
  # now, parameters have appropriate values (or NULL)
  # next section is where your model is coded
  
  # MODEL CODE HERE ---------------------------------------------------------
  
  # we do not differentiate alpha_intra from alpha_inter in this model
  # so, put together alpha_intra and alpha_inter, and the observations
  # with intraspecific ones at the beginning
  if(!is.null(alpha_intra)){
    alpha <- c(alpha_intra,alpha_inter)
    all_neigh_matrix <- cbind(neigh_intra_matrix,neigh_inter_matrix)
  }else{
    alpha <- alpha_inter
    all_neigh_matrix <- neigh_inter_matrix
  }
  
  # check germination rate against par names
  my.g <- germ[sort(names(alpha))]
  
  # sort all parameters according to their names
  alpha <- alpha[names(my.g)]
  all_neigh_matrix <- all_neigh_matrix[,names(my.g)]
  
  # alpha <- alpha * -1
  
  term = 0 #create the denominator term for the model
  for(z in 1:ncol(all_neigh_matrix)){
    term <- term - alpha[z]*log((my.g[z]*all_neigh_matrix[,z]) + 1) 
  }
  pred <- lambda * exp(term)
  
  # MODEL CODE ENDS HERE ----------------------------------------------------
  
  # the routine returns the sum of log-likelihoods of the data and model:
  # DO NOT CHANGE THIS
  llik <- dnorm(fitness, mean = (log(pred)), sd = (sigma), log=TRUE)
  return(sum(-1*llik))
}