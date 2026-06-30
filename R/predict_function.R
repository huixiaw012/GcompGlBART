#' @exportS3Method predict gc_gsoftbart_regression

predict.gc_gsoftbart_regression <- function(object, newdata_X, newdata_Z, iterations = NULL, ...) {


  opts <- object$opts

  Z <- as.matrix(newdata_Z)
  X <- as.matrix(newdata_X)

  if(opts$num_save == 1){
    beta <- object$beta[1:opts$num_save,]#/BM$sd_Y
  } else {
    beta <- colMeans(object$beta[1:opts$num_save,])#/BM$sd_Y
  }

    if(is.null(iterations))
    iterations <- seq(from = opts$num_burn+opts$num_thin,
                      to = opts$num_burn + opts$num_thin * opts$num_save,
                      by = opts$num_thin)


  for(i in 1:ncol(X)) {
    X[,i] <- object$ecdfs[[i]](X[,i])
  }

  pi <- function(i) {
    r <- as.numeric(object$forest$predict_iteration(X, i)) * object$sd_Y + object$mu_Y
    #  forest$predict_iteration(X, i) returns the predictions from a
    #   matrix X of predictors at iteration i. Requires that opts$cache_trees =
    #   TRUE in MakeForest(hypers, opts).
    eta <- as.numeric(Z %*% beta)

    return(r + eta)
  }

  mu <- t(sapply(iterations, pi))
  mu_hat <- colMeans(mu)

  return(list(mu = mu, mu_mean = mu_hat))
}

#' @exportS3Method predict gc_gsoftbart_probit
predict.gc_gsoftbart_probit <- function(object, pre_object = NULL, newdata_X, newdata_Z, iterations = NULL, ...) {


  opts <- object$opts

  Z <- as.matrix(newdata_Z)
  X <- as.matrix(newdata_X)

  if(is.null(pre_object)){
    if(opts$num_save == 1){
      beta <- object$beta[1:opts$num_save,]#/BM$sd_Y
    } else {
      beta <- colMeans(object$beta[1:opts$num_save,])#/BM$sd_Y
    }
  }else{

    ln.form <- object$linear_formula
    ln.vars <- all.vars(ln.form) # Extract variable names from linear formula

    pre.ln.form <- pre_object$linear_formula
    pre.ln.vars <- all.vars(pre.ln.form)

    beta <- object$beta#[1:opts$num_save,]
    pre.beta <- pre_object$beta#[1:opts$num_save,]

    colnames(beta) <- ln.vars
    colnames(pre.beta) <- pre.ln.vars

    not_in_formula <- setdiff(pre.ln.vars, ln.vars) # Find names in the previous beta that are NOT in the current beta
    common_vars <- intersect(ln.vars, pre.ln.vars)

    new.beta <- matrix(NA,nrow = nrow(pre.beta), ncol = ncol(pre.beta) )
    colnames(new.beta) <- pre.ln.vars
    new.beta[,common_vars] <- beta

    for (rth in 1:opts$num_save) {
      dat.beta  <- data.frame(pre.beta = pre.beta[rth,], beta = new.beta[rth,])
      mod.beta <- lm(beta~pre.beta, dat.beta[complete.cases(dat.beta),])
      dat.beta[not_in_formula,]$beta <- predict(mod.beta,dat.beta)[not_in_formula]
      new.beta[rth,] <- dat.beta$beta
    }
    if(opts$num_save == 1){
      beta <- new.beta[1:opts$num_save,]#/BM$sd_Y
    } else {
      beta <- colMeans(new.beta[1:opts$num_save,])#/BM$sd_Y
    }
  }



  if(is.null(iterations))
    iterations <- seq(from = opts$num_burn+opts$num_thin,
                      to = opts$num_burn + opts$num_thin * opts$num_save,
                      by = opts$num_thin)


  for(i in 1:ncol(X)) {
    X[,i] <- object$ecdfs[[i]](X[,i])
  }

  pi <- function(i) {
    r <- as.numeric(object$forest$predict_iteration(as.matrix(X), i)) + object$offset
    eta <- as.numeric(Z %*% beta) #* BM$sd_Y

    return(r + eta)
  }

  mu <- t(sapply(iterations, pi))
  p <- pnorm(mu)

  mu_hat <- colMeans(mu)
  p_hat <- colMeans(p)

  return(list(mu = mu, p = p, mu_mean = mu_hat, p_mean = p_hat))
}

