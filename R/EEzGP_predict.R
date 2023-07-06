#' @title The Prediction Function of \code{EEzGP} Model
#'
#' @description Predicts the output of the EEzGP model fitted by \code{\link[EzGP]{EEzGP_fit}}.
#'
#' @param X_new Matrix or vector containing the input(s) where the predictions are to be made. Each row is an input vector.
#' @param model The EEzGP model fitted by \code{\link[EzGP]{EEzGP_fit}}.
#' @param MSE_on A scalar indicating whether the uncertainty (i.e., mean squared error \code{MSE}) is calculated. Set to a non-zero value to calculate \code{MSE}.
#'
#' @import methods
#'
#' @return A prediction list containing the following components:
#' \itemize{
#' \item{\code{Y_hat}} {A vector containing the prediction values}
#' \item{\code{MSE}} {A vector containing the prediction uncertainty (i.e., the covariance or covariance matrix for the output(s) at prediction location(s)) }
#' }
#'
#' @export
#'
#' @references
#' \enumerate{
#' \item "EzGP: Easy-to-Interpret Gaussian Process Models for Computer Experiments with Both Quantitative and Qualitative Factors", Qian Xiao, Abhyuday Mandal, C. Devon Lin, and Xinwei Deng (\doi{10.1137/19M1288462})
#' }
#'
#' @seealso
#' \code{\link[EzGP]{EEzGP_fit}} to fit EEzGP model for the datasets.\cr
#'
#' @examples
#' # This function is used in a similar way as the use of EzGP_predict.
#' # See the examples in the documentation of the function EEzGP_fit.

EEzGP_predict <- function(X_new, model, MSE_on = 0){

  ## import model and check inputs
  if (!is(model,"EEzGP model")){
    stop('    The 2nd input should be a model of class "EEzGP model".')
  }
  if (length(MSE_on)!=1){
    stop('    MSE_on should be a scalar flag. Set it to 1 to turn it "on".')
  }

  X <- model$data$X
  Y <- model$data$Y
  p <- model$data$p
  q <- model$data$q
  m <- model$data$m
  n <- nrow(model$data$X)
  tau <- model$data$tau
  parv <- model$param
  covm <- cov_m(X, p, q, m, n, parv, tau, models = 1)

  y = as.matrix(Y)
  Tm = try(chol(covm), silent=TRUE)
  #round(t(Tm)%*%Tm,2) == round(covm,2)
  if ('try-error' %in% class(Tm)) {
    return(NULL)
  }
  m1 = as.matrix(c(rep(1,n)))
  invT = backsolve(Tm, diag(dim(Tm)[1]))
  invc = invT%*%t(invT)
  mu = as.numeric(1/sum(invc)*(t(m1) %*% invc %*% y))


  ##total number of parameters in the model
  npar = 1 + q + p + (sum(m)-q)

  ## a help function
  psum <- function(x1,x2, par2)
  {
    return(sum(-par2*(x1-x2)^2))
  }

  # calculating covariance between two inputs w1 and w2
  covx <- function(w1,w2, parv){
    #variance parameter sigma^2
    par1 = parv[1:(q+1)]
    #correlation parameter in G0
    par2 = parv[(q+2):(q+1+p)]
    #correlation parameter in G1 to Gq
    par3 = parv[(q+2+p): npar]
    x1 = w1[1:p]
    z1 = w1[(p+1):(p+q)]
    x2 = w2[1:p]
    z2 = w2[(p+1):(p+q)]
    res1 = par1[1]*exp(psum(x1,x2,par2))
    for (i in 1:q){
      if(z1[i] != z2[i]){
        res1 = res1+0
      }
      else{
        l = z1[i]
        if (l==1)
        {
          res1 = res1 + par1[i+1]*exp(psum(x1,x2,1))
        }
        else {
          res1 = res1 + par1[i+1]*exp(psum(x1,x2,par3[sum((m-1)[1:i]) - (m-1)[i] + l-1]))
        }
      }
    }
    return(res1)
  }

  prey <- function(wn){
    covv = matrix(0,n)
    for(i in 1:n){
      if (sum(round(wn,5)!=round(X[i,],5)) > 0)
      {
        covv[i] = covx(wn,X[i,],parv)
      } else {
        covv[i] = sum(parv[1:(q+1)]) + tau
      }
    }
    gamma = as.matrix(covv)

    # Y_hat
    Y_hat = mu  + t(gamma) %*% invc %*% (y - mu * m1)
    # MSE
    MSE = sum(parv[1:(q+1)]) - (t(gamma) %*% invc %*% gamma)
    + ((1 - t(m1) %*% invc %*% gamma)^2)/(t(m1) %*% invc %*% m1)
    result <- c(Y_hat, MSE)

    return(result)
  }

  prediction <- list()
  if (is.vector(X_new) == TRUE){

    p_all <- ncol(model$data$X)
    if (length(X_new) != p_all){
      stop('    The dimensionality of X_new is not correct!')
    }

    value <- prey(X_new)
    if (MSE_on){
      prediction$Y_hat <- value[1]
      prediction$MSE <- value[2]
    } else{
      prediction$Y_hat <- value[1]
    }

  } else if (is.matrix(X_new) == TRUE){

    p_all <- ncol(model$data$X)
    if (ncol(X_new) != p_all){
      stop('    The dimensionality of X_new is not correct!')
    }

    Y_hat = c()
    MSE = c()
    nn = nrow(X_new)
    for (i in 1:nn){
      value = prey(X_new[i,])
      Y_hat[i] = value[1]
      MSE[i] = value[2]
    }

    if (MSE_on){
      prediction$Y_hat <- Y_hat
      prediction$MSE <- MSE
    } else{
      prediction$Y_hat <- Y_hat
    }

  } else{
    stop('    X_new must be a matrix or a vector')
  }

  return(prediction)
}
