#' @title The Log-likelihood Function and The Analytical Gradients in \code{EzGP} Package
#'
#' @description Calculates the log-likelihood function value and the analytical gradients as described in \code{reference 1}.
#'
#' @param X Matrix or data frame containing the inputs of training data. Each row represents the input setting of a data point and the columns are values of quantitative variables and qualitative variables.
#' @param Y Vector containing the outputs of training data points.
#' @param p Number of quantitative factors in the given dataset \code{X}.
#' @param q Number of qualitative factors in the given dataset \code{X}.
#' @param m A vector containing numbers of levels in qualitative factors.
#' @param parv Parameters in the EzGP/EEzGP model.
#' @param tau Nugget if needed. The default nugget is 0, otherwise it has to be a non-negative real value.
#' @param models Model indicator that indicates which model the likelihoods and analytical gradients are applied to. 0 for EzGP model, 1 for EEzGP model.

#' @return A list of the following items:
#' \itemize{
#' \item{\code{objective}} {The log-likelihood function value.}
#' \item{\code{gradient}} {The analytical gradients.}
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
#' \code{\link[EzGP]{EzGP_fit}} to see how an EzGP model can be fitted to a training dataset.\cr
#' \code{\link[EzGP]{EzGP_predict}} to use the fitted EzGP model for prediction.\cr
#' \code{\link[EzGP]{EEzGP_fit}} to see how an EEzGP model can be fitted to a training dataset.\cr
#' \code{\link[EzGP]{EEzGP_predict}} to use the fitted EEzGP model for prediction.\cr
#' \code{\link[EzGP]{LEzGP_fit}} to see how a LEzGP model can be fitted to a training dataset.\cr

#' @examples
#' # see the examples in the documentation of the function EzGP_fit.

LLF_gradients <- function(X, Y, p, q, m, parv, tau = 0, models = 0) {

  # function input checking
  if (missing(X)){
    stop('    X must be provided.')
  }
  if (!is.matrix(X) && !is.data.frame(X)){
    stop('    X must be a matrix or a data frame.')
  }

  # Number of training data
  n = nrow(X)

  if (missing(Y)){
    stop('    Y must be provided.')
  }
  if (!all(is.finite(Y)) || !is.numeric(Y)){
    stop('    All the elements of Y must be finite numbers.')
  }
  if (is.vector(Y) == TRUE) {
    y <- as.matrix(Y)
  } else if (is.matrix(Y) == TRUE){
    if (ncol(Y) == 1){
      y <- Y
    } else{
      stop('    The response value (i.e. y) for each observation should has only one.')
    }
  }

  if (n != nrow(y)){
    stop('    The number of rows (i.e., observations) in X and Y should match!')
  }

  if (!is.finite(p) || !is.numeric(p)){
    stop('    p must be a finite number')
  }

  if (!is.finite(q) || !is.numeric(q)){
    stop('    q must be a finite number')
  }

  if (ncol(X) != (p+q)){
    stop('    The number of quantitative and qualitative factors does not match the input training data')
  }

  if (is.vector(m) == FALSE){
    stop('    m must be a vetcor')
  }

  if (length(m) != q){
    stop('    m must contain the number of levels of each qualitative factors, i.e. m should be a vector with length q')
  }

  ##total number of parameters in the model
  if (models == 0){
    npar = 1 + q + p + p*sum(m)
  } else if (models == 1){
    npar = 1 + q + p + (sum(m)-q)
  }

  ## a help function
  psum <- function(x1,x2, par2)
  {
    return(sum(-par2*(x1-x2)^2))
  }

  #check if all parameters are all positive
  if (min(parv) < 0)
  {
    return(list( "objective" = 10000000000,
                 "gradient" = rep(NA, npar)))
  }
  #expand grid to avoid for looping
  rcoord <- cbind(
    rep(seq_len(n - 1), times = rev(seq_len(n - 1))),
    unlist(lapply(
      X = rev(seq_len(n - 1)),
      FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = n)))

  #covariance matrix
  covm = cov_m (X, p, q, m, n, parv, tau, models)

  Tm = try(chol(covm), silent=TRUE)
  #round(t(Tm)%*%Tm,2) == round(covm,2)
  if ('try-error' %in% class(Tm)) {
    return(list( "objective" = 10000000000,
                 "gradient" = rep(NA, npar)))
  }
  m1 = as.matrix(c(rep(1,n)))
  invT = backsolve(Tm, diag(dim(Tm)[1]))
  invc = invT%*%t(invT)

  #likelihood function value
  MLE_result = 2*sum(log(diag(Tm))) + (t(y) %*% invc %*% y) - 1/sum(invc)*(t(m1) %*% invc %*% y)^2



  #calcaulte analytical gradients
  # trend estimator
  mu = as.numeric(1/sum(invc)*(t(m1) %*% invc %*% y))

  ### derivative for variance parameter sigma^2_0
  grad_var0 <- function()
  {
    ### derative function of sigma_0 for wi, wj
    gradf_var <- function(w12)
    {
      x1 = w12[1:p]
      x2 = w12[(p+q+1):(2*p+q)]
      #correlation parameter in G0
      par2 = parv[(q+2):(q+1+p)]
      gx = exp(psum(x1,x2,par2))
      return(as.numeric(gx))
    }
    der_m = matrix(0,n,n)
    Rtemp_var1 <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_var)
    der_m[rcoord] <- Rtemp_var1
    der_m <- der_m + t(der_m)
    diag(der_m) = 1
    result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
    return(result)
  }

  #grad_var0()

  ### derivative for variance parameter sigma^2_h h =1, ..., q
  grad_var <- function(h)
  {
    ### derative function of sigma_h for wi, wj
    gradf_var <- function(w12,h)
    {
      x1 = w12[1:p]
      z1 = w12[(p+1):(p+q)]
      x2 = w12[(p+q+1):(p+q+p)]
      z2 = w12[(p+q+p+1):(2*p+2*q)]

      if (models == 0){
        if(z1[h] != z2[h]){
          return(0)
        }
        else{
          l = as.numeric(z1[h])
          par3 = parv[(q+2+p): npar]
          gx = exp(psum(x1,x2, par3[(sum(m[1:h])*p - m[h]*p + (l-1)*p + 1) : (sum(m[1:h])*p - m[h]*p + (l-1)*p + p)]))
          return(as.numeric(gx))
        }
      } else if (models == 1){
        if(z1[h] != z2[h]){
          return(0)
        }
        else{
          l = as.numeric(z1[h])
          if (l==1)
          {
            gx = exp(psum(x1,x2,1))
          } else {
            par3 = parv[(q+2+p): npar]
            gx = exp(psum(x1,x2,par3[sum((m-1)[1:h]) - (m-1)[h] + l-1]))
          }
          return(as.numeric(gx))
        }
      }
    }
    der_m = matrix(0,n,n)
    Rtemp_var1 <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_var, h=h)
    der_m[rcoord] <- Rtemp_var1
    der_m <- der_m + t(der_m)
    diag(der_m) = 1
    result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
    return(result)
  }

  ### derivative for correlation parameters theta_(0)_s for s=1,...,p
  grad_cor0 <- function(s)
  {
    ### derative function of theta_0_s for wi, wj
    gradf_cor0 <- function(w12,s)
    {
      x1 = w12[1:p]
      #z1 = w12[(p+1):(p+q)]
      x2 = w12[(p+q+1):(p+q+p)]
      #z2 = w12[(p+q+p+1):(2*p+2*q)]
      #correlation parameter in G0
      par2 = parv[(q+2):(q+1+p)]
      gx = -parv[1] * (x1[s] - x2[s])^2 * exp(psum(x1,x2,par2))
      return(as.numeric(gx))
    }
    der_m = matrix(0,n,n)
    Rtemp_cor <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_cor0, s=s)
    der_m[rcoord] <- Rtemp_cor
    der_m <- der_m + t(der_m)
    diag(der_m) = 0
    result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
    return(result)
  }
  if (models == 0){
    ### derivative for correlation parameters theta_(h)_l_s for h=1,...,q, l=1,...,m_h, s=1,...,p
    grad_cor <- function(h,l,s)
    {
      ### derative function of theta_0_s for wi, wj
      gradf_corhs <- function(w12,h,l,s)
      {
        x1 = w12[1:p]
        z1 = w12[(p+1):(p+q)]
        x2 = w12[(p+q+1):(p+q+p)]
        z2 = w12[(p+q+p+1):(2*p+2*q)]
        if((z1[h] != l) | (z2[h] != l)){
          return(0)
        }
        else{
          #variance parameter sigma^2
          par1 = parv[1:(q+1)]
          #correlation parameter in G1 to Gq
          par3 = parv[(q+2+p): npar]
          gx = -par1[h+1] * (x1[s] - x2[s])^2 * exp(psum(x1,x2, par3[(sum(m[1:h])*p - m[h]*p + (l-1)*p + 1) : (sum(m[1:h])*p - m[h]*p + (l-1)*p + p)]))
          return(as.numeric(gx))
        }
      }
      der_m = matrix(0,n,n)
      Rtemp_cor <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_corhs, h=h, l=l, s=s)
      der_m[rcoord] <- Rtemp_cor
      der_m <- der_m + t(der_m)
      diag(der_m) = 0
      result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
      return(result)
    }
    ###analytical gradients results
    grad_res = grad_var0()
    for (i in 1:q) {
      grad_res = c(grad_res, grad_var(i))
    }
    for (i in 1:p)
    {
      grad_res = c(grad_res, grad_cor0(i))
    }
    for (h0 in 1:q)
    {
      for (l0 in 1:m[h0])
      {
        for (s0 in 1:p)
        {
          grad_res = c(grad_res, grad_cor(h0,l0,s0))
        }
      }
    }
  } else if (models == 1){
    ### derivative for correlation parameters theta_(h)_l for h=1,...,q, l=2,...,m_h,
    grad_cor1 <- function(h,l)
    {
      ### derative function of theta_0_s for wi, wj
      gradf_corhs <- function(w12,h,l)
      {
        x1 = w12[1:p]
        z1 = w12[(p+1):(p+q)]
        x2 = w12[(p+q+1):(p+q+p)]
        z2 = w12[(p+q+p+1):(2*p+2*q)]
        if((z1[h] != l) | (z2[h] != l)){
          return(0)
        }
        else{
          #variance parameter sigma^2
          par1 = parv[1:(q+1)]
          #correlation parameter in G1 to Gq
          par3 = parv[(q+2+p): npar]
          gx = - par1[h+1] * (sum((x1 - x2)^2)) * exp(psum(x1,x2, par3[sum((m-1)[1:h]) - (m-1)[h] + l-1]))
          return(as.numeric(gx))
        }
      }
      der_m = matrix(0,n,n)
      Rtemp_cor <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = gradf_corhs, h=h, l=l)
      der_m[rcoord] <- Rtemp_cor
      der_m <- der_m + t(der_m)
      diag(der_m) = 0
      result = sum(diag(invc %*% der_m)) - t(y-mu) %*% invc %*% der_m %*% invc %*% (y-mu)
      return(result)
    }
    ###analytical gradients results
    grad_res = grad_var0()
    for (i in 1:q) {
      grad_res = c(grad_res, grad_var(i))
    }
    for (i in 1:p)
    {
      grad_res = c(grad_res, grad_cor0(i))
    }
    for (h0 in 1:q)
    {
      for (l0 in 2:m[h0])
      {
        grad_res = c(grad_res, grad_cor1(h0,l0))
      }
    }
  }

  return( list( "objective" = MLE_result,
                "gradient" = grad_res))
}
