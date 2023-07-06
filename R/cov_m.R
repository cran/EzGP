#' @title The Function for Constructing the Covariance Matrix in \code{EzGP} Package
#'
#' @description Builds the covariance matrix for the given dataset according to different models.
#'
#' @param X Matrix or data frame containing the inputs of training data. Each row represents the input setting of a data point and the columns are values of quantitative variables and qualitative variables.
#' @param p Number of quantitative factors in the given dataset \code{X}.
#' @param q Number of qualitative factors in the given dataset \code{X}.
#' @param m A vector containing numbers of levels in qualitative factors.
#' @param n Number of training data points
#' @param parv Parameters in the EzGP/EEzGP model
#' @param tau Nugget if needed. The default nugget is 0, otherwise it has to be a non-negative real value.
#' @param models Model indicator that indicates which model the covariance matrix is built for. 0 for EzGP model, 1 for EEzGP model. The default setting is 0.
#'
#' @return The covariance matrix for the given dataset.
#'
#' @details \code{\link[EzGP]{EzGP_fit}}, \code{\link[EzGP]{EzGP_predict}}, \code{\link[EzGP]{EEzGP_fit}}, \code{\link[EzGP]{EEzGP_predict}}, \code{\link[EzGP]{LEzGP_fit}}, and \code{\link[EzGP]{LLF_gradients}} will call this function.
#'
#' @note This function is used inside other functions in this package and is \strong{NOT} exported once the EzGP package is loaded.
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
#'
#' @examples
#' # see the examples in the documentation of the function EzGP_fit.

cov_m <- function(X, p, q, m, n, parv, tau = 0, models = 0){

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

    if (models == 0){

      for (i in 1:q){
        if(z1[i] != z2[i]){
          res1 = res1+0
        }
        else{
          l = z1[i]
          res1 = res1 + par1[i+1]*exp(psum(x1,x2, par3[(sum(m[1:i])*p - m[i]*p + (l-1)*p + 1) : (sum(m[1:i])*p - m[i]*p + (l-1)*p + p)]))
        }
      }
      return(res1)

    } else if (models == 1){
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
  }

  ### a modified version of covx where w.12 = (w1, w2)
  covx.m <- function(w.12, parv){
    return( covx(w1 = w.12[1:(p+q)], w2 = w.12[(p+q+1):(2*p+2*q)], parv) )
  }

  #expand grid to avoid for looping
  rcoord <- cbind(
    rep(seq_len(n - 1), times = rev(seq_len(n - 1))),
    unlist(lapply(
      X = rev(seq_len(n - 1)),
      FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = n)))

  #covariance matrix
  covm = matrix(0,n,n)
  # first compute the vector of elements in covariance matrix
  Rtemp <- apply(cbind(X[rcoord[, 1], ], X[rcoord[, 2], ]), 1, FUN = covx.m, parv=parv)
  covm[rcoord] <- Rtemp
  covm <- covm + t(covm)
  diag(covm) <- sum(parv[1:(q+1)]) + tau

  return(covm)
}
