#' @title The Fitting Function of \code{LEzGP} Model
#'
#' @description Fits a Localized Easy-to-Interpret Gaussian process (LEzGP) model to a dataset as described in \code{reference 1}.
#'     The input variables are mixed (with both quantitative and qualitative inputs)
#'     The output variable is quantitative and scalar.
#'
#' @param X Matrix or data frame containing the inputs of training data. Each row represents the input setting of a data point and the columns are values of quantitative variables and qualitative variables.
#' @param Y Vector containing the outputs of training data points.
#' @param p Number of quantitative factors in the given dataset \code{X}.
#' @param q Number of qualitative factors in the given dataset \code{X}.
#' @param m A vector containing numbers of levels in qualitative factors.
#' @param tar_z A vector containing the qualitative part of the chosen target input (described in \code{reference 1}).
#' @param ns The chosen tuning parameter (described in \code{reference 1})
#' @param models The model for fitting the selected proper subset of the dataset \code{X}. 0 for EzGP model, 1 for EEzGP model.
#' @param tau Nugget if needed. The default nugget is 0, otherwise it has to be a non-negative real value.
#' @param lb Vector with lower bounds of the parameter estimation. "T" for applying the default setting of lb (a vector of length number of parameters whose elements are all 0.1), otherwise one must provide a vector with the length being the number of parameters.
#' @param ub Vector with upper bounds of the parameter estimation. "T" for applying the default setting of ub (a vector of length number of parameters whose first \code{q+1} elements are 100 and the rest \code{number of parameters - q - 1} elements are 10), otherwise one must provide a vector with length of the number of parameters.
#' @param x0 Vector with starting values for the optimization. "T" for applying the default setting of x0 (a vector made by \code{(lb + ub)/2}), otherwise one must provide a vector with length being the number of parameters.
#' @param xtol_rel Stopping criterion for relative change reached.
#' @param maxeval Termination condition by specifying a maximum number of function.
#' @param algorithm Optimization algorithm. See \href{https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/}{NLopt Algorithms} for more availiable algorithms.
#'
#' @import nloptr
#'
#' @return A model of class "LEzGP model" list of the following items:
#' \itemize{
#' \item{\code{param}} {A list containing the estimated parameters}
#' \item{\code{data}} {A list containing the fitted dataset and the information for fitting}
#' }
#'
#' @references
#' \enumerate{
#' \item "EzGP: Easy-to-Interpret Gaussian Process Models for Computer Experiments with Both Quantitative and Qualitative Factors", Qian Xiao, Abhyuday Mandal, C. Devon Lin, and Xinwei Deng (\doi{10.1137/19M1288462})
#' }
#'
#' @export
#'
#' @seealso
#' \code{\link[EzGP]{EzGP_predict}} to use the fitted EzGP model for prediction if your LEzGP model is fitted based on the EzGP model.\cr
#' \code{\link[EzGP]{EEzGP_predict}} to use the fitted EEzGP model for prediction if your LEzGP model is fitted based on the EEzGP model.\cr
#'
#' @examples
#' # Example with 9 quantitative and 9 qualitative variables (dataset included in the package):
#' #     Fit a LEzGP model based on the EEzGP/EzGP model(with default settings), and then
#' #     perform the prediction.
#' p = 9
#' q = 9
#' m=rep(3,9)
#' tau = 0
#' X = LEzGP_data[1:60, 1:(p+q)]
#' Y = LEzGP_data[1:60, p+q+1]
#' X_new = LEzGP_data[61:70, 1:(p+q)]
#' tar_z = X_new[1, (p+1):(p+q)]
#' ns = 7
#' # LEzGP Model Based on EEzGP Model
#' model <- LEzGP_fit(X, Y, p, q, m, tar_z, ns)
#' y_hat <- EEzGP_predict(X_new, model)
#' # Results showing
#' model
#' y_hat

LEzGP_fit <- function(X, Y, p, q, m, tar_z, ns, models = 1, tau = 0, lb = "T",
                      ub = "T",  x0 = "T", xtol_rel = 1.0e-5, maxeval = 100,
                      algorithm = "NLOPT_LD_LBFGS"){

  # function input checking
  if (is.vector(tar_z) == FALSE){
    stop('    tar_z must be a vector containing the qualitative part of the chosen target input, i.e., the length of the tar_z must equal to the number of qualitative factors.')
  }
  if (length(tar_z) != q){
    stop('    tar_z must be a vector containing the qualitative part of the chosen target input, i.e., the length of the tar_z must equal to the number of qualitative factors.')
  }
  if (!is.finite(ns) || !is.numeric(ns)){
    stop('    ns must be a finite number')
  }

  num = nrow(X)
  rows = c()
  i = 1
  for (j in 1:num){
    if (sum(X[j,(p+1):(p+q)] == tar_z) >= ns){
      rows[i] = j
      i = i + 1
    }
  }

  new_X <- X[rows, ]
  new_Y <- Y[rows]

  if (models == 1){
    return(EEzGP_fit(new_X, new_Y, p, q, m, tau, lb, ub, x0, xtol_rel, maxeval, algorithm))
  } else if (models == 0){
    return(EzGP_fit(new_X, new_Y, p, q, m, tau, lb, ub, x0, xtol_rel, maxeval, algorithm))
  }
}
