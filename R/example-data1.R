#' Dataset for the example in function 'EzGP_fit'
#'
#' Data are sampled from the modified math function based on Example 4.1 in the paper listed in \code{references}.
#' There are 3 quantitative factors and 3 qualitative factors each having 3 levels.
#' In this dataset, there are 1296 data points. For the simplicity of illustration, we take the first 81 rows as training data points, and the last 1215 rows as testing data points.
#'
#' @docType data
#'
#' @usage data(EzGP_data)
#'
#' @format A named list containing training data and testing data:
#' \describe{
#'   \item{"x1"}{1st quantitative factor}
#'   \item{"x2"}{2nd quantitative factor}
#'   \item{"x3"}{3rd quantitative factor}
#'   \item{"z1"}{1st qualitative factor, which has 3 levels}
#'   \item{"z2"}{2nd qualitative factor, which has 3 levels}
#'   \item{"z3"}{3rd qualitative factor, which has 3 levels}
#'   \item{"y"}{Response vector}
#' }
#'
#' @keywords example dataset
#'
#' @references
#' \enumerate{
#' \item "EzGP: Easy-to-Interpret Gaussian Process Models for Computer Experiments with Both Quantitative and Qualitative Factors", Qian Xiao, Abhyuday Mandal, C. Devon Lin, and Xinwei Deng (\doi{10.1137/19M1288462})
#' }
#'
#' @source The dataset can be generated with the code at the end of this description file.
#'
#' @examples
#' data(EzGP_data)
#' #Number of quantitative factors
#' p = 3
#' #Number of qualitative factors
#' q = 3
#' #Vector containing numbers of levels in qualitative factors
#' m=c(3,3,3)
#' # Nugget
#' tau = 0
#'
#' X = EzGP_data[1:81, 1:(p+q)]
#' Y = EzGP_data[1:81, p+q+1]
#' X_new = EzGP_data[82:1296, 1:(p+q)]
"EzGP_data"
