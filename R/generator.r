#' Generate survival samples
#'
#' This function generates survival data including response, covariates and status.
#' Data can be right-censored and heavy-tailed.
#' Morever, there may be a strong network structure and high-correlated pattern among covariates.
#' @param n The number of samples.
#' @param p The number of covariates.
#' @param d The number of nonzero values of true coefficients.
#' @param g The group size.
#' @param sig The two parameters of the uniform distribution. For example, \code{sig=c(a,b)}, which means the elements of true coefficients are from \code{U(a,b)}.
#' @param rho The strength of correlation among covariates.
#' @param error The error distribution. Quantitative for error="norm", error="contaminate", error="t2" or error="ev". For error="contaminate", y is from \code{0.7N(0,1)+0.3Cauchy}.
#' @param tran The transformation function. Quantitative for tran="lm" or tran="log"
#' @param censor.rate The censor rate of survival sample. Default is 0.
#' @param block The network structure of predictors. Quantitative for block="Block", block="Auto", block="BAND". For block="BAND", the number of groups \code{g} is equal to \code{p}.
#'
#' @return A list
#' \itemize{
#'   \item x - The \code{n x p} design matrix.
#'   \item y - The \code{n} repsonse.
#'   \item b - The true coefficients.
#'   \item status - The survival status. \code{1} is censored and \code{0} is uncensored.
#' }
#' @examples
#' library(mvtnorm)
#' library(Matrix)
#'
#' n = 400
#' p = 500
#' d = 15
#' g = 5
#' sig = c(0.5, 1.5)
#' rho = 0.9
#' error = "contaminate"
#' tran = "log"
#' censor.rate = 0.1
#' block = "Auto"
#'
#' dat = generator(n, p, d, g, sig, rho, error, tran, censor.rate, block)
#' x = dat$x
#' y = dat$y
#' b = dat$b
#' status = dat$status
#' @export

generator <- function(n, p, d, g, sig, rho, error = c("norm", "contaminate", "t2", "ev"), tran = c("lm", "log"),
                      censor.rate = 0, block = c("Block", "Auto", "BAND"))
{
  if(d%%g!=0){
    print("Warning:the number of nonzero coefficients d must be divisible by the group size!"); break
  }
  if(p%%g!=0){
    print("Warning:the number of coefficients p must be divisible by the group size!"); break
  }
  block = match.arg(block)
  error = match.arg(error)
  tran  = match.arg(tran)
  sig.min = sig[1]
  sig.max = sig[2]

  # x
  if(block=="Block"){
    s = matrix(rho,g,g)
    diag(s) = 1
    slist = NULL
    for (i in 1:(p/g)) { slist[[i]] = s }
    x = rmvnorm(n, sigma = as.matrix(bdiag(slist)))
  }
  if(block=="Auto"){
    s = outer(1:g, 1:g, FUN = function(x, y) rho^(abs(x - y)))
    slist = NULL
    for (i in 1:(p/g)) { slist[[i]] = s }
    x = rmvnorm(n, sigma = as.matrix(bdiag(slist)))
  }
  if(block=="BAND"){
    s = outer(1:p, 1:p, FUN = function(x, y) ifelse(abs(x-y)==1, rho^(abs(x - y)), 0))
    diag(s)=1
    x = rmvnorm(n, sigma = s)
  }

  # true beta
  b = runif(d, sig.min, sig.max)
  if(d%%2==0) { sgn = rep(c(1,-1), d/2); b = sgn*b }
  if(d%%2!=0) { sgn = c(rep(c(1,-1), d%/%2), 1); b = sgn*b }
  b = c(b, rep(0, p-d))

  # error
  if(error == "norm") e = rnorm(n, 0, 1)
  if(error =="contaminate") e = c(rnorm(0.7*n), rcauchy(0.3*n))
  if(error == "t2") e = rt(n, 2)
  if(error == "ev") e = log(rweibull(n, shape = 1, scale = 1)) - digamma(1)

  # transformation
  if(tran == "lm") g = function(x) x
  if(tran == "log") g = function(x) exp(x)

  # y
  y0 = g(c(x %*% b)+e)

  # censor rate
  if(censor.rate!=0){
    cens = quantile(y0, 1-censor.rate)
    y = pmin(y0, cens)
    status = 1 * (y0<cens)
  }else{
    y = y0
    status = rep(1, n)
  }

  # In the transformation model, y may be infinity.
  y[which(y=="NaN"|y=="Inf")] =  max(y[which(y!="NaN"&y!="Inf")])

  val = list(x = x,
             y = y,
             b = b,
             status = status)

  return(val)
}














