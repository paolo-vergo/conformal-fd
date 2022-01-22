#' Functional Multi Split conformal prediction intervals.
#'
#' Compute prediction intervals using functional multi split conformal inference.
#'
#' @param x The input variable, a list of n elements. Each element is composed by a list
#'  of p vectors(with variable length, since the evaluation grid may change).
#' @param t The grid points for the evaluation of function y_val. It has the same
#'  structure as x (list of list of vector).
#' @param y_val The response variable. It is either, as with x and t, a list of list of
#'  vectors or an fda object (of type fd, fData, mfData).
#' @param x0 The new points to evaluate, a list of n0 elements. Each element is composed
#'  by a list of p vectors(with variable length).
#' @param train.fun A function to perform model training, i.e., to produce an
#'   estimator of E(Y|X), the conditional expectation of the response variable
#'   Y given features X. Its input arguments should be x: list of features,
#'   and y: list of responses.
#' @param predict.fun A function to perform prediction for the (mean of the)
#'   responses at new feature values. Its input arguments should be out: output
#'   produced by train.fun, and newx: feature values at which we want to make
#'   predictions.
#' @param alpha Miscoverage level for the prediction intervals, i.e., intervals
#'   with coverage 1-alpha are formed. Default for alpha is 0.1.
#' @param split Indices that define the data-split to be used (i.e., the indices
#'   define the first half of the data-split, on which the model is trained).
#'   Default is NULL, in which case the split is chosen randomly.
#' @param seed Integer to be passed to set.seed before defining the random
#'   data-split to be used. Default is FALSE, which effectively sets no seed.
#'   If both split and seed are passed, the former takes priority and the latter
#'   is ignored.
#' @param randomized Should the randomized approach be used? Default is FALSE.
#' @param seed_beta The seed for the randomized version of the conformal.split.fun.
#' Default is FALSE.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#' @param training_size Vector containing the split proportion between training and calibration set.
#' It has B components. Default is 0.5.
#' @param s_type The type of modulation function.
#'  Currently we have 3 options: "identity","st-dev","alpha-max".
#' @param B Number of repetitions. Default is 100.
#' @param lambda Smoothing parameter. Default is 0.
#' @param tau It is a smoothing parameter:
#' tau=1-1/B  Bonferroni intersection method
#' tau=0 unadjusted intersection
#' Default is 1-(B+1)/(2*B).
#'
#' @return A list of the following: lo,up,grid_size,tn.
#'
#' @details The work is an extension of the univariate approach to Multi Split
#' conformal inference to a multivariate functional context.
#' @details This function is based on the package future.apply to
#'  perform parallelization. If this package is not installed, then the function
#'  will abort.
#'
#' @references "Multi Split Conformal Prediction" by Solari, Djordjilovic (2021) is
#' the baseline for the univariate case.
#'
#' @importFrom utils tail
#' @example inst/examples/ex.msplit.R
#' @export conformal.fun.msplit


conformal.fun.msplit = function(x,t, y_val, x0, train.fun, predict.fun, alpha=0.1,
                      split=NULL, seed=FALSE, randomized=FALSE,seed_beta=FALSE,
                      verbose=FALSE, training_size=NULL,s_type,B=100,lambda=0,
                      tau = 1-(B+1)/(2*B)) {



  if(is.null(training_size) || length(training_size)!=B)
    training_size=rep(0.5,B)

  if (!is.null(seed)) set.seed(seed)

  ### CONVERT DATA

  check.null.data(y_val)
  conv = convert2data(t,y_val,x,x0)
  tt = conv$t
  y = conv$y
  x = conv$x
  x0 = conv$x0


  n0=length(x0)
  grid_size=vapply(y[[1]],function(x) length(x),integer(1))
  dims=tail(cumsum(grid_size),1)


  future::plan(future::multisession)
  options(future.rng.onMisuse="ignore")


  Y_lo_up <- t(future.apply::future_sapply(1:B, function(bbb) {


    out<-conformal.fun.split(x,t, y_val, x0, train.fun, predict.fun, alpha*(1-tau) + (alpha*lambda)/B,
                             split, seed+bbb, randomized,seed_beta,
                             verbose, training_size[bbb] ,s_type)

    vect_s=c(unlist(out$s))
    rep_s=vect_s[rep(1:dims,n0)]

    pred=unlist(out$pred)

    lo<-pred-out$k_s*rep_s
    up<-pred+out$k_s*rep_s
    loup=cbind(t(lo),t(up))


    return(loup)
  }))




  full = ncol(Y_lo_up)/2
  Y = rbind(Y_lo_up[,1:full], Y_lo_up[,-(1:full)])

  tt=conformal.fun.split(x,t, y_val, x0, train.fun, predict.fun, alpha*(1-tau) + (alpha*lambda)/B,
                         split, seed+1, randomized,seed_beta,
                         verbose, training_size[1] ,s_type)$t


  finalInt = future.apply::future_sapply(1:full, function(kk) interval.build(Y[,kk],B,tr = tau*B + .001))
  lo<-split(finalInt[1,], rep(1:n0, each = dims))
  up<-split(finalInt[2,], rep(1:n0, each = dims))


  ## To avoid CRAN check errors
  ## R CMD check: make sure any open connections are closed afterward
  future::plan(future::sequential)





  return(list(lo=lo,up=up,grid_size=grid_size,tn=tt))
}

