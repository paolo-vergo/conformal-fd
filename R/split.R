#'Functional Split conformal prediction intervals.
#'
#' Compute prediction intervals using split conformal inference.
#'
#' @param x The input variable, a list of n elements. Each element is composed by a list
#'  of p vectors(with variable length, since the evaluation grid may change).
#'  If x is NULL, the function will sample it from a gaussian.
#' @param t The grid points for the evaluation of function y_val. It is a list of vectors.
#' If the y_val data type is "fData" or "mfData" is must be NULL.
#' @param y The response variable. It is either, as with x and t, a list of list of
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
#' @param seed_tau The seed for the randomized version.Default is FALSE.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#' @param training_size Split proportion between training and calibration set.
#' Default is 0.5.
#' @param s_type The type of modulation function.
#'  Currently we have 3 options: "identity","st-dev","alpha-max". Default is "std-dev".
#'
#' @return A list with the following components: t,pred,k_s,s_type,s,alpha,randomized,tau,
#'  extremes_are_included,average_width,product_integral. t and s are lists of vectors,
#' pred has the same interval structure of y_val, but the outside list is of length n0,
#' k_s, average_width and product_integral are all positive floats, alpha and tau are
#' positive floats less than 1, randomized and extremes_are_included are logical values,
#' while s_type is a string.
#'
#' @references The function structure is taken from "Conformal Prediction Bands
#' for Multivariate Functional Data" by Diquigiovanni, Fontana, Vantini (2021) and,
#' also, from "The Importance of Being a Band: Finite-Sample Exact Distribution-Free
#' Prediction Sets for Functional Data" by Diquigiovanni, Fontana, Vantini (2021).
#'
#' @example inst/examples/ex.split.R
#' @export


conformal.fun.split = function(x, t, y, x0, train.fun, predict.fun, alpha=0.1,
                             split=NULL, seed=FALSE, randomized=FALSE,seed_tau=FALSE,
                                verbose=FALSE, training_size=0.5,s_type="st-dev") {


  ############ DATA PREPARATION #############################
  check.null.data(y)
  conv = convert2data(t,y,x,x0)

  x = conv$x
  y = conv$y
  t = conv$t
  x0 = conv$x0


  n=length(y)
  p=length(y[[1]])
  grid_size=vapply(y[[1]],function(x) length(x),integer(1))

  # Check input arguments
  check.args(x=x,t_val=t,y=y,x0=x0,train.fun=train.fun,
             predict.fun=predict.fun, alpha=alpha, seed=seed, training_size
             =training_size, seed.tau=seed.tau, randomized=randomized)


  # Users may pass in a string for the verbose argument
  if (verbose == TRUE) txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }

  ## Splitting Data


  if(is.null(split)){

    if(ceiling(n*training_size) !=n )
      m=ceiling(n*training_size)
    else
      m=ceiling(n*training_size)-1

    l=n-m

   if(seed!=FALSE){set.seed(seed)}

    training=sample(1:n,m)

  }

  else{
    training=split
  }

  calibration=setdiff(1:n,training)

  if(randomized==FALSE) {tau=1} else{
    if(seed_tau!=FALSE){set.seed(seed_tau)}
    tau=stats::runif(n=1,min=0,max=1)
  }

  check.tau(alpha=alpha,tau=tau,l=l)


  # TRAINING & RESIDUALS COMPUTATION ######################

  if (verbose) {
    cat(sprintf("%sComputing models on first part ...\n",txt))
  }


  out = train.fun(x[training,drop=F],t,y[training])
  fit = predict.fun(out,x,t)

  pred = predict.fun(out,x0,t)


  if (verbose) {
    cat(sprintf("%sComputing residuals and quantiles on second part ...\n",txt))
  }


  vect_y=t(vapply(y,unlist,numeric(sum(grid_size))))
  vect_fit=t(vapply(fit,unlist,numeric(sum(grid_size))))

  vect_residuals_y=vect_y-vect_fit


  if (verbose) {
    cat(sprintf("%sComputing modulation function ...\n",txt))
  }

  s=computing_s_regression(vec_residual=vect_residuals_y[training,],type=s_type,
                           alpha=alpha,tau=tau,grid_size=grid_size)
  vect_s=unlist(s)



  ############## COMPUTATION OF THE NCS #################


  rho=apply(vect_residuals_y[calibration,],1,function(x) max(abs(x)/vect_s))
  k_s=sort(rho,decreasing=FALSE)[ceiling(l+tau-(l+1)*alpha)]


  if ((ceiling(l+tau-(l+1)*alpha))==1) v=0
  else{v=sum(sort(rho,decreasing=FALSE)[1:(ceiling(l+tau-(l+1)*alpha)-1)]==k_s)}
  if ((ceiling(l+tau-(l+1)*alpha))==l) r=0
  else{r=sum(sort(rho,decreasing=FALSE)[(ceiling(l+tau-(l+1)*alpha)+1):length(rho)]==k_s)}

  extremes_are_included= tau > (alpha*(l+1)-floor(alpha*(l+1)-tau)+r)/(r+v+2)
  average_width=mean(2*k_s*vect_s)
  product_integral=exp(mean(log(2*k_s*vect_s)))


  ##################### BUILD BOUNDS #######################



  return(structure(.Data=list(t,pred,k_s,s_type,s,alpha,randomized,tau,
                              extremes_are_included,average_width,product_integral),
                   names=c("t","pred","k_s","s_type","s","alpha","randomized","tau"
                           ,"extremes_are_included","average_width","product_integral")))

}

utils::globalVariables(c("seed.tau", "td", "vec_residual"))
