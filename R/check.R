#'CHECKS FOR CONFORMAL INFERENCE INTERVAL PREDICTION
#'
#' It contains all the check functions in the package
#' All the arguments are identical to the inputs of conformal.split.fun in the file
#' split.R
#' @param x The input variable, a list of n elements. Each element is composed by a list
#'  of p vectors(with variable length, since the evaluation grid may change).
#'  If x is NULL, the function will sample it from a gaussian.
#' @param t_val The grid points for the evaluation of function y_val. It is a list of vectors.
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
#' @param seed Integer to be passed to set.seed before defining the random
#'   data-split to be used. Default is FALSE, which effectively sets no seed.
#'   If both split and seed are passed, the former takes priority and the latter
#'   is ignored.
#' @param randomized Should the randomized approach be used? Default is FALSE.
#' @param seed.tau The seed for the randomized version.Default is FALSE.
#' @param training_size Split proportion between training and calibration set.
#' Default is 0.5.
#' @noRd

check.args=function(x,t_val,y,x0,train.fun,
               predict.fun, alpha, seed, training_size, seed.tau, randomized){


if(!is.null(x)){
  if ( is.list(x)==FALSE || is.data.frame(x)==TRUE || is.list(x[[1]])==FALSE
      || is.data.frame(x[[1]])==TRUE){



    stop("x must be a list of lists. Specifically, x must be a list of 'n' lists.
    Each of the 'n' lists must be made up of 'p' lists.
         'n' (i.e. the sample size) must be greater or equal than 2.
         'p'(i.e. the dimension of the multivariate function ) must be the same
         for all the multivariate functions.")}

}

  if (is.list(t_val)==FALSE || is.data.frame(t_val)==TRUE ){

    stop("t must be a list.
    Each of the 'p' lists must contain a numeric vector expressing the evaluation of the
    function on a grid (whose length can be different in the 'p' dimensions).
         'p'(i.e. the dimension of the multivariate function ) must be the same
         for all the multivariate functions.")}


  if (is.list(y)==FALSE || is.data.frame(y)==TRUE || is.list(y[[1]])==FALSE
      || is.data.frame(y[[1]])==TRUE){

    stop("y must be a list of lists. Specifically, y must be a list of 'n' lists.
    Each of the 'n' lists must be made up of 'p' lists.
    Each of the 'p' lists must contain a numeric vector expressing the evaluation of the
    function on a grid (whose length can be different in the 'p' dimensions).
         'n' (i.e. the sample size) must be greater or equal than 2.
         'p'(i.e. the dimension of the multivariate function ) must be the same
         for all the multivariate functions.")}


  if(!is.null(x0)){
  if (is.list(x0)==FALSE || is.data.frame(x0)==TRUE || is.list(x0[[1]])==FALSE
      || is.data.frame(x0[[1]])==TRUE){

    stop("x0 must be a list of lists. Specifically, x0 must be a list of 'm' lists.
    Each of the 'n' lists must be made up of 'p' lists.
         'm' (i.e. the sample size) must be greater or equal than 1.
         'p'(i.e. the dimension of the multivariate function ) must be the same
         for all the multivariate functions.")}
}

  if (length(unique(vapply(y,length,integer(1))))!=1)
    stop("'p'(i.e. the dimension of the multivariate function) must be
         the same for all the multivariate functions.")

  p=length(y[[1]])

 if(!(all(apply(t(vapply(y,function(x) vapply(x,length,integer(1)),integer(p))),2,
                function(y) length(unique(y))==1))))
  stop("The 'n' functions must be evaluated on the same p-variate grid.
       The grid can vary between the p domains.")


 if (!is.null(x) && length(x) != length(y)) {

   stop("length(x) and length(y) must match")
   }


  if (is.null(train.fun) || !is.function(train.fun))
    stop("train.fun must be a function")


  if (is.null(predict.fun) || !is.function(predict.fun))
    stop("predict.fun must be a function")


  check.num.01(alpha)


  if (is.null(seed)==TRUE || (seed!=FALSE & is.numeric(seed)==FALSE))
    stop("Argument 'seed' must be either FALSE or an integer.")



  if (is.null(training_size)==TRUE || (training_size!=FALSE & is.numeric(training_size)==FALSE)) stop("Argument 'training_size' must be either FALSE or an integer.")


  check.num.01(training_size)



  if (is.null(randomized)==TRUE || randomized %in% c("TRUE","FALSE")==FALSE)
    stop("Argument 'randomized' must be either TRUE or FALSE")

}




check.null.data=function(y){


  if (is.null(y)) stop("y must be either be punctual evaluation or fda
                                     object")
}



check.bool = function(b) {
  if (is.null(b) || length(b)!=1 || !is.logical(b))
    stop(paste(deparse(substitute(b)),"must be a Boolean"))
}

check.num = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a))
    stop(paste(deparse(substitute(a)),"must be a number"))
}

check.int = function(i) {
  if (is.null(i) || length(i)!= 1 || !is.numeric(i) || round(i) != i)
    stop(paste(deparse(substitute(i)),"must be an integer"))
}

check.pos.num = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a) || a<0)
    stop(paste(deparse(substitute(a)),"must be a positive number"))
}

check.pos.int = function(i) {
  if (is.null(i) || length(i)!= 1 || !is.numeric(i) || round(i) != i || i<1)
    stop(paste(deparse(substitute(i)),"must be a positive integer"))
}

check.num.01 = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a) || a<0 || a>1)
    stop(paste(deparse(substitute(a)),"must be a number between 0 and 1"))
}

check.tau=function(alpha,tau,l){

  if (alpha<tau/(l+1) & alpha>0)
    stop ("The prediction band obtained with such a small value of alpha is the entire space.
                                       If you are using the non randomized version of the algorithm, non-trivial prediction bands can be obtained by using an alpha greater or equal than 1/(l+1) and less than 1.
                                       If you are using the randomized version of the algorithm, non-trivial prediction bands can be obtained by using an alpha greater or equal than tau/(l+1) and less than (tau+l)/(l+1).")


  if (alpha>=(l+tau)/(l+1) || alpha<=0)
    stop("The alpha value is not admissible.
                                                   If you are using the non randomized version of the algorithm, non-trivial prediction bands can be obtained by using an alpha greater or equal than 1/(l+1) and less than 1.
                                                   If you are using the randomized version of the algorithm, non-trivial prediction bands can be obtained by using an alpha greater or equal than tau/(l+1) and less than (tau+l)/(l+1).")

}

check.s_regression=function(vec_residuals,type){
  #check on 'vec_residual' argument
  if( is.matrix(vec_residual)==FALSE & is.data.frame(vec_residual)==FALSE & (is.atomic(vec_residual)==FALSE || is.vector(vec_residual)==FALSE)) stop("vec_residual must be either a matrix, a dataframe or an atomic vector (naive case).")



  #check on 'type' argument
  possible_s_functions=c("identity","st-dev","alpha-max")
  if (is.null(type) || type %in% possible_s_functions==FALSE) {
    stop(c("The 'type' argument is not correct. Please select one of the following:",paste(possible_s_functions,collapse=", "),"."))
  }

}
