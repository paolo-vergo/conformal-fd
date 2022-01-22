#' COMPUTING THE MODULATION FUNCTION S
#'
#' It computes modulation functions which allows local scaling of the prediction bands.
#'
#' @param vec_residual A vector of the residuals obtained via functional modeling.
#' @param type A string indicating the type of modulation function chosen.
#' The alternatives are "identity","st-dev","alpha-max".
#' @param alpha The value of the confidence interval.
#' @param tau A number between 0 and 1 used for the randomized version of the algorithm.
#' @param grid_size A vector containing the number of grid points in each dimension.
#' @return It returns a the values of a modulation function in each dimension of the response.
#' @details More details can be found in the help of conformal.fun.split function.
#' @references The function structure is taken from "Conformal Prediction Bands
#' for Multivariate Functional Data" by Diquigiovanni, Fontana, Vantini (2021) and,
#' also, from "The Importance of Being a Band: Finite-Sample Exact Distribution-Free
#' Prediction Sets for Functional Data" by Diquigiovanni, Fontana, Vantini (2021).
#' @importFrom stats sd
#' @export



computing_s_regression=function(vec_residual,type,alpha,tau,grid_size){


  if( is.matrix(vec_residual)==FALSE & is.data.frame(vec_residual)==FALSE & (is.atomic(vec_residual)==FALSE || is.vector(vec_residual)==FALSE)) stop("vec_residual must be either a matrix, a dataframe or an atomic vector (naive case).")



  #check on 'type' argument
  possible_s_functions=c("identity","st-dev","alpha-max")
  if (is.null(type) || type %in% possible_s_functions==FALSE) {
    stop(c("The 'type' argument is not correct. Please select one of the following:",paste(possible_s_functions,collapse=", "),"."))}



  indicator_grid=NULL
  for (i in 1:length(grid_size)) indicator_grid=c(indicator_grid,rep(i,grid_size[i]))

  #----naive cases: just one observation and type %in% c("st-dev","alpha-max")

  if(is.atomic(vec_residual)==TRUE & is.vector(vec_residual)==TRUE & type=="st-dev") stop("st-dev can not be computed when the number of observations is equal to 1.")
  if(is.atomic(vec_residual)==TRUE & is.vector(vec_residual)==TRUE & type=="alpha-max") {

    out=split(abs(vec_residual),indicator_grid)
    names(out)=paste("s_",1:length(grid_size),sep="")
    return(out)
  }
  #---non-naive cases

  if (type=="identity"){
    out=split(rep(1,sum(grid_size)),indicator_grid)
    names(out)=paste("s_",1:length(grid_size),sep="")
    return(out)
  }

  if (type=="st-dev") {
    out=split(apply(vec_residual,2,sd),indicator_grid)
    names(out)=paste("s_",1:length(grid_size),sep="")
    return(out)
  }

  if (type=="alpha-max"){


    check.num.01(tau)

    abs_vec_residual=abs(vec_residual)

    #----------------------------------------------CHECKS ON alpha-----------------------------------------------

    if(ceiling(dim(abs_vec_residual)[1]+tau-(dim(abs_vec_residual)[1]+1)*alpha) >= dim(abs_vec_residual)[1]) {
      out=split(apply(abs_vec_residual,2,max),indicator_grid)
      names(out)=paste("s_",1:length(grid_size),sep="")
      return(out)}
    if(ceiling(dim(abs_vec_residual)[1]+tau-(dim(abs_vec_residual)[1]+1)*alpha) <= 0)           {
      out=split(rep(1,sum(grid_size)),indicator_grid)
      names(out)=paste("s_",1:length(grid_size),sep="")
      return(out)}

    #----------------------------------------------S ALPHA-MAX----------------------------------------------------

    sequence_sup=apply(abs_vec_residual,1,max)
    gamma=sort(sequence_sup,decreasing=FALSE)[ceiling(dim(abs_vec_residual)[1]+tau-(dim(abs_vec_residual)[1]+1)*alpha)]
    position_functions_in_H=which(sequence_sup <= gamma)
    out=split(apply(abs_vec_residual[position_functions_in_H,],2,max),indicator_grid)
    names(out)=paste("s_",1:length(grid_size),sep="")
    return(out)

  }


}

