#' CONCURRENT MODEL FOR REGRESSION OVER FUNCTIONAL DATA
#'
#' This model is a concurrent model, which may be fed to split.conf.fun.
#'
#' @return A training and a prediction function.
#'
#' @details For more details about the structure of the inputs go to split.R
#' @importFrom stats lm
#' @export


concurrent = function() {

  train.fun = function(x,t,y) {


    yy=lapply(y, rapply, f = c) # Now a list of n components (join the internal p lists)
    xx=lapply(x, rapply, f = c)
    yyy=do.call(rbind, yy) #Convert the previous yy to a matrix
    xxx=do.call(rbind, xx)

    full = ncol(yyy)
    full_x = ncol(xxx)

    if(full!=full_x)
      stop(" The concurrent model requires a value of x for each value of y.
           If the number of dimension is diffent, then use another model.
           For instance the mean_fun model is available.")


    coeff=vapply(1:full, function(i)
      lm(formula = yyy[,i] ~  xxx[,i ])$coefficients,numeric(2))

    return(list(coeff=coeff))
  }

  # Prediction function
  predict.fun = function(out,newx,t) {

    #Redefine structure

    new_xx=lapply(newx, rapply, f = c)
    new_xxx=do.call(rbind, new_xx)
    temp=out$coeff


    l=length(newx)
    ya=temp[1,]
    yaa=t(matrix(replicate(l,ya),nrow=length(ya)))
    sol=new_xxx*temp[2,]+yaa


    list_sol=lapply(seq_len(nrow(sol)), function(i) list(sol[i,]))

    return(list_sol)
  }



  return(list(train.fun=train.fun, predict.fun=predict.fun))
}
