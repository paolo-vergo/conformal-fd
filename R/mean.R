#' MEAN OF FUNCTIONAL DATA
#'
#' This model is a fed to a Functional Conformal Prediction function.
#'
#'
#' @return It outputs a training function and a prediction function.
#' @details For more details about the structure of the inputs go to split.fun.R
#'
#' @export


mean_lists = function() {

  # Training function
  train.fun = function(x,t,y,out=NULL) {

    n=length(y)
    p=length(y[[1]])
    grid_size=vapply(y[[1]],function(x) length(x),integer(1))
    acc_size=c(0,cumsum(grid_size))

    yy=lapply(y, rapply, f = c) # Now a list of n components
    #(join the internal p lists)
    yyy=do.call(rbind, yy) #Convert the previous yy to a matrix
    mmm=colMeans(yyy)

    m=lapply(1:p, function(i) {
      b=acc_size[i]+1
      e=acc_size[i+1]
      return(mmm[b:e])
    })


    return(list(m=m))
    }

  # Prediction function
  predict.fun = function(out,newx,t) {

    temp=out$m

    l=length(newx)
    sol=list(temp)[rep(1,l)]

    return(sol)
  }



  return(list(train.fun=train.fun, predict.fun=predict.fun))
}
