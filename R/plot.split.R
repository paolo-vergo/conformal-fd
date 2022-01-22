#' PLOT FUNCTIONAL SPLIT CONFORMAL CONFIDENCE BANDS
#'
#' The function plots the confidence bands provided by the conformal.fun.split function.
#'
#' @param out It's the output of the split function.
#' @return A list of of list composed by ggplots (output in position (ii,jj) is the
#' plot of the jj-th observation in the ii-th dimension of the multivariate functional response).
#' @details It exploits the package \code{\link[ggplot2]{ggplot}} and \code{\link[gridExtra]{grid.arrange}}
#' to better visualize the results. It outputs n0=length(x0) plots.
#' @details It plots, for each value in x0, the predicted functional value and bands in all the dimensions of the multivariate functional response.
#' @example inst/examples/ex.split.R
#' @export plot_fun

plot_fun=function(out){

  tn = out$t
  pred = out$pred
  k_s = out$k_s
  s = out$s

  n=length(pred)
  p=length(pred[[1]])
  grid_size=vapply(pred[[1]],function(x) length(x),integer(1))
  acc_size=c(0,cumsum(grid_size))
  vect_s=unlist(s)
  nCol <- floor(sqrt(p))
  g_list<-vector("list",n*p)


  for(jj in 1:n){
      for(ii in 1:p){

      s_index=(acc_size[ii]+1):(acc_size[ii+1])
      df=data.frame(td=as.numeric(tn[[ii]]),yg=as.numeric(pred[[jj]][[ii]]),y_min=as.numeric(pred[[jj]][[ii]]-k_s*vect_s[s_index]), y_max=as.numeric(as.numeric(pred[[jj]][[ii]])+k_s*vect_s[s_index]))

      g_list[[jj]][[ii]]<-ggplot2::ggplot()+ggplot2::geom_line(data=df,ggplot2::aes(x=td,y = yg),size=1) +
        ggplot2::geom_ribbon(data=df, ggplot2::aes(x=td,y = yg,ymin = y_min, ymax = y_max), fill = "red", alpha = 0.2)+
        ggplot2::labs(title = paste("Dimension ",as.character(ii), "Observation ", as.character(jj)),
             x = "t", y = "y")


    }

    do.call(gridExtra::"grid.arrange", c(g_list[[jj]], ncol=nCol))

  }

  return(g_list)
}

utils::globalVariables(c("y_max", "y_min", "yg"))

