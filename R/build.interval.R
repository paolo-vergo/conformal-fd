

#'JOIN MULTIPLE INTERVALS
#'
#'The function build the multisplit interval for a given x0) in each point of the evaluation
#'grid in t, by combining multiple splits according to the Multi Split algorithm.
#'
#'
#'@param yyy column vector of B lower bounds and B upper bounds
#'@param B number of replications
#'@param tr truncation threshold for the algorithm
#'@return A single interval.
#'@export



interval.build=function(yyy,B,tr){



  h=rep(1:0,each=B)

  o = order(yyy,2-h)

  ys <- yyy[o]
  hs <- h[o]

  count <- 0
  leftend <- 0
  lo<-up<-0


  for (j in 1:(2*B) ){
    if ( hs[j]==1 ) {
      count <- count + 1

      if ( count > tr && (count - 1) <= tr) {
        leftend <- ys[j]
      }

    }

    else {
      if ( count > tr && (count - 1) <= tr) {
        rightend <- ys[j]
        lo <- leftend
        up <- rightend
      }

      count <- count - 1
    }
  }




  return(c(lo,up))
}
