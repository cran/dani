figure.modify.margin <- function (simulations, nominal=NULL, col=NULL, pch=16, lty=NULL, xlim=NULL, ylim=NULL) {
  
  if (is.null(col)) col <- 1:(ncol(simulations[[1]])-1)
  if (length(col)==1) col <- rep(col, ncol(simulations[[1]])-1 )
  if (is.null(lty)) lty <- 1:(ncol(simulations[[1]])-1)
  if (length(lty)==1) lty <- rep(lty, ncol(simulations[[1]])-1 )
  if (length(pch)==1) pch <- rep(pch, ncol(simulations[[1]])-1 )
  if (is.null(xlim)) xlim <- c(min(simulations[[1]][,ncol(simulations[[1]])]),max(simulations[[1]][,ncol(simulations[[1]])]))
  if (is.null(ylim)) ylim <- c(min(simulations[[1]][,1:(ncol(simulations[[1]])-1)]),max(simulations[[1]][,1:(ncol(simulations[[1]])-1)]))
  plot(simulations[[1]][,ncol(simulations[[1]])], simulations[[1]][,1], type = "b", 
       pch=pch[1], col=col[1], lty=lty[1], xlim=xlim, ylim=ylim,
       xlab = expression(pi[0]), ylab= names(simulations)[1], 
       main = ifelse(simulations[[4]]=="RD", "Risk Difference", "log-Risk Ratio"))
  for (i in 2:(ncol(simulations[[1]])-1)) {
    lines(simulations[[1]][,ncol(simulations[[1]])], simulations[[1]][,i],  col=col[i], type = "b", pch=pch[i], lty=lty[i])
  }
  lines(simulations[[1]][,ncol(simulations[[1]])], simulations[[1]][,1],  col=col[1], type = "b", pch=pch[1], lty=lty[1])
  
  if (!is.null(nominal)) {
    abline(h=nominal, col="red")
  }  

}