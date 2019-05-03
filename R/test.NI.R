test.NI <- function(n0, n1, e0, e1, NIm, sig.level=0.025, scale="RD", print.out=TRUE ) {
  stopifnot((( scale == "RD" ) || ( scale == "RR" ) || ( scale == "AS" )), 
            sig.level<1, sig.level>0, n0>0, n1>0, e0>=0, e1>=0)
  if (scale=="RD") {
    se <- sqrt(e1/n1*(1-e1/n1)/n1+e0/n0*(1-e0/n0)/n0)
    estimate <- e1/n1-e0/n0
    Z <- (estimate - NIm)/se
    p <- pnorm(Z)
    CI.norm <- c(estimate-qnorm(1-sig.level)*se,estimate+qnorm(1-sig.level)*se)
    # Method "Newcombe 10" implemented in ci.pd gives better CI:
    CI <- ci.pd(matrix(c(e0, e1, n0-e0,n1-e1),ncol=2, byrow=TRUE), alpha = sig.level*2, print=F)[6:7]
    if (print.out==T) cat("Risk difference:\nEstimate: ", estimate, "\nStandard error: ", se, "\nConfidence interval (two-sided ", sig.level*200,"%): (", CI[1], ",", CI[2], ")\np-value:", p, ".\n" )
    
  } else if (scale == "RR") {
    se <- sqrt(1/e0-1/n0+1/e1-1/n1)
    estimate <- log((e1/n1)/(e0/n0))
    Z <- (estimate - NIm)/se
    p <- pnorm(Z)
    CI.norm <- CI <- c(estimate-qnorm(1-sig.level)*se,estimate+qnorm(1-sig.level)*se)
    if (print.out==T) cat("Log-risk ratio:\nEstimate: ", estimate, "\nStandard error: ", se, "\nConfidence interval (two-sided ", sig.level*200,"%): (", CI[1], ",", CI[2], ")\np-value:", p, ".\n" )
    
  } else if (scale == "AS") {
    se <- sqrt(1/(4*n0)+1/(4*n1))
    estimate <- asin(sqrt(e1/n1))-asin(sqrt(e0/n0))
    Z <- (estimate - NIm)/se
    p <- pnorm(Z)
    CI.norm <- CI <- c(estimate-qnorm(1-sig.level)*se,estimate+qnorm(1-sig.level)*se)
    if (print.out==T) cat("Arc-sine difference:\nEstimate: ", estimate, "\nStandard error: ", se, "\nConfidence interval (two-sided ", sig.level*200,"%): (", CI[1], ",", CI[2], ")\np-value:", p, ".\n" )
  }
  results <- list(estimate, se, Z, p, CI, CI.norm)
  names(results)<-c("estimate", "se", "Z", "p", "CI", "CI.norm")
  return(results)
}