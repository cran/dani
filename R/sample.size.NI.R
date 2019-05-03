sample.size.NI <-
  function(p0.expected, p1.expected, p1.tolerable, sig.level=0.025, power=0.9, r=1, scale="RD", print.out=TRUE) {
    stopifnot((( scale == "RD" ) || ( scale == "RR" ) || ( scale == "AS" )), p0.expected <=1, p1.expected<=1,
              p0.expected>=0, p1.expected>=0, sig.level<1, power<1, sig.level>0, power>0, r>0, p1.tolerable<=1, p1.tolerable>=0)
    if (scale == "RD") {
      NIm <- abs(p1.tolerable-p0.expected)
      n <- ( qnorm( 1-sig.level )+qnorm( power ))^2 * (( p0.expected * ( 1 - p0.expected ) + p1.expected * ( 1 - p1.expected )/r)/( p1.expected - p0.expected - NIm )^2)
      ss <- c( nC <- ceiling( n ), nA <- ceiling ( ceiling( n ) * r ))
      if ( print.out == T ) {
        cat("Power:", ( power ) * 100, "%\nOne-sided significance level:", sig.level*100, "%.\nExpected control event risk =", p0.expected*100, "%\nExpected active event risk =",p1.expected*100,"%\nThe sample size required to test non-inferiority within a", NIm*100, "% risk difference NI margin is:")
        cat("\n", ss[1], " individuals in the control group.\n", ss[2], " individuals in the active group.\n")
      }
     } else if (scale=="RR") {
       NIm <- abs(log(p1.tolerable/p0.expected))
       n= ( qnorm( 1 - sig.level ) + qnorm( power ))^2 * ((( 1 - p0.expected )/( p0.expected )+( 1 - p1.expected )/( r * p1.expected ))/( log( p1.expected / p0.expected ) - NIm )^2)
       ss <- c( nC <- ceiling( n ), nA <- ceiling ( ceiling( n ) * r ))
       if ( print.out == T ) {
        cat("Power:", ( power ) * 100, "%\nOne-sided significance level:", sig.level*100, "%.\nExpected control event risk =", p0.expected*100, "%\nExpected active event risk =",p1.expected*100,"%\nThe sample size required to test non-inferiority within a", NIm, "log-risk ratio NI margin is:")
        cat("\n", ss[1], " individuals in the control group.\n", ss[2], " individuals in the active group.\n")
       }
     } else if (scale=="AS") {
       NIm <- asin(sqrt(p1.tolerable))-asin(sqrt(p0.expected))
       n= ( qnorm( 1 - sig.level ) + qnorm( power ))^2 * ((( 1/( 4 * r ) + 1/4 ))/( asin( sqrt( p1.expected )) - asin( sqrt( p0.expected )) - NIm )^2)
       ss <- c( nC <- ceiling( n ), nA <- ceiling ( ceiling( n ) * r ))
       if ( print.out == T ) {
        cat("Power:", ( power ) * 100, "%\nOne-sided significance level:", sig.level*100, "%.\nExpected control event risk =", p0.expected*100, "%\nExpected active event risk =",p1.expected*100,"%\nThe sample size required to test non-inferiority within a", NIm, "arc-sine difference NI margin is:")
        cat("\n", ss[1], " individuals in the control group.\n", ss[2], " individuals in the active group.\n")
      }
    }
    return(ss)
  }

