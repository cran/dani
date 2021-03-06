\name{test.NI}
\alias{test.NI}
\title{
  Non-inferiority test on different scales
}
\description{
  A function for testing non-inferiority on either the risk difference, log-risk ratio or arc-sine difference scale. 
}
\usage{
  test.NI(n0, n1, e0, e1, NIm, sig.level=0.025, scale="RD", print.out=TRUE)  
}
\arguments{
  \item{e0}{Number of events in the control arm.}
  \item{e1}{Number of events in the active arm.}
  \item{n0}{Total sample size of the control arm.}
  \item{n1}{Total sample size of the active arm.}
  \item{NIm}{Non-inferiority margin on the specified scale.}
  \item{sig.level}{One-sided significance level for testing. Default is 0.025, i.e. 2.5\%.}
  \item{scale}{The scale on which we define the non-inferiority margin. Can be one of "RD" (Risk      difference), "RR" (log-risk ratio) or "AS" (arc-sine difference).}
  \item{print.out}{Logical. If FALSE, no output is printed.}
}
\details{
  This is a function to test non-inferiority of an active treatment against the control within a specific NI margin. The margin can be specified on a number of different scales. 
}

\value{
  The output is a list, containing the estimate, standard error, cofidence interval (two-sided 2*alpha level), Z statistic and p-value. For the riskdifference scale, the CI is provided both with the standard normal approximation and using newcombe 10 method, as per function ci.pd in package Epi.
  
}

\examples{
  
 n0<-1000
n1<-1000
e0<-0.05*n0
e1<-0.05*n1
NImRD=0.05
NImRR=log(2)
NImAS=asin(sqrt(0.1))-asin(sqrt(0.05))
alpha=0.025

test<-test.NI(n0=n0, n1=n1, e0=e0, e1=e1, NIm=NImRD, sig.level=alpha, scale = "RD")
test2<-test.NI(n0=n0, n1=n1, e0=e0, e1=e1, NIm=NImRR, sig.level=alpha, scale = "RR")
test3<-test.NI(n0=n0, n1=n1, e0=e0, e1=e1, NIm=NImAS, sig.level=alpha, scale = "AS")

  
}