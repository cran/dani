\name{sample.size.NI}
\alias{sample.size.NI}
\title{
  Sample size calculation tool for Non-Inferiority trials
}
\description{
A function for calculating sample size of a non-inferiority trial. 
}
\usage{
sample.size.NI(p0.expected, p1.expected, p1.tolerable, sig.level=0.025, power=0.9, r=1, 
  scale="RD", print.out=TRUE)  
}
\arguments{
\item{p0.expected}{Expected event risk in the control arm.}
\item{p1.expected}{Expected event risk in the active arm.}
\item{p1.tolerable}{Maximum tolerable active event risk for declaring non-inferiority.}
\item{sig.level}{One-sided significance level for testing. Default is 0.025, i.e. 2.5\%.}
\item{power}{Power of the trial, i.e. one minus type-II error of the study. Default is 0.9, i.e.90\%.  }
\item{r}{Allocation ratio, i.e. ratio between sample sizes in the active and control goups. Deafault is 1.}
\item{scale}{The scale on which we define the non-inferiority trial. Can be one of "RD" (Risk      difference), "RR" (log-risk ratio) or "AS" (arc-sine difference).}
\item{print.out}{Logical. If FALSE, no output is printed.}
}
\details{
 This is a function to calculate sample size needed to test non-inferiority of an active treatment against the control within a specific NI margin. The margin can be specified on a number of different scales, which give different answers. 
 }
 
\value{
The output is a vector ss, containing the sample sizes for the control and active arms respectively.
}

\examples{

  p0.expected<-0.05 # Expected control event rate
  p1.expected<-p0.expected # Same as expected active event rate
  p1.tolerable<-0.1  # Maximum tolerable active event rate 
  r<-1       # Allocation ratio
  power<-0.9 # Power
  alph<-0.025 # Significance level

sample.size.RD<-sample.size.NI(sig.level=alph, power=power, 
   p0.expected=p0.expected, p1.expected=p1.expected, 
   p1.tolerable=p1.tolerable, r=r)  # Risk difference scale
sample.size.lRR<-sample.size.NI(sig.level=alph, power=power, 
   p0.expected=p0.expected, p1.expected=p1.expected, 
   p1.tolerable=p1.tolerable, r=r, scale="RR")  # Log-risk ratio scale
sample.size.AS<-sample.size.NI(sig.level=alph, power=power, 
   p0.expected=p0.expected, p1.expected=p1.expected, 
   p1.tolerable=p1.tolerable, r=r, scale="AS")  # Arc-sine difference scale
  
}