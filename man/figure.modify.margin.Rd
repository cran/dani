\name{figure.modify.margin}
\alias{figure.modify.margin}
\title{
  Plot power and type 1 error of a "Modify margin" Non-inferiority trial design
}
\description{
  A function for plotting power and type 1 error of a "Modify margin" Non-inferiority trial design after running simulations . 
}
\usage{
  figure.modify.margin(simulations, nominal=NULL, col=NULL, 
              pch=16, lty=NULL, xlim=NULL, ylim=NULL) 
}

\arguments{
  \item{simulations}{A list obtained as an output from function simulations.modify.margin.}
  \item{nominal}{An optional argument to draw a red horizontal line specifying the nominal power or type 1 error.}
  \item{col}{A vector with colors of the different lines. See help file for the general plot function.}
  \item{pch}{ A vector with plotting characters of the different lines. See help file for the general plot function. }
  \item{lty}{A vector with line types of the different lines. See help file for the general plot function.}
  \item{xlim, ylim}{numeric vectors of length 2, giving the x and y coordinates ranges.}
 
 }

\details{
  This is a function to plot power (type 1 error) of the modify margin procedures for different values of actual control event risk. 
  
}


\examples{
  
  n.sim<-100 # In real applications use at least n.sim <- 10000
p0.expected<-0.05 # Expected control event rate
p1.expected<-p0.expected # Same as expected active event rate
p1.tolerable<-0.1  # Maximum tolerable active event rate 
NI.marg<-p1.tolerable-p0.expected
r<-1       # Allocation ratio
power<-0.9 # Power
alph<-0.025 # Significance level
alph.an<-c(0.01,0.015,0.02,0.025) # Significance level
range.of.p0<-seq(0.005,0.20,0.005)

# Power - risk difference

res<-simulations.modify.margin(p0.expected, p1.expected, p1.tolerable, 
thresholds=c(Inf, 0.0125,0.025,0.05),
range.of.p0=range.of.p0, sig.level.design=alph, sig.level.analysis=alph, power=power,
r=r, scale="RD", print.out=TRUE, ran.seed=1, n.sim=n.sim) 
  
figure.modify.margin(res, 0.9)

res2<-simulations.modify.margin(p0.expected, p1.expected, p1.tolerable, 
thresholds=c(Inf, 0.0125,0.025,0.05),
range.of.p0=range.of.p0, sig.level.design=alph, sig.level.analysis=alph, power=power,
r=r, scale="RD", print.out=TRUE, ran.seed=1, n.sim=n.sim, perf.measure="type1error") 
  
figure.modify.margin(res2, 0.025)
  
  
}