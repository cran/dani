\name{simulations.modify.margin}
\alias{simulations.modify.margin}
\title{
  Power and type 1 error of a "Modify margin" Non-inferiority trial design
}
\description{
  A function for running simulations to investigate power and type 1 error of a "Modify margin" Non-inferiority trial design. 
}
\usage{
  simulations.modify.margin(p0.expected, p1.expected, p1.tolerable, thresholds=c(0, Inf), 
            range.of.p0=NULL, sig.level.design=0.025, sig.level.analysis=0.025, 
            power=0.9, r=1, scale="RD", print.out=TRUE, ran.seed=1, n.sim=10000, 
            perf.measure="power") 
}

\arguments{
  \item{p0.expected}{Expected event risk in the control arm.}
  \item{p1.expected}{Expected event risk in the active arm.}
  \item{p1.tolerable}{Maximum tolerable event risk in the active arm.}
  \item{thresholds}{ A vector with the thresholds for modifying the non-inferiority margin to be compared. 
    Inf indicates the never modify margin procedure, while 0 the always modify margin. }
  \item{range.of.p0}{The values of true control event risk on which to explore the 
    inferential properties of the methods.}
  \item{sig.level.design}{One-sided significance level for testing used at the design stage. 
    Default is 0.025, i.e. 2.5\%.}
  \item{sig.level.analysis}{A vector with all the one-sided significance level(s) 
     used for testing at the analysis stage. 
    Default is a single level: 0.025, i.e. 2.5\%.}
  \item{power}{Power of the trial, i.e. one minus type-II error of the study. Default is 0.9, i.e.90\%.  }
  \item{r}{Allocation ratio, i.e. ratio between sample sizes in the active and control goups. Deafault is 1.}
  \item{scale}{The scale on which we define the non-inferiority trial. Can be one of "RD" (Risk      difference) or "RR" (log-risk ratio).}
  \item{print.out}{Logical. If FALSE, no output is printed.}
  \item{ran.seed}{The value to be used to initialise the random seed to obtain replicable results.}
  \item{n.sim}{Number of simulations to be performed. Defaults to 10000.}
  \item{perf.measure}{The performance measure to be used. It can be one of "power" or "type1error".}
}

\details{
  This is a function to perform simulations that can help the design of a non-inferiority trial by making it more resilient to unexpected control event risks. This is done by making use of the so-called "modify margin" procedure. The trial is designed as usual, using a certain value for the expected control event risk and defining the non-inferiority margin on either the risk difference or risk ratio scale. Then, when the trial is completed, the observed control event risk is compared with its a priori expectation and if the difference exceeds a certain threshold the margin is changed to a new value. The new value is chosen so that it reflects the so-called power-stabilising non-inferiority frontier (see paper in the references). 
At the design stage, this function can be used to find the optimal values of the threshold and of the significance level to be used in the analysis of the trial. First, it has to be run for different values of possible thresholds and for a fixed value of significance level for the analysis, equal to that used at the design stage. Then, if the preferred procedure inflates type 1 error in some areas, this function has to be re-run for different values of significance level, fixing the threshold to the one selected at the previous step. 
  
}

\value{
  The output is a list containing:
  
  - A matrix with the estimated power (or type-1 error) of each procedure in each scenario;
  
  - A matrix with the proportion of margins that have been modified using each procedure in each scenario;
  
  - A character variable, type, equal to "a" if the function has been used to find the optimal significance testing level and "t" if it has been used to find the optimal threshold;
  
  - A variable with the scale on which the non-inferiority margin was defined at the design stage;
  
  - A variable stating whetehr the function was used to estimate power or type 1 error.
}

\examples{
  
n.sim<-30 # Note in applications we would used higher values, i.e. >=10000
p0.expected<-0.05 # Expected control event rate
p1.expected<-p0.expected # Same as expected active event rate
p1.tolerable<-0.1  # Maximum tolerable active event rate 
r<-1       # Allocation ratio
power<-0.9 # Power
alph<-0.025 # Significance level
range.of.p0<-seq(0.005,0.20,0.005)

# Risk difference

res<-simulations.modify.margin(p0.expected, p1.expected, p1.tolerable, 
thresholds=c(Inf, 0.0125,0.025,0.05), range.of.p0=range.of.p0, sig.level.design=alph,
sig.level.analysis=alph, power=power, r=r, scale="RD", print.out=TRUE, 
ran.seed=1, n.sim=n.sim) 

alph.an<-c(0.01,0.015,0.02,0.025) # Significance level

res2<-simulations.modify.margin(p0.expected, p1.expected, p1.tolerable, thresholds=c(0.0125), 
range.of.p0=range.of.p0, sig.level.design=alph, sig.level.analysis=alph.an, 
power=power, r=r, scale="RD", print.out=TRUE, ran.seed=1, n.sim=n.sim) 

# Risk ratio

res3<-simulations.modify.margin(p0.expected, p1.expected, p1.tolerable, 
thresholds=c(Inf, log(1.25),log(1.5),log(2)), range.of.p0=range.of.p0,
sig.level.design=alph, sig.level.analysis=alph, power=power, r=r, 
scale="RR", print.out=TRUE, ran.seed=1, n.sim=n.sim) 

res4<-simulations.modify.margin(p0.expected, p1.expected, p1.tolerable, thresholds=log(1.25), 
range.of.p0=range.of.p0, sig.level.design=alph, sig.level.analysis=alph.an, 
power=power, r=r, scale="RR", print.out=TRUE, ran.seed=1, n.sim=n.sim) 


}