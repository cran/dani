simulations.modify.margin <-
  function(p0.expected, p1.expected, p1.tolerable, thresholds=c(0, Inf), range.of.p0=NULL, sig.level.design=0.025, sig.level.analysis=0.025, power=0.9, r=1, scale="RD", print.out=TRUE, ran.seed=1, n.sim=10000, perf.measure="power") {
    
    if ((length(thresholds)>1)&(length(sig.level.analysis)>1)) {
      stop("ERROR: only one between thresholds and sig.level.analysis should be a vector. The recommended procedure is to start by looking for the best threshold fixing sig.level.analysis=sig.level.design and then look for the best value of sig.level.design for the selected threshold only.")
    } else if ((length(thresholds)>1)) {
      type <- "t"
    } else if ((length(sig.level.analysis)>1)) {
      type <- "s"
    } else {
      type <- "o"
    }
    set.seed(ran.seed)
    stopifnot((( scale == "RD" ) || ( scale == "RR" )), (( perf.measure == "power" ) || ( perf.measure == "type1error" )),
              p0.expected <=1, p1.expected<=1, p1.tolerable<=1, p0.expected>=0, p1.expected>=0, p1.tolerable<=1, sig.level.design<1, power<1, 
              sig.level.design>0, power>0, r>0, n.sim>1, range.of.p0[1]>0, range.of.p0[length(range.of.p0)]<1)
    
    if (scale == "RD") {
      NI.marg <- p1.tolerable-p0.expected
    } else if ( scale == "RR" ) {
      NI.marg <- log(p1.tolerable/p0.expected)
    } 
    
    if (is.null(range.of.p0)) {
      range.of.p0<-seq(max(0.005,p0.expected-0.15),min(0.995,p0.expected+0.15),0.005)
    }
    n.pi0<-length(range.of.p0)
    
    if ( type == "t" ) {
      simulations<-array(NA, c(n.sim,5+length(thresholds)*2, n.pi0))
      Power.method<-matrix(NA,n.pi0,length(thresholds)+1)
      colnames(Power.method)<-c(paste(expression(epsilon), "=",thresholds), expression(pi[0]))
      n.modified<-matrix(NA,n.pi0,length(thresholds))
      colnames(n.modified)<-paste(expression(epsilon), "=",thresholds)
      sample.size<-ceiling(sample.size.NI( p0.expected = p0.expected, p1.expected = p1.expected, p1.tolerable = p1.tolerable, sig.level = sig.level.design, power = power, r = r, scale = scale, print.out = print.out ))  # Sample size
      
      for (j in 1:n.pi0) {
        pi0<-range.of.p0[j]  # True control event rate
        if (perf.measure == "power") {
          pi1<-pi0*(p1.expected/p0.expected)  # True active event rate
        } else {
          pi1<-sin(asin(sqrt(pi0))+asin(sqrt(p1.tolerable))-asin(sqrt(p0.expected)))^2  # True active event rate
        }
        
        # Data and results matrices
        
        simulations[,3,j]<-sample.size[1]
        simulations[,4,j]<-sample.size[2]
        
        for (i in 1:n.sim) {
          sample.size.obs<-sum(rbinom(1,sample.size[1]+sample.size[2],1/(1+r)))
          sample.size.obs<-c(sample.size.obs, sample.size[1]+sample.size[2]-sample.size.obs)
          simulations[i,1,j]<-rbinom(1,sample.size.obs[1],pi0)
          simulations[i,2,j]<-rbinom(1,sample.size.obs[2],pi1)
          simulations[i,5,j]<-test.NI(n0=sample.size.obs[1], n1=sample.size.obs[2], e0=simulations[i,1,j], e1=simulations[i,2,j], NIm=NI.marg, sig.level=sig.level.analysis, scale = scale, print.out=FALSE)$CI[2]
          pio0<-simulations[i,1,j]/simulations[i,3,j]
          for (t in 1:length(thresholds)) {
            if ( scale == "RD" ) {
              NI.marg2<-ifelse(abs(p0.expected-pio0)>thresholds[t],
                               sin(asin(sqrt(pio0))+asin(sqrt(p1.tolerable))-asin(sqrt(p0.expected)))^2-pio0,
                               NI.marg)
              simulations[i,5+t,j]<-(abs(simulations[i,5,j])<NI.marg2)
            } else if ( scale =="RR" ) {
              NI.marg2<-ifelse(abs(log(p0.expected/pio0))>thresholds[t],
                               log(sin(asin(sqrt(pio0))+asin(sqrt(p1.tolerable))-asin(sqrt(p0.expected)))^2/pio0),
                               NI.marg)
              simulations[i,5+t,j]<-(simulations[i,5,j]<NI.marg2)
            }
            simulations[i,5+length(thresholds)+t,j]<-(NI.marg!=NI.marg2)
          }
          if (is.na(simulations[i,6,j])) simulations[i,6:(5+length(thresholds)),j]<-0
        }
        for (t in 1:length(thresholds)) {
          Power.method[j,t]<-sum(simulations[,5+t,j])/n.sim
          n.modified[j,t]<-sum(simulations[,5+length(thresholds)+t,j])/n.sim
        }
        if (print.out==T) {
          cat("Scenario ", j, " (",intToUtf8(960),"0=", range.of.p0[j],") completed\n", sep="")
        }
      }
      
    } else {
      simulations<-array(NA, c(n.sim,5+length(sig.level.analysis)*2, n.pi0))
      Power.method<-matrix(NA,n.pi0,length(sig.level.analysis)+1)
      colnames(Power.method)<-c(paste(expression(alpha), "=",sig.level.analysis), expression(pi[0]))
      n.modified<-matrix(NA,n.pi0,length(sig.level.analysis))
      colnames(n.modified)<-paste(expression(alpha), "=",sig.level.analysis)
      sample.size<-ceiling(sample.size.NI( p0.expected = p0.expected, p1.expected = p1.expected, p1.tolerable = p1.tolerable, sig.level = sig.level.design, power = power, r = r, scale = scale, print.out = print.out ))  # Sample size
      
      for (j in 1:n.pi0) {
        pi0<-range.of.p0[j]  # True control event rate
        if (perf.measure == "power") {
          pi1<-pi0*(p1.expected/p0.expected)  # True active event rate
        } else {
          pi1<-sin(asin(sqrt(pi0))+asin(sqrt(p1.tolerable))-asin(sqrt(p0.expected)))^2  # True active event rate
        }
        
        # Data and results matrices
        
        simulations[,3,j]<-sample.size[1]
        simulations[,4,j]<-sample.size[2]
        
        for (i in 1:n.sim) {
          sample.size.obs<-sum(rbinom(1,sample.size[1]+sample.size[2],1/(1+r)))
          sample.size.obs<-c(sample.size.obs, sample.size[1]+sample.size[2]-sample.size.obs)
          simulations[i,1,j]<-rbinom(1,sample.size.obs[1],pi0)
          simulations[i,2,j]<-rbinom(1,sample.size.obs[2],pi1)
          pio0<-simulations[i,1,j]/simulations[i,3,j]
          for (s in 1:length(sig.level.analysis)) {
            simulations[i,5,j]<-test.NI(n0=sample.size.obs[1], n1=sample.size.obs[2], e0=simulations[i,1,j], e1=simulations[i,2,j], NIm=NI.marg, sig.level=sig.level.analysis[s], scale = scale, print.out=FALSE)$CI[2]
            if ( scale == "RD" ) {
              NI.marg2<-ifelse(abs(p0.expected-pio0)>thresholds,
                               sin(asin(sqrt(pio0))+asin(sqrt(p1.tolerable))-asin(sqrt(p0.expected)))^2-pio0,
                               NI.marg)
              simulations[i,5+s,j]<-(abs(simulations[i,5,j])<NI.marg2)
            } else if ( scale =="RR" ) {
              NI.marg2<-ifelse(abs(log(p0.expected/pio0))>thresholds,
                               log(sin(asin(sqrt(pio0))+asin(sqrt(p1.tolerable))-asin(sqrt(p0.expected)))^2/pio0),
                               NI.marg)
              simulations[i,5+s,j]<-(simulations[i,5,j]<NI.marg2)
            }
            simulations[i,5+length(sig.level.analysis)+s,j]<-(NI.marg!=NI.marg2)
          }
          if (is.na(simulations[i,6,j])) simulations[i,6:(5+length(sig.level.analysis)),j]<-0
          
        }
        for (s in 1:length(sig.level.analysis)) {
          Power.method[j,s]<-sum(simulations[,5+s,j])/n.sim
          n.modified[j,s]<-sum(simulations[,5+length(sig.level.analysis)+s,j])/n.sim
        }
        if (print.out==T) {
          cat("Scenario ", j, " (",intToUtf8(960),"0=", range.of.p0[j],") completed\n", sep="")
        }
      }
      
    }
    Power.method[,ncol(Power.method)]<-range.of.p0
    results<-list(Power.method, n.modified, type, scale, perf.measure)
    if (perf.measure == "power") {
      names(results)<-c("Power", "Modif.rates", "type", "scale", "perf.measure")
    } else {
      names(results)<-c("Type1Error", "Modif.rates", "type", "scale", "perf.measure")
    }
    return(results) 
  }