cdfnorm<-function(x) integrate(dnorm,-Inf,x)$value
UOCTheo<-function(S0,K,H,Tm,mu,sigma,r){
  lamda=mu/sigma-sigma/2
  lamda_=mu/sigma+sigma/2
  h=1/sigma*log(H/S0)
  k=1/sigma*log(K/S0)
  alpha=exp(2*lamda*h)
  alpha_=exp(2*lamda_*h)
  
  d1= lamda_*sqrt(Tm)-k/sqrt(Tm)
  d2= lamda *sqrt(Tm)-k/sqrt(Tm)
  d3= lamda_*sqrt(Tm)-h/sqrt(Tm)
  d4= lamda *sqrt(Tm)-h/sqrt(Tm)
  d5=-lamda *sqrt(Tm)-h/sqrt(Tm)
  d6=-lamda_*sqrt(Tm)-h/sqrt(Tm)
  d7=-lamda *sqrt(Tm)-(2*h-k)/sqrt(Tm)
  d8=-lamda_*sqrt(Tm)-(2*h-k)/sqrt(Tm)
  UOC<-function(S,K,H) exp((mu-r)*Tm)*S*(cdfnorm(d1)-cdfnorm(d3)-alpha_*(cdfnorm(d6)-cdfnorm(d8)))-exp(-r*Tm)*K*(cdfnorm(d2)-cdfnorm(d4)-alpha*(cdfnorm(d5)-cdfnorm(d7)))
  UOC(S0,K,H)
}

max_K<-function(S) max(S-K,0)
cBS<-function(t,x) {
  dPlus<-function(t,x) 1/(sigma*sqrt(t))*(log(x/K)+(r+sigma^2/2)*t)
  dMinus<-function(t,x) 1/(sigma*sqrt(t))*(log(x/K)+(r-sigma^2/2)*t)
  if (t==Tm){
    0
  }else{x*cdfnorm(dPlus(Tm-t,x))-K*exp(-r*(Tm-t))*cdfnorm(dMinus(Tm-t,x))}
}

library(quadprog)

K=100
S0=100
H=120
Tm=1
mu=1*0.01
sigma=20*0.01
r=0

pTheo<-UOCTheo(S0,K,H,Tm,mu,sigma,r)
cBSPrice0<-cBS(0,K)

nB = 20

nSteps=10000
nTimes=5000
dt=Tm/nSteps
dtSd=sqrt(dt)

MCPriceA<-rep(0,nB)
MCPriceB<-rep(0,nB)
MCPriceC<-rep(0,nB)
MCPriceD<-rep(0,nB)
timeTable<-matrix(rep(0,4*nB),nrow=nB)
ST<-rep(0,nTimes)
px<-rep(0,nTimes)
for(t in 1:nB) {
  W=rnorm(nSteps*nTimes)
  # a) basic Monte Carlo
  starttime=proc.time()
  
  X<-matrix(mu*dt+sigma*dtSd*W,nrow = nTimes)
  for(j in 1:nTimes){
    S=S0
    b=FALSE
    for (i in 1:nSteps){
      S=S*(1+X[j,i])
      if (S>H) 
      {b=TRUE}
    }
    ST[j]<-S
    if (b) {p=0} else {p=max(S-K,0)*exp(-r*Tm)}
    px[j]<-p
  }
  MCPriceA[t]<-mean(px)
  timeTable[t,1]<- (proc.time()-starttime)[3]
  
  # c) rescale by weight
  starttime=proc.time()
  Dmat=diag(1,nTimes)
  dvec=rep(0,nTimes)
  A<-matrix(c(rep(1,nTimes),sapply(ST, max_K)*exp(-r*Tm),diag(1,nTimes)),ncol=2+nTimes)
  bvec=c(1,cBSPrice0,rep(0,nTimes))
  w<-solve.QP(Dmat=Dmat, dvec=dvec, Amat=A, bvec=bvec, meq=2)$solution
  MCPriceC[t]<-px%*%w
  timeTable[t,3]<- (proc.time()-starttime)[3]+timeTable[t,1]
  
  # b) rescale by value
  starttime=proc.time()

  sdW=sd(W)
  meanW=mean(W)
  Z=(W-meanW)/sdW
  X<-matrix( mu*dt+sigma*dtSd*Z, nrow = nTimes )
  for (j in 1:nTimes){
    S=S0
    b=FALSE
    for (i in 1:nSteps){
      S=S*(1+X[j,i])
      if (S>H) 
      {b=TRUE}
    }
    ST[j]<-S
    if (b) {p=0} else {p=max(S-K,0)*exp(-r*Tm)}
    px[j]<-p
  }
  MCPriceB[t]<-mean(px)
  timeTable[t,2]<- (proc.time()-starttime)[3]

  
  # c) rescale by value and weight
  starttime=proc.time()
  Dmat=diag(1,nTimes)
  dvec=rep(0,nTimes)
  A<-matrix(c(rep(1,nTimes),sapply(ST, max_K)*exp(-r*Tm),diag(1,nTimes)),ncol=2+nTimes)
  bvec=c(1,cBSPrice0,rep(0,nTimes))
  w<-solve.QP(Dmat=Dmat, dvec=dvec, Amat=A, bvec=bvec, meq=2)$solution
  MCPriceD[t]<-px%*%w
  timeTable[t,4]<- (proc.time()-starttime)[3]+timeTable[t,2]
  
  print(timeTable[t,])
}
MC_average<-c(mean(MCPriceA),mean(MCPriceB),mean(MCPriceC),mean(MCPriceD))
MC_Diff<-c(mean(MCPriceA),mean(MCPriceB),mean(MCPriceC),mean(MCPriceD))-pTheo
MC_Sd<-c(sd(MCPriceA),sd(MCPriceB),sd(MCPriceC),sd(MCPriceD))
print(c(sum(timeTable[,1]),sum(timeTable[,2]),sum(timeTable[,3]),sum(timeTable[,4])))
