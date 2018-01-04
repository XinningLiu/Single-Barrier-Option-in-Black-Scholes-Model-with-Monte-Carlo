K=100
S0=100
H=120
Tm=1
mu=1*0.01
sigma=20*0.01
r=0

nSteps=1000
nTimes=30
dt=Tm/nSteps
dtSd=sqrt(dt)


plot(rep(H,nSteps+1),type="l",col="red",main="Basic Monte Carlo simulation")
W=rnorm(nSteps*nTimes)
px<-c()
X<-matrix(mu*dt+sigma*dtSd*W,nrow = nTimes)
for (j in 1:nTimes){
  S<-c(S0)
  b=FALSE
  for (i in 1:nSteps){
    S[i+1]<-S[i]*(1+X[j,i])
    if (S[i+1]>H) {
      b=TRUE
    }
  }
  if (b) {p=0} else {p=max(S[nSteps+1]-K,0)*exp(-r*Tm)}
  lines(S,type="l",col=rgb(runif(1),runif(1),runif(1)))
  px[j]<-p
}
mean(px)
sd(px)/sqrt(nTimes)