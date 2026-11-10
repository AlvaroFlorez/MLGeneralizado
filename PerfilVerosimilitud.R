library(faraway)
library(maxLik)
logdose <- c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839)
dead <- c(6, 13, 18, 28, 52, 53, 61, 60) # numbers dead
n <- c(59, 60, 62, 56, 63, 59, 62, 60) # binomial sample sizes
Datos=data.frame(logdose,n,dead)
DatosEsc = Datos
colnames(DatosEsc) = c('log. dosis','expuestos', 'muertos')

### funcion verosimilitud
### completa
llBinfull = function(beta,yval,n,xval){
  eta = beta[1] + xval*beta[2]
  prob = exp(eta)/(1+exp(eta))
  ll = sum(dbinom(yval,n,prob,log=T))
  ll
}
## perfil para beta1
llBin = function(beta0,beta1,yval,n,xval){
  eta = beta0 + xval*beta1
  prob = exp(eta)/(1+exp(eta))
  ll = sum(dbinom(yval,n,prob,log=T))
  ll
}



mle = glm(cbind(dead,n-dead)~logdose,family=binomial,data=DatosEsc)
summary(mle)
betaest = mle$coefficients

### grilla
beta00 = seq(-65,-55, length.out=100)
beta10 = seq(25,45, length.out=100)
ll00 = matrix(NA,length(beta00),length(beta10))
for(i in 1:length(beta00)){
  for(j in 1:length(beta10)){
    ll00[i,j] = llBinfull(c(beta00[i],beta10[j]),yval=dead,n=n,xval=logdose)
  }
}
max(ll00)

llBinfull(c(-70,40),yval=dead,n=n,xval=logdose)

filled.contour(beta00, beta10, ll00,nlevels = 100)

beta10.1 = seq(24,46,length.out=500)
beta00.1 = NULL
ll0.1 = NULL
for(i in 1:length(beta10.1)){
  mleEst0 = maxLik(llBin,beta1=beta10.1[i],yval=dead,n=n,xval=logdose,start=-60)
  beta00.1[i] = mleEst0$estimate
  ll0.1[i] = mleEst0$maximum
}
val =  2*(logLik(mle)-ll0.1)
cond = val < qchisq(0.95,1)
IC = range(beta10.1[cond])
IC
yl = logLik(mle)-qchisq(0.95,1)/2



filled.contour(beta00, beta10, ll00,nlevels = 50,
               plot.axes={
                 lines(beta00.1,beta10.1)
                 axis(side=1)
                 axis(side=2)
                 points(betaest[1],betaest[2],pch=20)
               },xlab=expression(beta[0]),ylab=expression(beta[1]))

plot(beta10.1,ll0.1,type='l',xlab=expression(beta[1]),
     ylab='log verosimilitud')
abline(v=betaest[2],col=2)
abline(v=IC,col=2,lty=2)
abline(h = yl,lty=2)



