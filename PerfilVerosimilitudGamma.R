library(maxLik)
x= c(2.3, 0.7, 8.0, 4.4, 15.6, 3.0, 1.9, 1.5, 1.7, 1.0, 2.4, 0.5, 3.3, 5.7, 6.9, 8.9, 1.9, 9.0, 4.5, 3.1)

loglik.gamma = function(parm,x){
  sum(dgamma(x,parm[1],1/parm[2],log=T))
}
loglik.gamma.beta = function(parm,x,alfa){
  sum(dgamma(x,parm,1/alfa,log=T))
}

mle.gamma = maxLik(loglik.gamma,start = c(1,1),x=x,method='NR')

alfa0.1 = seq(0.5,10,length.out=500)
beta0.1 = NULL
ll0.1 = NULL
for(i in 1:length(alfa0.1)){
  if(i ==1){beta00 =2}else{beta00 = beta0.1[i-1]}
  mleEst0 = maxLik(loglik.gamma.beta,x=x,alfa=alfa0.1[i],start=beta00)
  beta0.1[i] = mleEst0$estimate
  ll0.1[i] = mleEst0$maximum
}

alp = seq(0.5,6,length.out=100)
bet = seq(0.5,6,length.out=100)
griddf <- expand.grid(alp,bet)
llr = apply(griddf,1,function(y){loglik.gamma(y,x)})
llrm = matrix(llr,nrow=100,byrow = F)

filled.contour(bet,alp,llrm,nlevels=25,
               plot.axes={
                 lines(beta0.1,alfa0.1,col=1,type='l')
                 points(mle.gamma$estimate[1],mle.gamma$estimate[2],pch=19)
                 axis(side=1)
                 axis(side=2)
               },xlab=expression(theta[2]),ylab=expression(theta[1])              
)

val =  2*(logLik(mle.gamma)-ll0.1)
cond = val < qchisq(0.95,1)
IC = range(beta0.1[cond])
IC
confint(mle.gamma)
yl = logLik(mle.gamma)-qchisq(0.95,1)/2

plot(beta0.1,ll0.1,type='l',xlim=c(min(beta0.1),3),ylim=c(-52.8,max(ll0.1)),
     xlab=expression(theta[1]),ylab='log-versomilitud')
abline(h=yl,lty=2)
abline(v=IC,lty=2)
abline(v=mle.gamma$estimate[1])

#############
