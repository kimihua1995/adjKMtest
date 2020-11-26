# sub-function for el.est()
#######################
jacob.fun<-function(treat,x,psix,a,tau,d){
  n <- length(treat)
  n1 <- sum(treat)
  n2 <- n-n1
  psi=as.matrix(psix-t(array(rep(a,n),dim=c(d,n))))
  midd1<-(n1+psi%*%tau)^2;  midd11<-(n1+psi%*%tau)
  midd2<-(n2-psi%*%tau)^2;midd22<-(n2-psi%*%tau)
  E=diag(rep(1,d))
  a1<-a2<-a3<-a4<-0
  for (i in 1:n)
  {
    di<- treat[i]
    psii<-as.matrix(psix[i,]-a)
    midd1i<-midd1[i];midd12i<-midd11[i]
    midd2i<-midd2[i];midd23i<-midd22[i]
    a1<-a1+(-E*midd12i+psii%*%t(tau))*di/midd1i
    a2<-a2+(-psii%*%t(psii))*di/midd1i
    a3<-a3+(-E*midd23i-psii%*%t(tau))*(1-di)/midd2i
    a4<-a4+(psii%*%t(psii))*(1-di)/midd2i
  }
  a1<-a1;a2<-a2
  a3<-a3;a4<-a4
  Jacob<-rbind(cbind(a1,a2),cbind(a3,a4))
}

############This is the function for estimating a and tau
###in EL likelihood function
est.a.tau<-function(treat,x,psix,a,tau,d)
{
  n <- length(treat)
  n1 <- sum(treat)
  n2 <- n-n1
  psi<-psix-t(array(rep(a,n),dim=c(d,n)))
  rep_treat<-array(rep(treat,d),dim=c(n,d))
  f1<-apply(rep_treat*psi/matrix(rep((n1+psi%*%tau),d),ncol=d),2,sum)
  f2<-apply((1-rep_treat)*psi/matrix(rep((n2-psi%*%tau),d),ncol=d),2,sum)
  f<-matrix(c(f1,f2),ncol=1)
}

##########This is the estimating function to obtain the estimator
###of a and tau by netwon raphson method
estamtor<-function(treat,x,psix,a,tau,d)
  #tau=tau0;a=a0;d=dd
{
  n <- length(treat)
  n1 <- sum(treat)
  n2 <- n-n1
  para<-rbind(a,tau);
  a0<-as.matrix(para[1:d]);tau0<-as.matrix(para[-(1:d)])
  fvalue<-est.a.tau(treat,x,psix,a0,tau0,d)
  converge<-F
  if(max(abs(fvalue))< 1.0e-06)
  {
    converge <- T
    return(list(a=as.matrix(para[1:d]),tau=as.matrix(para[-(1:d)]),converge=converge))
  }

  dfvalue<-jacob.fun(treat,x,psix,a0,tau0,d)


  for (i in 1:100)
  {
    #print(i)
    para1<-para-ginv(dfvalue)%*%fvalue

    error<-max(abs(para1-para))
    if (error<1.0e-6|max(abs(fvalue))< 1.0e-06)
    {
      converge <- T
      return(list(a=as.matrix(para1[1:d]),tau=as.matrix(para1[-(1:d)]),converge=converge))
    }

    para<-para1
    a0<-as.matrix(para[1:d]);tau0<-as.matrix(para[-(1:d)])
    fvalue<-est.a.tau(treat,x,psix,a0,tau0,d)
    dfvalue<-jacob.fun(treat,x,psix,a0,tau0,d)
  }

}


estimator.pi<-function(y,delta,treat,x,psix)
{
  dd=ncol(psix)
  tau0=array(rep(0,dd),dim=c(dd,1))
  a0=array(rep(0,dd),dim=c(dd,1))

  n <- length(y)
  n1 <- sum(treat)
  n2 <- n-n1
  estnusi=estamtor(treat,x,psix,a0,tau0,dd)   #the estimator for a and tau
  pi=rep(1,n)

  if (estnusi$converge==T&&length(estnusi$converge)!=0)
  {
    a<-estnusi$a;tau<-estnusi$tau
    psi<-psix-t(array(rep(a,n),dim=c(dd,n)))
    pi1<-1/(n1+psi%*%tau)
    pi2<-1/(n2-psi%*%tau)
    pi[treat==1]=pi1[treat==1];pi[treat==0]=pi2[treat==0]
  }
  return(pi)
}

el.est.old <- function(y,delta,pi,t){
  require(MASS)
  n=length(y)
  m=length(t)
  z0=y[delta==1]
  pi.1=pi[delta==1]
  N=length(z0);

  if (N>0)
  {
    II0=(matrix(rep(z0,m),ncol=m)<=matrix(rep(t,N),ncol=m,byrow=T))*matrix(rep(pi.1,m),ncol=m)
    I0=as.numeric(matrix(rep(y,N),ncol=N)>=matrix(rep(z0,n),ncol=N,byrow=T))*matrix(rep(pi,N),ncol=N)
    ssum=matrix(rep(apply(I0,2,sum),m),ncol=m)
    temp=1-II0/(ssum+0.000001);
  }  else
    temp=matrix(1,n,m)

  suvdf=apply(temp,2,prod)
}


##########################################
# estimate S(t) and standard deviation by EL
el.est <- function(y,delta,treat,x,psix_moment=c("first","second"),treat.select,t,
                   get.sd=TRUE,Nboot=500,standardize=FALSE){
  require(MASS)
  if (standardize) {x <- scale(x)}
  if (psix_moment == "first"){
    psix <- x
  }else if (psix_moment == "second"){
    p <- ncol(x)
    combb <- combn(p,2)
    res <- apply(x, 1, function(x) { apply(combb, 2, function(y) prod(x[y])) })
    psix <- cbind(x^2, t(res))
  }

  Y1=y;Delta1=delta;X1=x;Treat1=treat;Psix1=psix;  # save for bootstrap

  pi<-estimator.pi(y,delta,treat,x,psix)

  y=y[treat==treat.select]
  delta=delta[treat==treat.select]
  pi=pi[treat==treat.select]

  n=length(y)
  m=length(t)
  z0=y[delta==1]
  pi.1=pi[delta==1]
  N=length(z0);

  if (N>0)
  {
    II0=(matrix(rep(z0,m),ncol=m)<=matrix(rep(t,N),ncol=m,byrow=T))*matrix(rep(pi.1,m),ncol=m)
    I0=as.numeric(matrix(rep(y,N),ncol=N)>=matrix(rep(z0,n),ncol=N,byrow=T))*matrix(rep(pi,N),ncol=N)
    ssum=matrix(rep(apply(I0,2,sum),m),ncol=m)
    temp=1-II0/(ssum+0.000001);
  }  else
    temp=matrix(1,n,m)

  suvdf=apply(temp,2,prod)

  if (get.sd){
  s1=NULL;s0=NULL

  n=length(Y1)
  for (i in 1:Nboot)
  {
    U=ceiling(n*runif(n,min=0,max=1))
    y=Y1[U];delta=Delta1[U];x=X1[U,];treat=Treat1[U]
    psix=as.matrix(Psix1[U,])
    n1=sum(treat);n0=n-n1

    y1=y[treat==1];delta1=delta[treat==1]
    y0=y[treat==0];delta0=delta[treat==0]


    pi<-estimator.pi(y,delta,treat,x,psix)
    aa=range(pi)[1];aa1=range(pi)[2]
    if(aa<=0||aa1>=1) next
    pi1<-pi[treat==1];
    s1<-rbind(s1,el.est.old(y1,delta1,pi1,t));

    pi0<-pi[treat==0];
    s0<-rbind(s0,el.est.old(y0,delta0,pi0,t))
  }
  s1.sd=apply(s1,2,sd); s0.sd=apply(s0,2,sd);
  return(list(St = suvdf, sd = treat.select*s1.sd + (1-treat.select)*s0.sd))
  } else {
    return(list(St = suvdf, sd = NA))
    }


}


