# estimated S(t) and its standard deviation by IPTW

km.est.pi<-function(y,delta,treat,x,treat.select,t)
  #  y: observed time;
  #  delta: indicator for event
  #  treat: indicator for treatment
  #  x: covariance matrix
  #  treat.select: 1 or 0
  #  t: time points want to estimate (vector)
  #  return a list with S(t) and sd
{
  n=length(y)
  m=length(t)
  z0=y[delta==1]
  ex=ps(treat,x)
  pi=(treat.select == 1)*treat/ex + (treat.select == 0)*(1-treat)/(1-ex)
  # pi: treatment weight by propensity score: treat/ex or (1-treat)/(1-ex)
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

  suvdf=apply(temp,2,prod)  # point estimate


  y1=matrix(rep(y,n),ncol=n)
  y2=matrix(rep(y,n),ncol=n,byrow=T)

  if (treat.select == 1){
    mid1=apply(matrix(rep(treat/ex^2,n),ncol=n)*(y1>=y2),2,sum)/n
    yw1=apply(matrix(rep(treat/ex,n),ncol=n)*(y1>=y2),2,sum)/n
    s=sapply(t, function(t) sum(treat/(ex)*delta*(y<=t)*(yw1+0.000001)^(-3)*mid1)/n)
  }else if (treat.select == 0){
    mid1=apply(matrix(rep((1-treat)/(1-ex)^2,n),ncol=n)*(y1>=y2),2,sum)/n
    yw1=apply(matrix(rep((1-treat)/(1-ex),n),ncol=n)*(y1>=y2),2,sum)/n
    s=sapply(t, function(t) sum((1-treat)/(1-ex)*delta*(y<=t)*(yw1+0.000001)^(-3)*mid1)/n)
  }

  sd.s1=sqrt(s/n)*suvdf

  return(list(St = suvdf, sd = sd.s1))
}
