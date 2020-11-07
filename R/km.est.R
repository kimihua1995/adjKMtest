# estimate S(t) and its standard deviation by standard KM method

km.est<-function(y,delta,treat,treat.select,t)
  #  y: observed time;
  #  delta: indicator for event
  #  treat: indicator for treatment
  #  treat.select: 1 or 0
  #  t: time points want to estimate (vector)
  #  return a list with S(t) and sd
{
  y=y[treat==treat.select];delta=delta[treat==treat.select]

  n=length(y)
  m=length(t)
  z0=y[delta==1]
  N=length(z0);

  if (N>0){
    II0=(matrix(rep(z0,m),ncol=m)<=matrix(rep(t,N),ncol=m,byrow=T))
    I0=(matrix(rep(y,N),ncol=N)>=matrix(rep(z0,n),ncol=N,byrow=T))
    ssum=matrix(rep(apply(I0,2,sum),m),ncol=m)
    temp=1-II0/(ssum+0.000001);
  }else
    temp=matrix(1,n,m)

  suvdf=apply(temp,2,prod)  # point estimate

  y1=matrix(rep(y,n),ncol=n)
  y2=matrix(rep(y,n),ncol=n,byrow=T)
  ybar=apply((y1>=y2),2,sum)
  vt=sapply(t,function(t) sum(delta*(y<=t)/ybar/(ybar)))
  sd=sqrt(vt)*suvdf        # standard deviation

  return(list(St = suvdf, sd = sd))

}
