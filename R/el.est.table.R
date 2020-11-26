boot.sd<-function(Y1,Delta1,X1,Treat1,psix_moment=c("first","second"),
                  Nboot,wt,standardize=FALSE,alpha=0.05)
{#Y1=y;Delta1=delta;X1=x;Treat1=treat;Psix1=psix;

  require(MASS)
  if (standardize) {X1 <- scale(X1)}
  if (psix_moment == "first"){
    Psix1 <- X1
  }else if (psix_moment == "second"){
    p <- ncol(X1)
    combb <- combn(p,2)
    res <- apply(X1, 1, function(x) { apply(combb, 2, function(y) prod(x[y])) })
    Psix1 <- cbind(X1, X1^2, t(res))
  }


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
    s1<-rbind(s1,el.est.old(y1,delta1,pi1,wt));

    pi0<-pi[treat==0];
    s0<-rbind(s0,el.est.old(y0,delta0,pi0,wt))
  }
  s1.sd=apply(s1,2,sd); s0.sd=apply(s0,2,sd);
  Dif=s1-s0; Dif.sd=apply(Dif,2,sd)
  s1.CI=apply(s1,2,quantile,probs=c(alpha/2, 1-alpha/2))
  s0.CI=apply(s0,2,quantile,probs=c(alpha/2, 1-alpha/2))
  Dif.CI=apply(Dif,2,quantile,probs=c(alpha/2, 1-alpha/2))
  return(list(sd=Dif.sd, s1.CI=s1.CI, s0.CI=s0.CI, Dif.CI=Dif.CI))
}



# generate the table with S(t) with sd at selected time points by EL
#    for treatment group, control group, and their difference

el.est.table <- function(y,delta,treat,x,psix_moment=c("first","second"),t,
                         Nboot=500,standardize=FALSE){

  S1 <- el.est(y,delta,treat,x,psix_moment,1,t,TRUE,Nboot,standardize)
  S0 <- el.est(y,delta,treat,x,psix_moment,0,t,TRUE,Nboot,standardize)

  S1.est <- S1$St; S1.sd <- S1$sd
  S0.est <- S0$St; S0.sd <- S0$sd

  Dif.est <- S1.est - S0.est; Dif.sd <- boot.sd(y,delta,x,treat,psix_moment,Nboot,t,standardize)$sd

  table.est <- round(rbind(S1.est,S0.est,Dif.est),3)
  table.sd <- round(rbind(S1.sd,S0.sd,Dif.sd),3)


  require(formattable)

  table <- matrix(NA,nrow = length(t), ncol = 4)
  colnames(table) <- c("Method:EL","treatment","control","Difference")
  for (i in 1:length(t)){
    table[i,] <- c(paste0("t=",t[i]),
                   paste(format(table.est[,i],nsmall = 3),"(",
                         format(table.sd[,i],nsmall = 3),")"))
  }
  table <- as.data.frame(table)



  return(formattable(table))

}
