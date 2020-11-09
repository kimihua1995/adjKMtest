boot.sd<-function(Y1,Delta1,X1,Treat1,Psix1,Nboot,wt)
{#Y1=y;Delta1=delta;X1=x;Treat1=treat;Psix1=psix;


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
  s1.CI=apply(s1,2,quantile,probs=c(.025,.975))
  s0.CI=apply(s0,2,quantile,probs=c(.025,.975))
  Dif.CI=apply(Dif,2,quantile,probs=c(.025,.975))
  return(list(sd=Dif.sd, s1.CI=s1.CI, s0.CI=s0.CI, Dif.CI=Dif.CI))
}



# generate the table with S(t) with sd at selected time points by EL
#    for treatment group, control group, and their difference

el.est.table <- function(y,delta,treat,x,psix,t,Nboot=500){

  S1 <- el.est(y,delta,treat,x,psix,1,t)
  S0 <- el.est(y,delta,treat,x,psix,0,t)

  S1.est <- S1$St; S1.sd <- S1$sd
  S0.est <- S0$St; S0.sd <- S0$sd

  Dif.est <- S1.est - S0.est; Dif.sd <- boot.sd(y,delta,x,treat,psix,Nboot,t)$sd

  table.est <- round(rbind(S1.est,S0.est,Dif.est),3)
  table.sd <- round(rbind(S1.sd,S0.sd,Dif.sd),3)


  require(formattable)

  table <- matrix(NA,nrow = 3,ncol = length(t)+2)
  table[,1] <- rep("EL",3)
  table[,2] <- c("treatment","control","Difference")
  colnames(table) <- c("Method","Parameter",paste0("t=",t))
  for (i in 1:length(t)){
    table[,i+2] <- paste(table.est[,i],"(",table.sd[,i],")")
  }
  table <- as.data.frame(table)



  return(formattable(table))

}
