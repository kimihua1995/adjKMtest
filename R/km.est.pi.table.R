# generate the table with S(t) with sd at selected time points by IPTW
#    for treatment group, control group, and their difference

km.est.pi.table <- function(y,delta,treat,x,t,standardize=FALSE){

  S1 <- km.est.pi(y,delta,treat,x,1,t,standardize)
  S0 <- km.est.pi(y,delta,treat,x,0,t,standardize)

  S1.est <- S1$St; S1.sd <- S1$sd
  S0.est <- S0$St; S0.sd <- S0$sd
  Dif.est <- S1.est - S0.est; Dif.sd <- sqrt(S1.sd^2 + S0.sd^2)

  table.est <- round(rbind(S1.est,S0.est,Dif.est),3)
  table.sd <- round(rbind(S1.sd,S0.sd,Dif.sd),3)


  require(formattable)

  table <- matrix(NA,nrow = length(t), ncol = 4)
  colnames(table) <- c("Method:IPTW","treatment","control","Difference")
  for (i in 1:length(t)){
    table[i,] <- c(paste0("t=",t[i]),
                   paste(format(table.est[,i],nsmall = 3),"(",
                         format(table.sd[,i],nsmall = 3),")"))
  }
  table <- as.data.frame(table)



  return(formattable(table))
}
