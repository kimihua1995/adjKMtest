# generate the table with S(t) with sd at selected time points by standard KM methods
#    for treatment group, control group, and their difference

km.est.table <- function(y,delta,treat,t){

  S1 <- km.est(y,delta,treat,1,t)
  S0 <- km.est(y,delta,treat,0,t)

  S1.est <- S1$St; S1.sd <- S1$sd
  S0.est <- S0$St; S0.sd <- S0$sd
  Dif.est <- S1.est - S0.est; Dif.sd <- sqrt(S1.sd^2 + S0.sd^2)

  table.est <- round(rbind(S1.est,S0.est,Dif.est),3)
  table.sd <- round(rbind(S1.sd,S0.sd,Dif.sd),3)

  require(formattable)

  table <- matrix(NA,nrow = 3,ncol = length(t)+2)
  table[,1] <- rep("KM",3)
  table[,2] <- c("treatment","control","Difference")
  colnames(table) <- c("Method","Parameter",paste0("t=",t))
  for (i in 1:length(t)){
    table[,i+2] <- paste(table.est[,i],"(",table.sd[,i],")")
  }
  table <- as.data.frame(table)



  return(formattable(table))
}
