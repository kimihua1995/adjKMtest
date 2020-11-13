# plot the estimates of survival functions for each method
# add S1.est-S0.est
# add 95% CI
# add shade
# choose alpha level

plot.st <- function(y,delta,treat,x,psix=x,t,method = c("KM","IPTW","EL"),Nboot=500){
  if (method == "KM"){
    S1 <- km.est(y,delta,treat,1,t)
    S0 <- km.est(y,delta,treat,0,t)
    S1.est <- S1$St
    S0.est <- S0$St
    Dif.est <- S1.est - S0.est
    S1.sd <- S1$sd
    S0.sd <- S0$sd
    Dif.sd <- sqrt(S1.sd^2 + S0.sd^2)
  }else if (method == "IPTW"){
    S1 <- km.est.pi(y,delta,treat,x,1,t)
    S0 <- km.est.pi(y,delta,treat,x,0,t)
    S1.est <- S1$St
    S0.est <- S0$St
    Dif.est <- S1.est - S0.est
    S1.sd <- S1$sd
    S0.sd <- S0$sd
    Dif.sd <- sqrt(S1.sd^2 + S0.sd^2)
  }else if (method == "EL"){
    S1 <- el.est(y,delta,treat,x,psix,1,t,FALSE)
    S0 <- el.est(y,delta,treat,x,psix,0,t,FALSE)
    S1.est <- S1$St
    S0.est <- S0$St
    Dif.est <- S1.est - S0.est
    boot <- boot.sd(y,delta,x,treat,psix,Nboot,t)
    S1.CI <- boot$s1.CI
    S0.CI <- boot$s0.CI
    Dif.CI <- boot$Dif.CI
  }

  if (method == "EL"){
  layout(matrix(c(1,1,2,2),nrow=1,ncol=4,byrow=T))
  plot(x=t,y=S1.est, ylim = c(0,1), xlim = c(min(t),max(t)), col = rgb(0,0,1,1/2),
       main = paste0("Estimates of S(t) and ", expression(Delta(t))," by ",method),
       xlab = "t",ylab = "S(t)",
       type = "l",lwd = 3,lty=1,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
  lines(x=t,y=S1.CI[1,], lwd=2, lty=3, col = "light blue")
  lines(x=t,y=S1.CI[2,], lwd=2, lty=3, col = "light blue")
  lines(x=t,y=S0.est, col = rgb(1,0,0,1/2), lwd = 3, lty=1)
  lines(x=t,y=S0.CI[1,], lwd=2, lty=3, col = "light coral")
  lines(x=t,y=S0.CI[2,], lwd=2, lty=3, col = "light coral")
  legend("bottomleft", c("treatment","control"),
         fill = rgb(0:1,0,1:0,1/2), bty = 'n', border = NA,cex=1.5)

  plot(x=t,y=Dif.est,xlim = c(min(t),max(t)),ylim = c(-0.5,0.5),col = "orange",
       ylab = expression(Delta(t)), xlab = "t",
       type = "l",lwd = 2,lty=1,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
  lines(x=t,y=Dif.CI[1,], lwd=2, lty=3, col = "light salmon")
  lines(x=t,y=Dif.CI[2,], lwd=2, lty=3, col = "light salmon")
  abline(h=0, lty=2,lwd=2)
  }else{
    layout(matrix(c(1,1,2,2),nrow=1,ncol=4,byrow=T))
    plot(x=t,y=S1.est, ylim = c(0,1), xlim = c(min(t),max(t)), col = rgb(0,0,1,1/2),
         main = paste0("Estimates of S(t) and ", expression(Delta(t))," by ",method),
         xlab = "t",ylab = "S(t)",
         type = "l",lwd = 3,lty=1,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
    lines(x=t,y=S1.est+qnorm(.975)*S1.sd, lwd=2, lty=3, col = "light blue")
    lines(x=t,y=S1.est-qnorm(.975)*S1.sd, lwd=2, lty=3, col = "light blue")
    lines(x=t,y=S0.est, col = rgb(1,0,0,1/2), lwd = 3, lty=1)
    lines(x=t,y=S0.est+qnorm(.975)*S0.sd, lwd=2, lty=3, col = "light coral")
    lines(x=t,y=S0.est-qnorm(.975)*S0.sd, lwd=2, lty=3, col = "light coral")
    legend("bottomleft", c("treatment","control"),
           fill = rgb(0:1,0,1:0,1/2), bty = 'n', border = NA,cex=1.5)

    plot(x=t,y=Dif.est,xlim = c(min(t),max(t)),ylim = c(-0.5,0.5),col = "orange",
         ylab = expression(Delta(t)), xlab = "t",
         type = "l",lwd = 2,lty=1,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
    lines(x=t,y=Dif.est+qnorm(.975)*Dif.sd, lwd=2, lty=3, col = "light salmon")
    lines(x=t,y=Dif.est-qnorm(.975)*Dif.sd, lwd=2, lty=3, col = "light salmon")
    abline(h=0, lty=2,lwd=2)
  }
}
