# plot the estimates of survival functions for each method
# add S1.est-S0.est
# add 95% CI

plot.st <- function(y,delta,treat,x,psix=x,t,method = c("KM","IPTW","EL")){
  if (method == "KM"){
    S1 <- km.est(y,delta,treat,1,t)
    S0 <- km.est(y,delta,treat,0,t)
    S1.est <- S1$St
    S0.est <- S0$St
  }else if (method == "IPTW"){
    S1 <- km.est.pi(y,delta,treat,x,1,t)
    S0 <- km.est.pi(y,delta,treat,x,0,t)
    S1.est <- S1$St
    S0.est <- S0$St
  }else if (method == "EL"){
    S1 <- el.est(y,delta,treat,x,psix,1,t,FALSE)
    S0 <- el.est(y,delta,treat,x,psix,0,t,FALSE)
    S1.est <- S1$St
    S0.est <- S0$St
  }

  plot(x=t,y=S1.est, ylim = c(0,1), xlim = c(min(t),max(t)), col = rgb(0,0,1,1/4),
       main = paste0("Estimates of S(t) by ",method),xlab = "t",ylab = "S(t)",
       type = "l",lwd = 3,lty=1)
  lines(x=t,y=S0.est, col = rgb(1,0,0,1/4), lwd = 3, lty=1)
  legend("topright", c("treatment","control"),
         fill = rgb(0:1,0,1:0,1/4), bty = 'n', border = NA)
}
