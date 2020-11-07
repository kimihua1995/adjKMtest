# calculating propensity score
ps <- function(treat,x){
  glm.sol1<-glm((treat)~x, family=binomial(link="logit"))
  beta=glm.sol1$coef
  pre<-cbind(1,x.new)%*%beta
  ex<-(1-1/(1+exp(pre)))

  return(ex)
}
