rm(list=ls())
library(randtoolbox)
library(MASS)

############# QUESTION 1 ###############
#on doit simuler plusieurs fois la grandeur dans l'espérance
#le seul élément stochastique est le mouvement Brownien
# pour les variables antithétiques: le mouvement brownien est très stable



#simulation du MB
k<-20 #longueur du MB
n<-100 #nombre de simulations
W <- mvrnorm(n=n,mu=rep(0,k),Sigma=diag(k)/k)
W<- apply(W,1,cumsum)

plot(W[,1],type='l',xlim = c(0,20),ylim = c(-3,3))
for (i in 2:n){
  lines(W[,i])
}

#simulation de la valeur qui nous intéresse
value<-function(T=1,r=0.05,k=20,K=0.5,sigma=0.3,mu=0){
  t <- (1:k)/k
  W <- rnorm(n=k-1, mean=0,sd=1/k)
  W<- c(0,cumsum(W))
  S <- exp((mu-sigma**2/2)*t)
  S<-S*exp(sigma*W)
  plot(S,type='l')
  return(exp(-r*T)*max((mean(S)-K),0))
}


#plusieurs simulations
v<-c()
for (i in 1:100){
  v<-c(v,value())
}

plot(v,type='l')
