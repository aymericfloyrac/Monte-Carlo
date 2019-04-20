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
value<-function(T=1,r=0.05,k=20,K=5,sigma=0.3,mu){
  if (missing(mu)){
    mu<-r
    } #on suppose que, si non spécifié, le drift et le taux d'intérêt sont égaux
  t <- (1:k)/k
  W <- rnorm(n=k-1, mean=0,sd=1/k)
  W<- c(0,cumsum(W))
  S <- K*exp((mu-sigma**2/2)*t) #est-ce qu'on suppose que la valeur du BS au temps 0 vaut K? 
  S<-S*exp(sigma*W)
  return(exp(-r*T)*max(mean(S)-K,0))
}
value()

#plusieurs simulations
v<-c()
for (i in 1:1000){
  v<-c(v,value())
}

hist(v,type='l')
sum(v>0)
sd(v)

#variables antithétiques 
# concernant le MB, on sait que -W a même loi que W. 

antithetic<-function(T=1,r=0.05,k=20,K=5,sigma=0.3,mu){
  if (missing(mu)){
    mu<-r
  } #on suppose que, si non spécifié, le drift et le taux d'intérêt sont égaux
  t <- (1:k)/k
  W <- rnorm(n=k-1, mean=0,sd=1/k)
  W1 <- c(0,cumsum(W))
  W2 <- c(0,-cumsum(W))
  S <- K*exp((mu-sigma**2/2)*t) #est-ce qu'on suppose que la valeur du BS au temps 0 vaut K? 
  S1<-S*exp(sigma*W1)
  S2<-S*exp(sigma*W2)
  return(0.5*exp(-r*T)*(max(mean(S1)-K,0)+max(mean(S2)-K,0)))
}

antithetic()
v<-c()
for (i in 1:1000){
  v<-c(v,antithetic())
}
sd(v) #on observe bien une réduction de variance environ d'un facteur 2

#variables de contrôle
# idée (à voir): prendre le Black Scholes comme variable de contrôle


