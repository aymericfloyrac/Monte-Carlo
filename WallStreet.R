rm(list=ls())
library(MASS)

"
une idée générale: pour le rendu graphique, on peut plotter l'évolution de la moyenne
à chaque période (k périodes) avec des error bars
"

########################################
############# QUESTION 1 ###############
########################################

#on doit simuler plusieurs fois la grandeur dans l'espérance
#le seul élément stochastique est le mouvement Brownien
# pour les variables antithétiques: le mouvement brownien est très stable



############# SIMULATION DU MB ###############
#dans un premier temps on fait une construction random walk 
k<-20 #longueur du MB
n<-1000 #nombre de simulations
W <- rnorm(n=n*k,mean = 0,sd=1/k)
W<-matrix(W,nrow=n,ncol=k)
W<- t(apply(W,1,cumsum))
dim(W)

plot(W[1,],type='l',xlim = c(0,20),ylim = c(-3,3))
for (i in 2:n){
  lines(W[i,])
}

plot(colMeans(W),type='l')

#dans un deuxieme temps on peut faire une construction brownian bridge 


############# MC STANDARD ###############

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

mean(v)
sd(v)

############# ANTITHETIQUE ###############
# concernant le MB, on sait que -W a même loi que W. 

monteCarlo<-function(N=1000,T=1,r=0.05,k=20,K=5,sigma=0.3,mu,antithet=T){
  if (missing(mu)){
    mu<-r
  } #on suppose que, si non spécifié, le drift et le taux d'intérêt sont égaux
  t <- (1:k)/k
  W <- rnorm(n=N*(k-1), mean=0,sd=1/k)
  W <- matrix(W,nrow = N,ncol = k-1)
  W<-t(apply(W,1,cumsum))
  W1 <- c(rep(0,N),W)
  W1<-matrix(W1,nrow = N,ncol = k)
  S <- rep(K*exp((mu-sigma**2/2)*t),N) 
  S <-matrix(S,nrow=N,ncol=k)
  S1 <-S*exp(sigma*W1)
  if (antithet){
    W2 <- -W1
    S2 <-S*exp(sigma*W2)
    C <- 0.5*exp(-r*T)*(pmax(rowMeans(S1)-K,0)+pmax(rowMeans(S2)-K,0))
  }
  else{
    C <- exp(-r*T)*pmax(rowMeans(S1)-K,0)
  }
  #beta <- cov(W1[,k],C)/var(W1[,k])
  #C <- C-beta*W1[,k]
  #return (list('mean'=mean(C),'sd'=sd(C)))
  return (c(mean(C),sd(C)))
}

monteCarlo(antithet=F)
monteCarlo(antithet = T)

############# VARIABLE DE CONTROLE ###############
# idée (à voir): prendre le Black Scholes comme variable de contrôle, ou bien juste le BM 

#on change la façon de simuler parce qu'on va avoir besoin 
#d'estimateurs de cov et de variances

control<-function(N=1000,T=1,r=0.05,k=20,K=5,sigma=0.3,mu){
  if (missing(mu)){
    mu<-r
  } #on suppose que, si non spécifié, le drift et le taux d'intérêt sont égaux
  t <- (1:k)/k
  W <- rnorm(n=N*(k-1), mean=0,sd=1/k)
  W <- matrix(W,nrow = N,ncol = k-1)
  W<-t(apply(W,1,cumsum))
  W1 <- c(rep(0,N),W)
  W1<-matrix(W1,nrow = N,ncol = k)
  W2 <- -W1
  S <- rep(K*exp((mu-sigma**2/2)*t),N) 
  S <-matrix(S,nrow=N,ncol=k)
  S1 <-S*exp(sigma*W1)
  S2 <-S*exp(sigma*W2)
  C <- 0.5*exp(-r*T)*(pmax(rowMeans(S1)-K,0)+pmax(rowMeans(S2)-K,0))
  beta <- cov(W1[,k],C)/var(W1[,k])
  C <- C-beta*W1[,k]
  #return (list('mean'=mean(C),'sd'=sd(C)))
  return (c(mean(C),sd(C)))
}

control()



############# QUASI MC ###############
#si on devait coder un générateur QMC nous mêmes on pourrait faire du Box Muller
# avec des uniformes en treillis
#mais on va utiliser la librairie déjà faite
library(randtoolbox)
# on va utiliser les fonctions halton, torus et sobol pour générer des gaussiennes
# pour faire le RQMC on utilise l'option scrambling de la fonction sobol

qmc<-function(N=1000,T=1,r=0.05,k=20,K=5,sigma=0.3,mu,funcStr='halton',antithet=T,randomized=F){
  if (missing(mu)){
    mu<-r
  } #on suppose que, si non spécifié, le drift et le taux d'intérêt sont égaux
  t <- (1:k)/k
  if (randomized){
    W <- (1/k)*sobol(normal=T,n = N,dim = k-1,scrambling = 3)
  }
  else{
    switch (funcStr,
      halton={
        W <- (1/k)*halton(normal=T,n = N,dim = k-1)
      },
      torus={
        W <- (1/k)*torus(normal=T,n = N,dim = k-1)
      },
      sobol={
        W <- (1/k)*sobol(normal=T,n = N,dim = k-1)
      },
      print('qmc function not specified')
    )
  }
  
  W<-t(apply(W,1,cumsum))
  W1 <- c(rep(0,N),W)
  W1<-matrix(W1,nrow = N,ncol = k)
  S <- rep(K*exp((mu-sigma**2/2)*t),N) 
  S <-matrix(S,nrow=N,ncol=k)
  S1 <-S*exp(sigma*W1)
  if (antithet){
    W2 <- -W1
    S2 <-S*exp(sigma*W2)
    C <- 0.5*exp(-r*T)*(pmax(rowMeans(S1)-K,0)+pmax(rowMeans(S2)-K,0))
  }
  else{
    C <- exp(-r*T)*pmax(rowMeans(S1)-K,0)
  }
  beta <- cov(W1[,k],C)/var(W1[,k])
  C <- C-beta*W1[,k]
  #return (list('mean'=mean(C),'sd'=sd(C)))
  return (c(mean(C),sd(C)))
}
qmc(funcStr = 'sobol')


############# RECAPITULONS ###############
results<-data.frame(matrix(nrow=7,ncol=2),row.names = c('mc','antithetic','ctrl','qmc','qmc antithetic','rqmc','rqmc anti'))
colnames(results)<-c('moyenne','ecart-type')
results[1,]<-monteCarlo(antithet=F,N=10000)
results[2,]<-monteCarlo(antithet = T)
results[3,]<-control()
results[4,]<-qmc(antithet = F,N=10000)
results[5,]<-qmc(antithet = T)
results[6,]<-qmc(antithet = F,randomized = T)
results[7,]<-qmc(antithet = T,randomized = T)
results

#on a donc pas la même moyenne selon les méthodes, ce qui est gênant, 
# par contre les méthodes antithétiques donnent de meilleures variances
# on pourrait rajouter l'option antithétique ou non à la méthode contrôle 

########################################
############# QUESTION 2 ###############
########################################

multiLevel<-function(N,epsilon,T=1,r=0.05,sigma=0.3,K=5,mu){
  if (missing(mu)){mu<-r}
  h0<-T
  L<- as.integer(-log(epsilon/T)/log(2))
  M0 <- as.integer(-log(epsilon)/epsilon**2)
  g <- rnorm(n = N*M0)
  g <- matrix(g,nrow=N,ncol=M0)
  Sg<-K*exp((mu-sigma**2/2)*h0+sigma*sqrt(h0)*g) #approximation lineaire
  estim <- rowMeans(pmax(Sg-K,0)) #à voir si c'est vraiment un rowmean
  for (l in 1:L){
    Ml <- as.integer(M0/2**l) #nombre de simulations au niveau l 
    h<-h0/2**l #pas de discrétisation du niveau l 
    sig <- sigma*sqrt(h) #ecart type correspondant
    Sf <- matrix(1,nrow=N,ncol=Ml) #schema d'euler de pas fin
    Sg <- matrix(1,nrow=N,ncol=Ml) #schema d'euler de pas grossier
    
    for (k in 1:2**(l-1)){
      g1 <- rnorm(n=N*Ml)
      g2 <- rnorm(n=N*Ml)
      #evolution de deux schemas: pas fin et pas grossier
      Sf <- Sf*exp((mu-sigma**2/2)*h+sig*g1)
      Sf <- Sf*exp((mu-sigma**2/2)*h+sig*g2)
      Sg <- Sg*exp(2*(mu-sigma**2/2)*h+sig*(g1+g2))
    }
    Cf<-pmax(Sf-K,0)
    Cg<-pmax(Sg-K,0)
    estim<-estim+rowMeans(Cf-Cg)
    print(length(estim))
  }
  return (exp(-r*T)*estim)
}

e<-multiLevel(N=1000,epsilon = 0.1)
mean(e)
sd(e)
