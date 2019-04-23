#coucou
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
W<- apply(W,2,cumsum)
dim(W)

plot(W[,1],type='l',xlim = c(0,20),ylim = c(-3,3))
for (i in 2:n){
  lines(W[,i])
}

plot(colMeans(W))

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
  beta <- cov(W1[,k],C)/var(W1[,k])
  C <- C-beta*W1[,k]
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

qmc<-function(N=1000,T=1,r=0.05,k=20,K=5,sigma=0.3,mu,funcStr='halton',antithet=T){
  if (missing(mu)){
    mu<-r
  } #on suppose que, si non spécifié, le drift et le taux d'intérêt sont égaux
  t <- (1:k)/k
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

############# RANDOMIZED QMC ###############
#aucune idée 

############# RECAPITULONS ###############
results<-data.frame(matrix(nrow=5,ncol=2),row.names = c('mc','antithetic','ctrl','qmc','qmc antithetic'))
colnames(results)<-c('moyenne','ecart-type')
results[1,]<-monteCarlo(antithet=F)
results[2,]<-monteCarlo(antithet = T)
results[3,]<-control()
results[4,]<-qmc(antithet = F)
results[5,]<-qmc(antithet = T)

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
  Sg<-1*(1+mu*h0+sigma*sqrt(h0)*g)
  estim <- colMeans(pmax(S,0)) #à voir si c'est vraiment un colmean
  
  for (l in 1:L){
    Ml <- as.integer(M0/2**l)
    h<-h0/2**l
    sig <- sigma*sqrt(h)
    Sf <- matrix(1,nrow=N,ncol=Ml)
    Sg <- matrix(1,nrow=N,ncol=Ml)
    
    for (k in 1:2**(l-1)){
      g1 <- rnorm(n=N*Ml)
      g2 <- rnorm(n=N*Ml)
      
      Sf <- Sf*(1+mu*h+sig*g1)
      Sf <- Sf*(1+mu*h+sig*g2)
      Sg <- Sg*(1+2*mu*h+sig*(g1+g2))
    }
    Cf<-pmax(Sf-K,0)
    Cg<-pmax(Sg-K,0)
    estim<-estim+mean(Cf-Cg)
  }
  return (estim)
}

multiLevel(N=1000,epsilon = 0.01)

