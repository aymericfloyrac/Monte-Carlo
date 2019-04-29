# -*- coding: utf-8 -*-
from numpy import sqrt, log
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

T = 1.
sig = 0.2
r = 0.05
S0 = 100.
K = 110

def estimMultiniveaux(eps,P):
    """
    Fonction qui renvoie un echantillon de P tirages iid de
    l'estimateur multi-niveaux ayant une precision cible
    d'ordre eps
    """    
    h_0 = T #pas initial
    L = int(-log(eps/T)/log(2.)) #nombre de niveaux
    M_0 = int(-log(eps)/(eps*eps)) #nombre de tirages indep au niveau 0
    
    g = np.random.randn(M_0,P)
    
    Se_gr = S0*(1. + r*h_0 + sig*sqrt(h_0)*g)
    
    estim_multi = np.mean((K-Se_gr)*(K>Se_gr), axis=0)
    print(estim_multi.shape)
    for l in range(1,L+1):
        M_l = int(M_0/2.**l)
    
        ##########################################
        # Parametres du schema d'Euler de pas fin
        ##########################################
        h_fin = h_0/2.**l
        sig_fin = sig*sqrt(h_fin)
        
        Se_fin = S0*np.ones((M_l,P)) #schema d'Euler de pas fin
        Se_gr = S0*np.ones((M_l,P)) #schema d'Euler de pas grossier
        
        for k in range(2**(l-1)):
            g1 = np.random.randn(M_l, P)
            g2 = np.random.randn(M_l, P)
            ###########################################################
            # Evolution des deux schemas avec pas fin et pas grossier
            ###########################################################
            Se_fin = Se_fin * (1. + r*h_fin + sig_fin*g1)
            Se_fin = Se_fin * (1. + r*h_fin + sig_fin*g2)
            
            Se_gr = Se_gr * (1. + r*2*h_fin + sig_fin*(g1+g2))
            
        payoff_fin = (K-Se_fin)*(K > Se_fin)
        payoff_gr = (K-Se_gr)*(K > Se_gr)
        
        estim_multi = estim_multi + np.sum(payoff_fin - payoff_gr, axis=0)/M_l
        
    return estim_multi

eps = 0.2 #precision cible
P = 1000 #nombre de tirages de l'estimateur

###################################################
# Formule explicite pour l'esperance E[(K-S_T)^+]
###################################################
d = (np.log(S0/K) + r*T) / (sig*np.sqrt(T)) + sig*np.sqrt(T)/2.
d2 = d - sig*np.sqrt(T)
esperance = K*norm.cdf(-d2) - S0*np.exp(r*T)*norm.cdf(-d)

plt.clf()

for i in range(3):
    estim_multi = estimMultiniveaux(eps,P)
    
    plt.hist(estim_multi, normed="True", bins=int(sqrt(P)), label='eps = %1.2f' %(eps))
    
    #########################################################
    ## Estimation empirique de l'erreur quadratique
    erreur_quadratique = np.sum( (estim_multi - esperance)**2) / P
    #########################################################
    print("erreur quadratique (eps=%1.2f) : %1.2f" %(eps, np.sqrt(erreur_quadratique)))
    
    eps = eps/2

plt.axvline(esperance, linewidth=2.0, color='k',label="esperance")
plt.legend(loc="best")