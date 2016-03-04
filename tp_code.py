#Question 1: Formule de Black & Scholes

from scipy.stats import norm

def BlackScholes(CallPutFlag,S,K,T,r,v):
    d1 = (log(float(S)/K)+((r)+v*v/2.)*T)/(v*sqrt(T))
    d2 = d1-v*sqrt(T)
    if CallPutFlag=='c':
        return S*exp(-d*T)*norm.cdf(d1)-K*exp(-r*T)*norm.cdf(d2)
    else:
        return K*exp(-r*T)*norm.cdf(-d2)-S*exp(-d*T)*norm.cdf(-d1)

# prix en fonction de T
p=[]
for i in xrange(0.5,10,0.1):
   p.append([i,BlackScholes(CallPutFlag,S,K,i,r,v)])

# prix en fonction de S0
p=[]
for i in xrange(1,200,1):
    p.append([i,BlackScholes(CallPutFlag,i,K,T,r,v)])

# prix en fonction de sigma
p=[]
for i in xrange(0.01,1,0.01):
   p.append([i,BlackScholes(CallPutFlag,S,K,T,r,i)])

#export to csv
with open('fichier.csv', 'w', newline='') as fp:
    a = csv.writer(fp, delimiter=',')
    data = p
    a.writerows(data)
	
#Question 2: Grecques d'un risk-reversal

def Black_Scholes_Greeks_Call(S, K, r, v, T):
    T_sqrt = sqrt(T)
    d1 = (log(float(S)/K)+((r)+v*v/2.)*T)/(v*T_sqrt)
    d2 = d1-v*T_sqrt
    Delta = norm.cdf(d1)
    Gamma = norm.pdf(d1)/(S*v*T_sqrt)
    Theta =- (S*v*norm.pdf(d1))/(2*T_sqrt) - r*K*exp( -r*T)*norm.cdf(d2)
    Vega = S * T_sqrt*norm.pdf(d1)
    Rho = K*T*exp(-r*T)*norm.cdf(d2)
    return Delta, Gamma, Theta, Vega, Rho
    
def Black_Scholes_Greeks_Put(S, K, r, v, T):
    T_sqrt = sqrt(T)
    d1 = (log(float(S)/K)+r*T)/(v*T_sqrt) + 0.5*v*T_sqrt
    d2 = d1-(v*T_sqrt)
    Delta = norm.cdf(d1)
    Gamma = norm.pdf(d1)/(S*v*T_sqrt)
    Theta = -(S*v*norm.pdf(d1)) / (2*T_sqrt)- r*K * exp(-r*T) * norm.cdf(d2)
    Vega = S * T_sqrt * norm.pdf(d1)
    Rho = -K*T*exp(-r*T) * norm.cdf(d2)
    return Delta, Gamma, Theta, Vega, Rho

def Risk_Reversal_Greeks(S,K1,K2,r,v,T):
    Delta = Black_Scholes_Greeks_Put(S, K1, r, v, T)[0] - Black_Scholes_Greeks_Call(S, K2, r, v, T)[0]
    Gamma = Black_Scholes_Greeks_Put(S, K1, r, v, T)[1] - Black_Scholes_Greeks_Call(S, K2, r, v, T)[1]
    Theta = Black_Scholes_Greeks_Put(S, K1, r, v, T)[2] - Black_Scholes_Greeks_Call(S, K2, r, v, T)[2]
    Vega = Black_Scholes_Greeks_Put(S, K1, r, v, T)[3] - Black_Scholes_Greeks_Call(S, K2, r, v, T)[3]
    Rho= Black_Scholes_Greeks_Put(S, K1, r, v, T)[4] - Black_Scholes_Greeks_Call(S, K2, r, v, T)[4]
    return Delta, Gamma, Theta, Vega, Rho
	
p=[]
for i in drange(0.1,200,0.5):
   p.append([i,Risk_Reversal_Greeks(i,K1,K2,r,v,T)[0],Risk_Reversal_Greeks(i,K1,K2,r,v,T)[1],Risk_Reversal_Greeks(i,K1,K2,r,v,T)[2],Risk_Reversal_Greeks(i,K1,K2,r,v,T)[3],Risk_Reversal_Greeks(i,K1,K2,r,v,T)[4]])

#Question 3: Arbre Binomial

from math import exp
import numpy as np

def BinomialTree3(type,S0, K, r, sigma, T, N=2000,american="false"):   
    deltaT = float(T) / N
	u = np.exp(sigma * np.sqrt(deltaT))
    d = 1.0 / u
    #init the arrays using numpy
    fs =  np.asarray([0.0 for i in range(N + 1)])   
    #stock tree for calculations of expiration values
    fs2 = np.asarray([(S0 * u**j * d**(N - j)) for j in range(N + 1)])
    #we vectorize the strikes as well so the expiration check will be faster
    fs3 =np.asarray( [float(K) for i in range(N + 1)])
	
    a = np.exp(r * deltaT)
    p = (a - d)/ (u - d)
    oneMinusP = 1.0 - p
 
   
    # Compute the leaves, f_{N, j}
    if type =="C":
        fs[:] = np.maximum(fs2-fs3, 0.0)
        
    else:
        fs[:] = np.maximum(-fs2+fs3, 0.0)
        
    #calculate backward the option prices
    for i in range(N-1, -1, -1):
       fs[:-1]=np.exp(-r * deltaT) * (p * fs[1:] + oneMinusP * fs[:-1])
       fs2[:]=fs2[:]*u
      
       if american=='true':
           #Check if we should exerce option
           if type =="C":
                fs[:]=np.maximum(fs[:],fs2[:]-fs3[:])
           else:
                fs[:]=np.maximum(fs[:],-fs2[:]+fs3[:])
                if fs[:].any()<(-fs2[:]+fs3[:]).any():
                    print('h')
                
    return fs[0]

	
#Question 4: Calcul du delta

def deltaDifferencesFinis(S,epsilon,N):
    return (binomialCallAmerican(S+epsilon,K,dt,r,sigma,N)-binomialCallAmerican(S,K,dt,r,sigma,N))/(epsilon)

#Question 5: Monte-Carlo

from math import *
from random import gauss

def ST(Zi):
    return S0 * exp((r-(sigma**2)/2)*T+ sigma* sqrt(T)*Zi)
    

def mc(n):
    resultat=0.0
    for i in xrange(1,n,1):
        resultat += max(ST(gauss(0,1))-K,0)
    #Actualisation : multiplier par exp(-rT)
    return resultat/n * exp(-r*T)
	

# Pour le calcul de delta
def ST_so(Zi,S0):
    return S0 * exp((r-(sigma**2)/2)*T+ sigma* sqrt(T)*Zi)
    
def delta():
    return (ST_s0(gauss(0,1),S0+epsilon)-ST_s0(gauss(0,1),S0))/(epsilon)
    
def mcdelta(n):
    resultat=0.0
    for i in xrange(1,n,1):
        resultat += delta()
    return resultat/n
	
#Question 6

#Monte-Carlo
def STBarriere(Zi,barriere):
    ST = S0 * exp((r-(sigma**2)/2)*T+ sigma* sqrt(T)*Zi)
    if ST > barriere:
        return 0
    else:
        return ST
        
def mcbarriere(n,barriere):
    resultat=0.0
    for i in range(1,n,1):
        resultat += max(STBarriere(gauss(0,1),barriere)-K,0)
    #Actualisation : multiplier par exp(-rT)
    return resultat/n * exp(-r*T)
	
#Arbre Binomial
def binomialBarriere(s,x,T,r,sigma,n,barriere):
    deltaT=T/n
    u=exp(sigma * sqrt(deltaT))
    d=1.0 / u
    a = exp(r*deltaT)
    p=(a-d) / (u-d)
    # create array to stock values
    v = [[0.0 for j in range(i + 1)] for i in range (n+1)]
    #calculate the payoff
    for j in range(n+1):
        new_s = s* u**j * d**(n-j)
        if new_s > barriere:
            new_s = 0
        v[n][j] = max(new_s-x,0.0)
    #recursively deduce the previous values
    for i in range(n-1,-1,-1):
        for j in range(i + 1):
            v1=exp(-r*deltaT)*(p*v[i+1][j+1]+(1.0-p)*v[i+1][j])
            v[i][j]=v1
    resultat = v[0][0]
    return resultat   

    

