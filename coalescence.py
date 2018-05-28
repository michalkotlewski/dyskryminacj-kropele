# -*- coding: utf-8 -*-
import numpy as np

#stale - SI!
V=10**6    
n0=100*10**6 #100 cm^-3
N=1000  #liczba superkropelek
dt=0.01
ro_w=999.7  #dla temp 10 C
ro_a=1.262 #dla wilgotnosci 60%
g=9.80665 #przyspieszenie graw
nu_a=17.08*10**(-6) # lepkosc kinematyczna powietrza [Pa·s]
d_0=9.06 #stała
c_1=0.0902 #stala
mont=1000#liczba iteracji w Monte Carlo

def efficiency(j, k):
    #do rozbudowy, na razie stała
    return 0.5

def prob_real_droplets(r, u, j, k):
    return  efficiency(j, k) * np.pi * (r[j]+r[k])**2 *np.absolute(u(j)-u(k)) * dt/V   

def prob_super_droplets(r, u, eta, j, k):
    return np.maximum(eta[j],eta[k]) * prob_real_droplets(r, u, j, k)

def prob_super_droplets_linear(r, u, eta, j, k):
    return N*(N-1)/(2*(N//2)) * prob_super_droplets(r, u, eta, j, k)

def linear_sampling():
    one = np.arange(n0)
    random = np.random.permutation(one) 
    return np.reshape(random,(-1,2))

def probability(r, u, eta, j,k):
    #real droplets
    p_real=efficiency(j, k) * np.pi * (r[j]+r[k])**2 *np.absolute(u(j)-u(k)) * dt/V
    p_super=np.maximum(eta[j],eta[k]) * p_real
    p_ls=N*(N-1)/(2*(N//2)) * p_super
    return p_ls
    
def coal(r, u, eta):
    
    ar = linear_sampling()
    eps=np.reshape(eta,(50,2))
    
    for i in range(mont):
        for j in range(int(n0/2)):
            psi=np.random.uniform(0,1)
            pdb=prob_super_droplets_linear(r, u, eta, ar[j,0],ar[j,1])
            if pdb<psi:
                # if eps[j,0]==1 and if eps[j,1]==1:
                    
                if eps[j,0]==eps[j,1]:
                    eps[j,1]=eps[j,0]/2
                    eps[j,0]=eps[j,0]/2
                    
                    r[ar[j,0]]= (r[ar[j,0]]**3+r[ar[j,0]]**3)**(1/3)
                    r[ar[j,1]]= r[ar[j,0]]
                    
                if eps[j,0]>eps[j,1]:
                    eps[j,0]=eps[j,1]-eps[j,0]
                    r[ar[j,1]]= (r[ar[j,0]]**3+r[ar[j,0]]**3)**(1/3)
    return r, eps