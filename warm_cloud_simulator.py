# -*- coding: utf-8 -*-
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as pl
from scipy.optimize import fsolve
"""
def inverse_transform_sampling(data, n_bins=40, n_samples=10000):
    hist, bin_edges = np.histogram(data, bins=n_bins, density=True)
    cum_values = np.zeros(bin_edges.shape)
    cum_values[1:] = np.cumsum(hist*np.diff(bin_edges))
    inv_cdf = interpolate.interp1d(cum_values, bin_edges)
    r = np.random.uniform(0,1,n_samples)
    #r = np.random.rand(n_samples)
    return inv_cdf(r)
"""



##Slawkowe wtf
#def inverse_transform_sampling(data, n_bins=40, n_samples=1000):
#    hist, bin_edges = np.histogram(data, bins=n_bins, density=True)
#    cum_values = np.zeros(bin_edges.shape)
#    cum_values[1:] = np.cumsum(hist*np.diff(bin_edges))
#    inv_cdf = interpolate.interp1d(cum_values, bin_edges)
#    r = np.random.rand(n_samples)
#    return inv_cdf(r)


#def f(data, v_mean):
#    return 1/v_mean*np.exp(-data/v_mean)



def f(data, v_mean):
    return v_mean*np.log(1/v_mean-data)

###Inverse transform sampling
#karmi się liczbą superkropelek i srednim promieniem
#wypluwa listę
#r_sr = 30.531*10**(-6)

def inverse_transform_sampling(n, r_sr):
    
    v_sr = 4 * np.pi * r_sr**3 /3
    list = np.random.uniform(0,1,n)
    
    radius=[]
    for i in range(n):
        #function for fsolve, f=0
        def f(x):
            return 1-np.exp(-x/v_sr)-list[i]
        #finding volume corresponding to generated number
        y = fsolve(f,v_sr)
    
        #radius from volume
        r = (3*y[0]/(4 * np.pi))**(1/3)
    
        #appending
        radius.append(r)
    return radius


#def u(r,alfa,beta):
#    return alfa*r**beta

V=10**6    #wszystko jest w cm3
n0=100
N=1000 
dt=0.01
ro_w=999.7  #dla temp 10 C
ro_a=1.262 #dla wilgotnosci 60%
g=9.80665 #przyspieszenie graw
nu_a=17.08*10**(-6) # lepkosc kinematyczna powietrza [Pa·s]
d_0=9.06 #stała
c_1=0.0902 #stala
mont=1000#liczba iteracji w Monte Carlo



r_mean=0.000030531
v_mean=4*np.pi/3*r_mean**3
#<<<<<<< HEAD
eta=np.repeat(V*N/n0,n0)
#=======

eta=np.repeat(V*N/n0,n0)
print(np.shape(eta))
#>>>>>>> e166a11f9aa094fe72d47f5fbebd9c65968c4d20
v=inverse_transform_sampling(f(np.arange(0,1,0.001),v_mean),100,N)
print(v.shape)
#v=f(np.random.uniform(1000),v_mean)
#pl.hist(v)
r=(3/(4*np.pi)*v)**(1/3)
pl.hist(r)
pl.show()
#print(f(np.arange(0.004),v_mean))
#print(v) 
#print(r)







#prędkosc licze w SI    
def u(i): #terminal velocities
    X=32*((r[i]*0.01)**3)*(ro_w-ro_a)*g/(ro_a*nu_a**2)
    
    b_RE=1/2*c_1*X**0.5*(((1+c_1*X**0.5)**0.5-1)**(-1))*(1+c_1*X**0.5)**(-0.5)
    a_RE=d_0**2/4*((1+c_1*X**0.5)**0.5-1)**2/(X**b_RE)
    
    A_nu=a_RE*nu_a**(1-2*b_RE)*((4*ro_w*g)/(3*ro_a))**b_RE
    B_nu=3*b_RE-1
    
    return A_nu*(2*r[i]*0.01)**B_nu  

def u2(x): #terminal velocities not from index, but from the value
    X=32*((x)**3)*(ro_w-ro_a)*g/(ro_a*nu_a**2)
    
    b_RE=1/2*c_1*X**0.5*((1+c_1*X**0.5)**0.5-1)**(-1)*(1+c_1*X**0.5)**(-0.5)
    a_RE=d_0**2/4*((1+c_1*X**0.5)**0.5-1)**2/(X**b_RE)
    
    A_nu=a_RE*nu_a**(1-2*b_RE)*((4*ro_w*g)/(3*ro_a))**b_RE
    B_nu=3*b_RE-1
    
    return A_nu*(2*x)**B_nu  


#values of terminal velocities for plot
rad=[]
for i in range(1,301):
    rad.append(i/1000000)
vel=[]
for i in range(300):
    vel.append(u2(rad[i]))
    
pl.plot(rad,vel,'ro')




def efficiency(j, k):
    #do rozbudowy, na razie stała
    return 0.5

def prob_real_droplets(j, k):
    return  efficiency(j, k) * np.pi * (r[j]+r[k])**2 *np.absolute(u(j)-u(k)) * dt/V   

def prob_super_droplets(j, k):
    return np.maximum(eta[j],eta[k]) * prob_real_droplets(j, k)

def prob_super_droplets_linear(j, k):
    return N*(N-1)/(2*(N//2)) * prob_super_droplets(j, k)

def linear_sampling():
    one = np.arange(n0)
    random = np.random.permutation(one) 
    return np.reshape(random,(-1,2))





ar = linear_sampling()
eps=np.reshape(eta,(50,2))


for i in range(mont):
    for j in range(int(n0/2)):
        psi=np.random.uniform(0,1)
        pdb=prob_super_droplets_linear(ar[j,0],ar[j,1])
        if pdb<psi:
            #if eps[j,0]==1 and if eps[j,1]==1:
                
                
            if eps[j,0]==eps[j,1]:
                eps[j,1]=eps[j,0]/2
                eps[j,0]=eps[j,0]/2
                
                r[ar[j,0]]= (r[ar[j,0]]**3+r[ar[j,0]]**3)**(1/3)
                r[ar[j,1]]= r[ar[j,0]]
                
            if eps[j,0]>eps[j,1]:
                eps[j,0]=eps[j,1]-eps[j,0]
                r[ar[j,1]]= (r[ar[j,0]]**3+r[ar[j,0]]**3)**(1/3)
                
                
        
       
        
    #else:
     #   delta=eps[i,0]-eps[i,1]
      #  eps[i,1]=int(pdb*eps[i,1])
       # eps[i,0]=int(pdb*eps[i,0])+delta
        
   # eps_array=np.reshape(eps, (1,N))
    
    #den=1/V*np.sum(eps_array) #droplet number density  
    
    

        
        
        
        
    
    
        
        
    
    
    


    
    
    
    







    
    
    
    

