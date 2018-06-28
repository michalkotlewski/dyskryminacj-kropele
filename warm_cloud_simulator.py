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
    return v_mean * np.log(1 / v_mean - data)

###Inverse transform sampling
#karmi się liczbą superkropelek i srednim promieniem
#wypluwa listę


def inverse_transform_sampling(n, r_sr):
    
    v_sr = 4 * np.pi * r_sr**3 /3
    list = np.random.uniform(0, 1, n)
    
    radius=[]
    
    for i in range(n):
        #function for fsolve, f=0
        def f(x):
            return 1 - np.exp(-x/v_sr) - list[i]
        #finding volume corresponding to generated number
        y = fsolve(f, v_sr)
    
        #radius from volume
        r = (3 * y[0] / (4 * np.pi))**(1/3)
    
        #appending
        radius.append(r)
    return radius


V = 10**6    #wszystko jest w cm3
n0 = 100
N = 1000 
dt = 0.01
ro_w = 999.7  #dla temp 10 C
ro_a = 1.262 #dla wilgotnosci 60%
g = 9.80665 #przyspieszenie graw
nu_a = 17.08*10**(-6) # lepkosc kinematyczna powietrza [Pa·s]
d_0 = 9.06 #stała
c_1 = 0.0902 #stala
mont = 1000#liczba iteracji w Monte Carlo
r_mean = 0.000030531
v_mean = 4 * np.pi/3 * r_mean**3
r_sr = 30.531*10**(-6)

#<<<<<<< HEAD
r_mean=0.000030531
v_mean=4*np.pi/3*r_mean**3
#<<<<<<< HEAD
eta=np.repeat(V*N/n0,n0)
#=======
#=======
#>>>>>>> fb5a89bca6ce4d93303edc06d97d7d9f4d92217c

eta = np.repeat(V * N/n0, n0) #funkcja krotnosci
print(np.shape(eta))
#<<<<<<< HEAD
#>>>>>>> e166a11f9aa094fe72d47f5fbebd9c65968c4d20
#v=inverse_transform_sampling(f(np.arange(0,1,0.001),v_mean),100,N)
#print(v.shape)
#v=f(np.random.uniform(1000),v_mean)
#pl.hist(v)
#r=(3/(4*np.pi)*v)**(1/3)
#=======

r = inverse_transform_sampling(N, r_sr) #wektor promieni
#print(radii)
#>>>>>>> fb5a89bca6ce4d93303edc06d97d7d9f4d92217c
pl.hist(r)
pl.show()


#prędkosc licze w SI, czyli w m/s  
def u(i): #terminal velocities not from index, but from the value
    X = 32*((r[i])**3)*(ro_w-ro_a)*g/(ro_a*nu_a**2)
    
    b_RE = 1/2*c_1*X**0.5*((1+c_1*X**0.5)**0.5-1)**(-1)*(1+c_1*X**0.5)**(-0.5)
    a_RE = d_0**2/4*((1+c_1*X**0.5)**0.5-1)**2/(X**b_RE)
    
    A_nu = a_RE*nu_a**(1-2*b_RE)*((4*ro_w*g)/(3*ro_a))**b_RE
    B_nu = 3*b_RE-1
    
    return A_nu*(2*r[i])**B_nu  


def u2(x): #terminal velocities not from index, but from the value
    X = 32*((x)**3)*(ro_w-ro_a)*g/(ro_a*nu_a**2)
    
    b_RE = 1/2*c_1*X**0.5*((1+c_1*X**0.5)**0.5-1)**(-1)*(1+c_1*X**0.5)**(-0.5)
    a_RE = d_0**2/4*((1+c_1*X**0.5)**0.5-1)**2/(X**b_RE)
    
    A_nu = a_RE*nu_a**(1-2*b_RE)*((4*ro_w*g)/(3*ro_a))**b_RE
    B_nu = 3*b_RE-1
    
    return A_nu*(2*x)**B_nu  


#values of terminal velocities for plot 
rad=[]
for i in range(1,301):
    rad.append(i/1000000)
vel=[]
for i in range(300):
    vel.append(u2(rad[i]))
    
pl.plot(rad,vel,'ro')



                    # TURBULENCJA - DEFINICJE

en_dis = 0.01  # energy dissipation
kolm_length = (nu_a**3/en_dis)**(1/4)
kolm_time=(nu_a/en_dis)**(1/2)
kolm_vel=(nu_a*en_dis)**(1/4)
Rey=200 #Reynolds


          

def relax_time(k): #obliczenia czasu  relaksacji - zalezy od promienia czastki
    return ro_w*ro_a*r[k]**2/(18*nu_a)

def Stokes_num(j):
    return relax_time(j)/kolm_time

def C_1(i): #stala potrzebna do dalszych obliczen
    f1=0.1886*np.exp(20.306/Rey)
    f2=-0.1988*Stokes_num(i)**4 + 1.5275*Stokes_num(i)**3 - 4.2942*Stokes_num(i)**2 + 5.3406*Stokes_num(i)
    return f2/(g*kolm_time/kolm_vel)**f1

def g_12(j, k):
    a_o=(11+7*Rey)/(205+Rey)
    a_og=a_o+np.pi/8*(g**2*kolm_time/kolm_vel)**2
    Fa=20.115*(a_og/Rey)**(1/2)
    r_d=kolm_length*(abs(Stokes_num(j)-Stokes_num(k))*Fa)**(1/2)#length scale of the acceleration diffusion
    return ((kolm_length**2 + r_d**2)/((max(r[j], r[k]))**2 + r_d**2))**(max(C_1(j),C_1(k))/2)

def K(i,j): #koncowy czlon turbulencyjny
    w=(200/(15/en_dis/nu_a)**(1/2))**(1/2)
    return 2*np.pi*(r[i]+r[j])**2*w*g_12(i,j)



#Czlon gravitacyjny

def efficiency(j, k):
    #do rozbudowy, na razie stała
    return 0.5

def prob_real_droplets_g(j, k):  # graw
    return  efficiency(j, k) * np.pi * (r[j]+r[k])**2 *np.absolute(u(j)-u(k)) * dt/V 


def prob_real_droplets_t(j, k):  # turb
    return  efficiency(j, k) * np.pi * (r[j]+r[k])**2 *np.absolute(u(j)-u(k)) * dt/V 
  

def prob_super_droplets_g(j, k):
    return np.maximum(eta[j],eta[k]) * prob_real_droplets_g(j, k)

def prob_super_droplets_linear_g(j, k):
    return N*(N-1)/(2*(N//2)) * prob_super_droplets_g(j, k)

def linear_sampling():
    one = np.arange(n0)
    random = np.random.permutation(one) 
    return np.reshape(random,(-1,2))






ar = linear_sampling()
eps=np.reshape(eta,(50,2))


for i in range(mont):
    for j in range(int(n0/2)):
        psi=np.random.uniform(0,1)
        pdb=prob_super_droplets_linear_g(ar[j,0],ar[j,1])
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
    
    
#Turbulencja - dokładnie to samo, co z gravitacją, tylko zamiana jednego członu
def efficiency(j, k):
    #do rozbudowy, na razie stała
    return 0.5

def prob_real_droplets_g(j, k):  # graw
    return  efficiency(j, k) * K(j,k) * dt/V 


def prob_real_droplets_t(j, k):  # turb
    return  efficiency(j, k) * np.pi * (r[j]+r[k])**2 *np.absolute(u(j)-u(k)) * dt/V 
  

def prob_super_droplets_g(j, k):
    return np.maximum(eta[j],eta[k]) * prob_real_droplets_g(j, k)

def prob_super_droplets_linear_g(j, k):
    return N*(N-1)/(2*(N//2)) * prob_super_droplets_g(j, k)

def linear_sampling():
    one = np.arange(n0)
    random = np.random.permutation(one) 
    return np.reshape(random,(-1,2))

for i in range(mont):
    for j in range(int(n0/2)):
        psi=np.random.uniform(0,1)
        pdb=prob_super_droplets_linear_g(ar[j,0],ar[j,1])
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
                       
        
        
        
    
    
        
        
    
    
    


    
    
    
    







    
    
    
    

