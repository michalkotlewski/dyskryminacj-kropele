# -*- coding: utf-8 -*-
import numpy as np
import scipy.interpolate as interpolate
import matplotlib.pyplot as pl
from scipy.optimize import fsolve

V = 10**6    #wszystko jest w cm3
n0 = 100*10**6
N = 100000 #docelowo 100000 
dt = 0.1 #moze byc 0.1

ro_w = 999.7  #dla temp 10 C
ro_a = 1.262 #dla wilgotnosci 60%
g = 9.80665 #przyspieszenie graw
nu_a = 17.08*10**(-6) # lepkosc kinematyczna powietrza [Pa·s]
d_0 = 9.06 #stała
c_1 = 0.0902 #stala
mont = 18000#liczba iteracji w Monte Carlo
r_mean = 0.000030531
v_mean = 4 * np.pi/3 * r_mean**3
r_sr = 30.531*10**(-6)

eta = np.repeat(V * n0/N, N)

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

###to chyba nie dziala!
#prędkosc licze w SI, czyli w m/s  
def u(i): #terminal velocities
    X = 32*((r[i]*0.01)**3)*(ro_w-ro_a)*g/(ro_a*nu_a**2)
    
    b_RE = 1/2*c_1*X**0.5*(((1+c_1*X**0.5)**0.5-1)**(-1))*(1+c_1*X**0.5)**(-0.5)
    a_RE = d_0**2/4*((1+c_1*X**0.5)**0.5-1)**2/(X**b_RE)
    
    A_nu = a_RE*nu_a**(1-2*b_RE)*((4*ro_w*g)/(3*ro_a))**b_RE
    B_nu = 3*b_RE-1
    
    return A_nu*(2*r[i]*0.01)**B_nu  

#A to jak najbardziej
def u2(x): #terminal velocities not from index, but from the value
    X = 32*((x)**3)*(ro_w-ro_a)*g/(ro_a*nu_a**2)
    
    b_RE = 1/2*c_1*X**0.5*((1+c_1*X**0.5)**0.5-1)**(-1)*(1+c_1*X**0.5)**(-0.5)
    a_RE = d_0**2/4*((1+c_1*X**0.5)**0.5-1)**2/(X**b_RE)
    
    A_nu = a_RE*nu_a**(1-2*b_RE)*((4*ro_w*g)/(3*ro_a))**b_RE
    B_nu = 3*b_RE-1
    
    return A_nu*(2*x)**B_nu  

#oblicanie jądro turbulencyjnego
    
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






def efficiency(j, k):
    #do rozbudowy, na razie stała
    return 0.5

def prob_real_droplets_g(j, k):  # graw
    K_grav=np.pi*(r[j]+r[k])**2*np.absolute(u2(r[j])-u2(r[k])) 
    return  efficiency(j, k) * (K**2+K_grav**2)**0.5 * dt/V 


def prob_real_droplets_t(j, k):  # turb
    return  efficiency(j, k) * np.pi * (r[j]+r[k])**2 *np.absolute(u2(r[j])-u2(r[k])) * dt/V 
  

def prob_super_droplets_g(j, k):
    return np.maximum(eta[j],eta[k]) * prob_real_droplets_g(j, k)

def prob_super_droplets_linear_g(j, k):
    return N*(N-1)/(2*(N//2)) * prob_super_droplets_g(j, k)


def linear_sampling(n):
    one = np.arange(n)
    random = np.random.permutation(one) 
    return np.reshape(random,(-1,2))

##Test for probabilities
#ar = linear_sampling(N)
#index=ar[1]
#
#r=inverse_transform_sampling(N,r_mean)
#probability=[]
#
#for i in range(len(ar)):
#    index=ar[i]
#    probability.append(prob_super_droplets_linear_g(index[0],index[1]))
#    
#fig1 = pl.hist(r,25)
#fig2 = pl.hist(probability,25)
#
#probability[1]

############################

#GAMMA for multiple coalescence

def gamma(probability, psi):
    pf=np.floor(probability)
    cond=probability-pf
    if psi<cond:
        result=pf+1
        return result
    else:
        return pf



    
##################
##INITIALISATION

r=inverse_transform_sampling(N,r_mean)
r0=[]
r1=[]
r2=[]

for i in range(len(r)):
   r0.append(r[i]) 
#fig1=pl.hist(r,25)
#pl.show()
#fig11=pl.hist(eta,25)
eta0=[]
eta1=[]
eta2=[]

for i in range(len(eta)):
    eta0.append(eta[i])
    
####
##LOOK CAREFULLY!!!    
##MAGIC STARTS HERE:
    
for i in range(mont):
    
    #creating a sample:
    ls=linear_sampling(len(r))
    probab=[]
    for j in range(len(ls)):
        
        #getting index
        if eta[ls[j,0]]>eta[ls[j,1]]:
            k=ls[j,0]
            m=ls[j,1]
        else:
            k=ls[j,1]
            m=ls[j,0]
        
        #drawing a random number:
        psi=np.random.uniform(0,1)
        
        #calculating probability for linear sample:
        pdb=prob_super_droplets_linear_g(ls[j,0],ls[j,1])
        probab.append(pdb)
        #multiple coalescence:
        gm=gamma(pdb,psi)
        if gm !=0:
        #    print(gm)
        
            multi=np.minimum(gm,np.floor(eta[m]/eta[k]))
                
            if eta[k]-multi*eta[m]>0:
                eta[k]=eta[k]-multi*eta[m]
                
                r[m]=(multi*(r[k])**3+(r[m])**3)**(1/3)
                
            else:
                eta[k]=np.floor(eta[m]/2)
                eta[m]=eta[m]-np.floor(eta[m]/2)
                
                r[k]=(multi*(r[k])**3+(r[m])**3)**(1/3)
                r[m]=r[k]
                
    #fig=pl.hist(probab,25)
    #pl.show()
    #czyszczenie
    for l in range(len(r)):
        if eta[l]==0:
            index=eta.index(max(eta))
            eta[l]=np.floor(eta[index]/2)
            r[l]=r[index]
            
            eta[index]=eta[index]-np.floor(eta[index]/2)
    #zapisywanie
    if i==np.floor(mont/3):
        for a in range(len(r)):
            r1.append(r[a])
            eta1.append(eta[a])
    if i==np.floor(2*mont/3):
        for a in range(len(r)):
            r2.append(r[a])
            eta2.append(eta[a])     
            
    ##loop counter (every 100 iteration)
    if np.floor(i%100)==0:
        print(i)
#####################################

#fig1=pl.hist(r,25)
#pl.show()               
#r=[]
#x=[]
#v=[]
#for i in range(300):
#    r.append((i+1)/1000000)
#    v.append(u2(r[i]))
#os x w mikrometrach
#for i in range(len(r)):
#    x.append(r[i]*1000000)
#pl.plot(x,v, 'r+')

#for i in range(10):
#    print(i)
             
pl.hist([r0,r1,r2,r],25, label=["t=0s","t=600s","t=1200s","t=1800s"])
pl.legend()
pl.show()  

####
#Zapisywanie do plikow
with open("r0.txt", 'w') as f:
    for s in range(len(r0)):
        f.write(str(r0[s]) + '\n')

with open("r1.txt", 'w') as f:
    for s in range(len(r1)):
        f.write(str(r1[s]) + '\n')

with open("r2.txt", 'w') as f:
    for s in range(len(r2)):
        f.write(str(r2[s]) + '\n')
        
with open("r3.txt", 'w') as f:
    for s in range(len(r)):
        f.write(str(r[s]) + '\n')

with open("e0.txt", 'w') as f:
    for s in range(len(eta0)):
        f.write(str(eta0[s]) + '\n')

with open("e1.txt", 'w') as f:
    for s in range(len(eta1)):
        f.write(str(eta1[s]) + '\n')

with open("e2.txt", 'w') as f:
    for s in range(len(eta2)):
        f.write(str(eta2[s]) + '\n')
        
with open("e3.txt", 'w') as f:
    for s in range(len(eta)):
        f.write(str(eta[s]) + '\n')
s
