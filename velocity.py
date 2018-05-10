# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pl
V=10**6    
n0=100*10**6
N=1000
dt=0.01
ro_w=999.7  #dla temp 10 C
ro_a=1.262 #dla wilgotnosci 60%
g=9.80665 #przyspieszenie graw
nu_a=17.08*10**(-6) # lepkosc kinematyczna powietrza [Pa·s]
d_0=9.06 #stała
c_1=0.0902 #stala
mont=1000#liczba iteracji w Monte Carlo

def u(r,i): #terminal velocities
    X=8*(V*10**(-6))*(ro_w-ro_a)*g/(np.pi*((2*r[i]*0.01)**2)*ro_a*nu_a**2)
    
    b_RE=1/2*c_1*X**0.5*((1+c_1*X**0.5)**0.5-1)**(-1)*(1+c_1*X**0.5)**(-0.5)
    a_RE=d_0**2/4*((1+c_1*X**0.5)**0.5-1)**2/(X**b_RE)
    
    A_nu=a_RE*nu_a**(1-2*b_RE)*((4*ro_w*g)/(3*ro_a))**b_RE
    B_nu=3*b_RE-1
    
    return A_nu*(2*r[i]*0.01)**B_nu  


def velocities(list):
    velo = []
    for i in range(len(list)):
        velo.append(u(list,i))
        
    return velo

##dalej - wykres
r=[]
x=[]
v=[]
for i in range(3000):
    r.append((i+1)/10000000)
    v.append(u(r,i))
#os x w mikrometrach
for i in range(len(r)):
    x.append(r[i]*1000000)
pl.plot(x,v, 'r+')

for i in range(10):
    print(i)