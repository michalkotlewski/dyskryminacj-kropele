# -*- coding: utf-8 -*-


###Main project file

#Import bibliotek
import matplotlib.pyplot as pl



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

r_mean = 30.531*10**(-6)

#Tworzenie probki (n,r)
#n - liczba superkropelek
#r - srednie r
#r_mean = 30.531*10**(-6)
#output - lista
from inverse_transform import inverse_transform_sampling

#Predkosci superkropelek
#wykorzystuja stale!
from velocity import velocities

lista = inverse_transform_sampling(1000000, r_mean)
vel = velocities(lista)


#zderzenia
from coalescence import coal


#rysunki testowe

fig = pl.hist(lista,25)

pl.plot(lista,vel, 'r+')