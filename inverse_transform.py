# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import fsolve

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

#for testing
#lista=inverse_transform_sampling(10000,30.531*10**(-6))
#fig = pl.hist(lista,25)