# -*- coding: utf-8 -*-
import numpy as np
import scipy.interpolate as interpolate

def inverse_transform_sampling(data, n_bins=40, n_samples=1000):
    hist, bin_edges = np.histogram(data, bins=n_bins, density=True)
    cum_values = np.zeros(bin_edges.shape)
    cum_values[1:] = np.cumsum(hist*np.diff(bin_edges))
    inv_cdf = interpolate.interp1d(cum_values, bin_edges)
    r = np.random.rand(n_samples)
    return inv_cdf(r)


def f(data, v_mean):
    return 1/v_mean*np.exp(-data/v_mean)

def u(r,alfa,beta):
    return alfa*r**beta

V=10**12    #wszystko jest w cm3
n0=100
N=100
dt=0.01
r_mean=0.0030531
v_mean=4*np.pi/3*r_mean**3
eta=np.repeat(V*n0/N,N)
v=inverse_transform_sampling(f(np.arange(0,1,0.01),v_mean),100,N)
r=(3/(4*np.pi)*v)**(1/3)
#print(f(np.arange(0.004),v_mean))
print(v)
print(r)

def efficiency(j, k):
    #do rozbudowy, na razie stała
    return 0.5

#do ustalenia wartosci stałych!!!
alfa = 0.5
beta = 0.5

def prob_real_droplets(j, k):
    return  efficiency(j, k) * np.pi * (r[j]+r[k])**2 *np.absolute(u(r[j],alfa, beta)-u(r[k], alfa, beta)) * dt/V   

def prob_super_droplets(j, k):
    return np.maximum(eta[j],eta[k]) * prob_real_droplets(j, k)

def prob_super_droplets_linear(j, k):
    return N*(N-1)/(2*(N//2)) * prob_super_droplets(j, k)

def linear_sampling():
    one = np.arange(N)
    random = np.random.permutation(one) 
    return np.reshape(random,(-1,2))

