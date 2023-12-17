# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 13:29:42 2023

@author: nazli
"""
import numpy as np
import random as r
import matplotlib.pyplot as plt

#Parameters
N = 100 # number of the sites on the latice
l = [0] * N
a = 0.2  # initiation rate
b = 0.8 # termination rate
p = 1   # elongation rate
light_intensity = []
rhos_compared=[]
rhos=[]
codons = {i:[]for i in range(N)}

# Simulate TASEP
def TASEP(t):
    forced = False
    S = a * (1 - l[0]) + b * l[N - 1]
    for i in range(0, N - 1):
        S += p * l[i] * (1 - l[i + 1])
    k = r.uniform(0, S)
    Sp = 0
       
    if l[0] == 0:#initiation alpha first site
        Sp += a
        if k <= Sp and not forced:
            l[0] = 1
            forced = True
                
                
    if l[0] == 1 and l[1] == 0:#site 1
        Sp += p
        if k <= Sp and not forced:
            l[0] = 0
            l[1] = 1
            forced = True
        
    if l[N - 2] == 1 and l[N - 1] == 0:#site N-1 
        Sp += p
        if k <= Sp and not forced:
            l[N - 1] = 1
            l[N - 2] = 0
            forced = True
          
    if l[N - 1] == 1:#termination beta site N
        Sp += b
        if k <= Sp and not forced:
            l[N - 1] = 0
            forced = True
           
    if not forced:
        
        for i in range(1, N-2):#elongation the part that is invariant by translation
            if l[i] == 1 and l[i+1] == 0:
                Sp += p
                if k <= Sp and not forced:
                    l[i] = 0
                    l[i+1] =1
                    forced = True 
    if forced:
        d = 0
        for i in range(N):
            if l[i] == 1:
                d += i+1
                codons[i].append(1)
        light_intensity.append(d)
        

    ribosomes = sum(l)
    rho = ribosomes / N  # density
    rhos.append(rho)  

    return(rhos)
t = 10000
for i in range(t):
    TASEP(t)

#codon density
codon = list(codons.keys())
average_occupancy = [sum(values) / t for values in codons.values()]

#Plotting
def plot_results(x, theo, moyenne):
    plt.plot(x, rhos , 'b', label='Density vs. Time')
    plt.plot(x, theo, 'r-', linewidth= 0.8, label='Theoritical value')
    plt.xlabel('Time')
    plt.ylabel('Density')
    plt.title('TASEP Simulation-Low Density')
    plt.annotate(f'a={a}, b={b}, p={p}, Ave.Den={moyenne}', xy=(0.55, 0.3), xycoords='axes fraction', ha='center', va='center')

    plt.legend()
    plt.show()   
def plot_light(x, y2, expected):
    plt.plot(x, y2, 'g', label = 'Light Intensity vs. Time')
    plt.plot(x, expected, 'r', linewidth = 0.8, label='Theoritical value')
    plt.xlabel('Time')
    plt.ylabel('Light Intensity')
    plt.title('TASEP Simulation-Suntag')
    plt.legend()
    plt.show() 
def plot_codon_occupancy(x1, y):
    plt.plot(codon, average_occupancy, marker= 'o', linestyle='-', color='#FF5733')
    plt.xlabel('Position')
    plt.ylabel('Average codon occupancy')
    plt.ylim(0, 0.6)
    plt.show()
def plot_density_distribution(h):
    
    plt.hist(h, bins=None, range=None, density=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=0.8, log=False, color=None, label='Density distribution', stacked=False)
def steady_state(stat):
    average_in_steady_state = sum(rhos[stat:]) / (t - stat)
    print("Average Density in Steady State:", average_in_steady_state)

moyenne = sum(rhos) / t
print("Average Density:", moyenne)
print('Final state',l)

x=np.linspace(0,t,len(rhos))
y2 = light_intensity
expected = [1000]*len(light_intensity)
theo = [a/p] * len(rhos)
x1 = np.linspace(0, N, 1)
y = np.linspace(0,1,len(average_occupancy))
h = rhos
plot_results(x, theo, moyenne)
plot_light(x, y2, expected)
plot_codon_occupancy(x1, y)
plot_density_distribution(h)
steady_state(1000)

#comparing multiple lattices
# s=50
# for i in range(s):
#     rhos=[]
#     for j in range(t):
#         TASEP(t)
#     rhos = TASEP(t)
#     moyenne = sum(rhos) / t
#     rhos_compared.append(moyenne)
#     # codon = list(codons.keys())
#     # average_occupancy = [sum(values) / t for values in codons.values()]

# def plot_density_distribution_compared(h1):
    
#     plt.hist(h1, bins=None, range=None, density=False, weights=None, cumulative=False, bottom=None, histtype='bar', align='mid', orientation='vertical', rwidth=0.8, log=False, color=None, label='Density distribution', stacked=False)
    
# h1=rhos_compared
# plot_density_distribution_compared(h1)
# print(rhos_compared)
# moyennes = sum(rhos_compared) / s
# print("Average Density of multiple lattices:", moyennes)