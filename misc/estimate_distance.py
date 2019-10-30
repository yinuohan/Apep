'''
Applies Kepler's laws to estimate the distance of Apep. 
'''

import numpy as np
import matplotlib.pyplot as plt

kpc = 2.063e+8 # AU

M = np.array([10, 15, 20, 30, 50]) # Solar mass
d = np.arange(0.75,2.5,0.01) # kpc
t = 80 # years

r_obs = 47 # mas
e = 0.45
theta = 25 # degrees
a_min = r_obs/(1+e)
a_max = r_obs/(1-e)#/np.cos(theta*np.pi/180)
a = 35


plt.figure(1,figsize=(8,6))
#plt.title('Separation angle for ' + str(t) + '-year period', fontname='serif', fontsize=15)
plt.xlabel('Distance (kpc)', fontname='serif', fontsize=15)
plt.ylabel('Semi-major axis (mas)', fontname='serif', fontsize=15)

for m in M:
    r = m**(1/3) * t**(2/3) # AU
    theta1 = r/(d*kpc) * 180/np.pi * 3600 * 1000 # mas
    theta2 = r/d # some unit
    plt.plot(d, theta1)

plt.plot([d[0],d[-1]],[r_obs, r_obs],'--')
plt.plot([d[0],d[-1]],[a, a],'--')
#plt.plot([d[0],d[-1]],[a_min, a_min],':')
#plt.plot([d[0],d[-1]],[a_max, a_max],':')
plt.fill_between([d[0],d[-1]],[a_min, a_min],[a_max, a_max], alpha=0.2)

import matplotlib.font_manager as font_manager
font = font_manager.FontProperties(family='serif')
plt.legend([str(M[i]) + ' $M_{\odot}$' for i in range(len(M))] + ['NIR separation', 'Inferred semi-major axis'], frameon=True, prop=font)
#plt.legend([str(M[i]) + ' solar mass' for i in range(len(M))] + ['Current visual separation'], frameon=False)
#plt.xscale('log')
#plt.yscale('log')
plt.show()

plt.savefig('C://users//Yinuo//Desktop//High_Distance.png', dpi = 800)