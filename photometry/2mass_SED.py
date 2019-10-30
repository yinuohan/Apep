'''
Author: Yinuo Han

Constructs the SED of 2MASS J16004953-5142506 and allows
for testing with dereddening and blackbody spectra.
'''


import numpy as np
import scipy.interpolate as intp
import matplotlib.pyplot as plt
import os
from mag2jy import mag2jy


## Wavelengths and fluxes
# Surveys on 2mass star
wav_surveys = np.array([0.532, 0.673, 0.797, 1.235, 1.662, 2.159]) # Gaia mean wavelengths (g_blue, g_mean, g_red) from C. Jordi et al 2010
mag_surveys = np.array([13.654579, 12.457168, 11.378937, 9.586,8.681,8.443])

# Ratios from VISIR photometry
wav_VISIR = np.array([8.72, 11.52])
ratio_VISIR = np.array([0.045, 0.012]) * 1e-2
#wav_VISIR = np.array([8.72, 11.52, 19.5])
#ratio_VISIR = np.array([0.045, 0.010023842, 0.006322624]) * 1e-2

# Surveys on Apep (picked relevant lambdas in SIMBAD, discarded outlying points)
wav_apep = np.array([5.73, 8.61, 11.6, 18.4, 22.1])
flux_apep = np.array([13.3, 37.2, np.mean([71.2, 56.7]), 94.2, np.mean([82.8, 84.1])])

# Fit spline through surveys on Apep and sample
tck0 = intp.splrep(wav_apep[0:10],flux_apep[0:10],s=1e1)
flux_apep_sample = intp.splev(wav_VISIR,tck0)

'''
plotx = np.linspace(np.min(wav_apep),np.max(wav_apep),100)
ploty = intp.splev(plotx,tck0)
plt.subplots()
plt.title('Spline through Apep survey fluxes')
plt.scatter(wav_apep, flux_apep)
plt.plot(plotx, ploty)
plt.scatter(wav_VISIR, flux_apep_sample)
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-8,1e4])
plt.show()
'''

# Putting 2mass survey and (Apep survey + flux ratio) together
wavelengths = np.concatenate((wav_surveys, wav_VISIR))
fluxes = np.concatenate((mag2jy(wav_surveys, mag_surveys), flux_apep_sample * ratio_VISIR)) # Jy

'''
plt.subplots()
plt.title('Fluxes of 2mass star - now in Jy')
plt.scatter(wavelengths, fluxes)
plt.xscale('log')
plt.yscale('log')
plt.show()
'''

## Settings - to be fitted later!
#plt.close('all')

wav = 0.55
A_wav = 6.4

T = 30e3
area = 3e-20 # steradians

#T = 30e3
#area = 3e-20 # steradians

out_cloud = 0


## Defaults for dereddening
lam=[2, 4, 23, 41, 73, 91, 120, 130, 150, 180, 200, 218, 240, 260,  280, 330, 365, 440, 550, 700, 900, 1250, 1650, 2200, 3400, 5000, 7000,  9000, 9700, 10000, 12000, 15000, 18000, 20000]
lam=np.array(lam)/1000.0 # now in microns
red=[1.35, 3.39, 7.31, 9.15, 19.1, 17.2, 12.71, 11.09, 9.44, 8.93, 10.08, 11.29, 9.03, 7.63, 6.90, 5.87, 5.53, 4.70, 3.55, 2.66, 1.70, 1.00, 0.624, 0.382, 0.182, 0.095, 0.070, 0.157, 0.208, 0.192, 0.098, 0.053, 0.083, 0.075]

if out_cloud:
    lam=[120, 130, 150, 180, 200, 218, 240, 260, 280, 330, 365, 440, 550, 700, 900, 1250, 1650, 2200, 3400, 5000, 7000, 9000, 9700, 10000,  12000, 15000, 18000, 20000]
    lam=np.array(lam)/1000.0 # now in microns
    red=[5.32, 4.89, 4.57, 4.66, 5.32, 6.03, 5.13, 4.59, 4.34, 4.12, 4.07, 3.67, 3.06, 2.43, 1.70, 1.00, 0.624, 0.382, 0.182, 0.095, 0.070, 0.157, 0.208, 0.192, 0.098, 0.053, 0.083, 0.075]

if(max(wavelengths) > max(lam) or min(wavelengths) < min(lam)):
   print('### ERROR - lam out of interpolation range ###')


## Deredden flux
tck = intp.splrep(lam,red)

red_wav0 = intp.splev(wav,tck)
reddening = intp.splev(wavelengths,tck) * A_wav / red_wav0

star_flx = fluxes/10**(reddening/(-2.5))


## Plot
c = plt.rcParams['axes.prop_cycle'].by_key()['color']

plt.subplots(figsize=(6.5,5))
plt.plot(wavelengths, fluxes, '+:', c=c[0])
plt.plot(wavelengths, star_flx, 'x:', c=c[1])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Wavelength ($\mu$m)', fontname='serif', fontsize=12)
plt.ylabel('Flux (Jy)', fontname='serif', fontsize=12)
#plt.title('SED of 2MASS star at $A_v$ = ' + str(A_wav), fontname='serif', fontsize=12)

plt.tick_params('both', which='both', direction='in', top=True, right=True)
#plt.show()


## Blackbody
from astropy import units as u
from astropy.modeling.blackbody import blackbody_nu

wav2 = np.logspace(np.log10(0.5), np.log10(12)) * u.micrometer
flux = blackbody_nu(wav2,T*u.K).to(u.jansky/u.sr) * area

#plt.plot(wav, flux1)
plt.plot(wav2.to(u.micrometer), flux, '--', c=c[2])

import matplotlib.font_manager as font_manager
font = font_manager.FontProperties(family='serif')

plt.legend(['Visual flux', 'Dereddened flux', 'Blackbody at T = ' + str(int(T)) + ' K'], frameon=False, prop=font)
plt.show()
#plt.text(0.47,0.007,'Area = ' + str(area))

