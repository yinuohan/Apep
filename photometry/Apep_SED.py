'''
Author: Yinuo Han

Constructs the dereddened SED of Apep and plots 
blackbody spectra. 
Requires Apep_photometry.xlsx.
'''


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.interpolate import splev, splrep
from scipy.interpolate import interp1d
from astropy.io import fits
import os
from mag2jy import mag2jy
from unred import unred

## Settings
plot_log = 1


## Read survey data
path = 'Where is the data?'
df = pd.read_excel(path+'Apep_photometry.xlsx')
df = df.sort_values('wavelength(micron)')

nrows = df.shape[0]
ncols = df.shape[1]

# Get physical flux
wavelengths = df['wavelength(micron)']
raw = df['pflux(Jy)']
weight = df['weight']

# Mask
rawmask = np.isfinite(raw) * df['Total'] == 1

# Deredden
Av = 9.8
raw_unred = unred(raw, wavelengths, Av, 0.55)
spl = splrep(wavelengths[rawmask], raw_unred[rawmask], w=weight[rawmask], s=5)

# Interpolate missing points
intp_unred = raw_unred.copy()
for i in range(intp_unred.shape[0]):
    if np.isnan(intp_unred.iloc[i]) and df['Total'].iloc[i] == 1:
        print (i)
        intp_unred.iloc[i] = splev(wavelengths.iloc[i], spl)

# Get component flux
total = intp_unred * df['Total']
C = intp_unred * df['C']
C_upper = intp_unred * df['C_upper']
N = intp_unred * df['N']
Dust = intp_unred * df['Dust']
Dust_upper = intp_unred * df['Dust_upper']

## Masks
totalmask = np.isfinite(total)
Cmask = np.isfinite(C)
C_uppermask = np.isfinite(C_upper)
Nmask = np.isfinite(N)
Dustmask = np.isfinite(Dust)
Dust_uppermask = np.isfinite(Dust_upper)

## Begin plot
fig, ax = plt.subplots(1,figsize=(7.5,6))

c = plt.rcParams['axes.prop_cycle'].by_key()['color']

#c = ['b', 'y', 'g', 'r']


## Plot blackbody
from astropy import units as u
from astropy.modeling.blackbody import blackbody_nu

# T in Kelvins, Area in steradians
smooth_xshooter = 1

# Generate blackbody for WR
T3 = 70e3
area3 = 1.2e-20

T4 = 11e2
area4 = 2e-16

wav = np.logspace(np.log10(0.5), np.log10(5)) * u.micrometer
flux3 = blackbody_nu(wav,T3*u.K).to(u.jansky/u.sr) * area3 * u.steradian
flux4 = blackbody_nu(wav,T4*u.K).to(u.jansky/u.sr) * area4 * u.steradian

plt.plot(wav.to(u.micrometer), (flux3+flux4).to(u.jansky), '-', c=c[1], label='_nolegend_')

plt.plot(wav.to(u.micrometer), flux3.to(u.jansky), ':', c=c[1], label='_nolegend_')

plt.plot(wav.to(u.micrometer), flux4.to(u.jansky), ':', c=c[1], label='_nolegend_')

# Generate blackbody for O8
T1 = 34e3
area1 = 5e-20

wav = np.logspace(np.log10(0.5), np.log10(20)) * u.micrometer
flux = blackbody_nu(wav,T1*u.K).to(u.jansky/u.sr) * area1 * u.steradian

plt.plot(wav.to(u.micrometer), flux.to(u.jansky), ':', c=c[2], label='_nolegend_')

# Generate blackbody for Dust
T2 = 2.7e2
area2 = 3.5e-13

wav = np.logspace(np.log10(2.8), np.log10(30)) * u.micrometer
flux = blackbody_nu(wav,T2*u.K).to(u.jansky/u.sr) * area2 * u.steradian

plt.plot(wav.to(u.micrometer), flux.to(u.jansky), '-.', c=c[3], label='_nolegend_')


## Plot components

# Total
xt = np.array(wavelengths[totalmask])
yt = np.array(total[totalmask])
wt = weight[totalmask]

st = splrep(np.log(xt), np.log(yt), w=wt, s=2)
xt2 = np.linspace(np.log(xt)[0], np.log(xt)[-1], 1000)
yt2 = splev(xt2, st)
xt2 = np.exp(xt2)
yt2 = np.exp(yt2)

plt.plot(xt, yt, 'X', c=c[0])
plt.plot(xt2, yt2, '--', c=c[0], label='_nolegend_')


# North
xn = np.array(wavelengths[Nmask])
yn = np.array(N[Nmask])

plt.plot(xn, yn, 's', c=c[2])

# Central
xc = np.array(wavelengths[Cmask])
yc = np.array(C[Cmask])

plt.plot(xc, yc, '*', c=c[1], markersize=8)

# Dust
xd = np.array(wavelengths[Dustmask])
yd = np.array(Dust[Dustmask])

xdu = np.array(wavelengths[Dust_uppermask])
ydu = np.array(Dust_upper[Dust_uppermask])

plt.plot(xd, yd, 'o', c=c[3])
plt.plot(xdu, ydu, 'v', c=c[3])




## Xshooter

path = 'Where is the data?'

C_file = 'median_apep_vis_central_eng'
N_file = 'median_apep_vis_north_com'

C = fits.getdata(path + C_file + '.fits')
N = fits.getdata(path + N_file + '.fits')

XC = np.arange(533.66,533.66+0.02*len(C),0.02)*1e-3
XN = np.arange(533.66,533.66+0.02*len(N),0.02)*1e-3

C = unred(C, XC, 11.4, 0.55)
N = unred(N, XN, 11.4, 0.55)

inunit = u.erg/u.s/(u.cm)**2/u.angstrom
cc = 3e8 * u.m/u.s

XC = XC*u.um
XN = XN*u.um

C = (C*inunit * XC**2/cc).to(u.Jy)
N = (N*inunit * XN**2/cc).to(u.Jy)

def sliding_mean(data_array, window=5):
    data_array = np.array(data_array)
    new_list = []
    for i in range(len(data_array)):
        indices = range(max(i - window + 1, 0),  
                        min(i + window + 1, len(data_array)))
        avg = 0
        for j in indices:
            avg += data_array[j]
        avg /= float(len(indices))
        new_list.append(avg)
    return np.array(new_list)

if smooth_xshooter == 1:
    C = sliding_mean(C,window=200)
    N = sliding_mean(N,window=200)

lim = [10000, 13000, 18000]

plt.plot(XC[lim],C[lim],'*--',c=c[1], markersize=8)
plt.plot(XN[lim],N[lim],'s--',c=c[2])


## Show plot
if plot_log == 1:
    plt.xscale('log')
    plt.yscale('log')

#plt.title('SED of Apep dereddened with $A_v$ = '+str(Av), fontname='serif', fontsize=15)
plt.xlabel('Wavelength ($\mu$m)', fontname='serif', fontsize=12)
plt.ylabel('Flux (Jy)', fontname='serif', fontsize=12)

plt.tick_params('both', which='both', direction='in', top=True, right=True)

plt.ylim([1e-3, 1e3])

import matplotlib.font_manager as font_manager
font = font_manager.FontProperties(family='serif')

plt.legend(['Total', 'Northern companion', 'Central binary', 'Dust', 'Dust upper limit'], frameon=False, prop=font)

plt.show()


## Test unred
'''
wav = np.arange(0.5,70,0.01)
flux = np.ones(len(wav))

plt.figure(2)
plt.plot(wav,flux)
plt.plot(wav,unred(flux,wav,11.4,0.55))
plt.xscale('log')
plt.yscale('log')
plt.show()
'''
