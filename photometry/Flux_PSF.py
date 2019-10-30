'''
Author: Yinuo Han

Calibrates relative fluxes from NACO photometry into
units of Jy based on PSF fluxes in units of Jy. 
Requires calibrator_flux.csv.
'''


import os
import math
import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import scipy.ndimage as img
import scipy.interpolate as intp
from datetime import date
import glob

## Settings
invert = 1

colourmap = 'viridis'

plot_log = 1

#plt.close('all')

## Configurations
wavelengths = [
'IB2.24',
'NB3.74',
'NB4.05',
'H7',
'H9',
'J7',
'J9',
'K7',
'K9',
'L9'
]


#wavelengths = [
#'IB2.24',
#]


## Get files
path = 'Where is the data?'
all_files = glob.glob(path+'stacked_cube'+'*'+'.fits')

image_files = {}
PSF_files = {}

for wav in wavelengths:
    image_files[wav] = [file for file in all_files if (wav in file and ('_PSF' not in file))]
    PSF_files[wav] = [file for file in all_files if (wav in file and ('_PSF' in file))]


## Analysis
image_fluxes = {}
PSF_fluxes = {}
ratios = {}

for wav in wavelengths:
    #print('\n' + wav)
    image_fluxes[wav] = []
    PSF_fluxes[wav] = []
    
    for file in image_files[wav]+PSF_files[wav]:
        # Determine if target or PSF
        if file in image_files[wav]:
            PSF_mode = 0
            #print('\nTarget ' + file[-6])
        elif file in PSF_files[wav]:
            PSF_mode = 1
            #print('\nPSF ' + file[-6])
        
        
        # Read data
        image = fits.getdata(file)
        
        # Get size of image
        c0, c1 = image.shape[0], image.shape[1]
        
        # Define coordinates
        xx = np.tile(np.arange(-c1/2,c1/2),[c0,1])
        yy = np.tile(np.arange(-c1/2,c1/2).reshape(c0,1),[1,c1])
        rr = (xx**2 + yy**2)**(1/2)
        
        # Make windows
        '''
        [Name] = grid of region with 1s and 0s
        [Name]_r = radius of region
        [Name]_c = circle surrounding region for plotting
        '''
        if wav == 'L9':
            everything_r = 50
            bg_r = 60
        else:
            everything_r = 100
            bg_r = 120
        
        everything = rr < everything_r
        everything_c = plt.Circle((c0/2, c1/2), everything_r, color='g', fill=False)
        
        
        bg = (rr < bg_r) * (rr >= everything_r)
        bg_c = plt.Circle((c0/2, c1/2), bg_r, color='g', fill=False)
        
        # Plot and confirm
        '''
        view_image = input("Enter 1 to view image: ")
        if view_image == '1':
            plt.subplots()
            plt.imshow(image*invert, origin='lower', cmap=colourmap)
            plt.gcf().gca().add_artist(everything_c)
            plt.gcf().gca().add_artist(bg_c)
            plt.title('Image of ' + file)
            plt.show()
            
            #time.sleep(5)
            input("Press Enter to continue...")
            plt.close('all')
        '''
        
        # Subtract background from image
        background = np.sum(image*bg)/np.sum(bg)
        image2 = np.array(image - background)
        #print("Background of image: ", background)
        #print("Background of image2: ", np.sum(image2*bg)/np.sum(bg))
        
        # Calculate fluxes
        total_flux = np.sum(image2 * everything)
        #print('Flux: ' + str(total_flux))
        
        if PSF_mode == 0:
            image_fluxes[wav] += [total_flux]
        elif PSF_mode == 1:
            PSF_fluxes[wav] += [total_flux]
    
    ratios[wav] = np.array(image_fluxes[wav]) / np.array(PSF_fluxes[wav])
    print('\nTarget/PSF for', str(wav), str(ratios[wav]))


## Test
'''
wav = 'H9'
num = '_1'
PSF = 'PSF_'
#PSF = ''

image = fits.getdata(path+'stacked_cube_'+PSF+wav+num+'.fits')

c0, c1 = image.shape[0], image.shape[1]
xx = np.tile(np.arange(-c1/2,c1/2),[c0,1])
yy = np.tile(np.arange(-c1/2,c1/2).reshape(c0,1),[1,c1])
rr = (xx**2 + yy**2)**(1/2)

if wav == 'L9':
    everything_r = 40
    bg_r = 55
else:
    everything_r = 100
    bg_r = 120

everything = rr < everything_r
everything_c = plt.Circle((c0/2, c1/2), everything_r, color='g', fill=False)
bg = (rr < bg_r) * (rr >= everything_r)
bg_c = plt.Circle((c0/2, c1/2), bg_r, color='g', fill=False)
        
plt.subplots()
plt.imshow(image*invert, origin='lower', cmap=colourmap)
plt.gcf().gca().add_artist(everything_c)
plt.gcf().gca().add_artist(bg_c)
plt.title('Image of ' + wav + num)
plt.show()

background = np.sum(image*bg)/np.sum(bg)
image2 = np.array(image - background)
print("Background of image: ", background)
print("Background of image2: ", np.sum(image2*bg)/np.sum(bg))
        
total_flux = np.sum(image2 * everything)
print('Flux: ' + str(total_flux))
'''

## Physical flux
cal_names = {
'IB2.24': 'HD 144648',
'NB3.74': 'HD 144648',
'NB4.05': 'HD 144648',
'H7': '[W71b] 113-03',
'H9': 'HD 142489',
'J7': '[W71b] 113-03',
'J9': 'HD 142489',
'K7': '[W71b] 113-03',
'K9': 'HD 142489',
'L9': 'IRAS 15539-5219'
}

microns = {
'IB2.24': 2.24,
'NB3.74': 3.740,
'NB4.05': 4.051,
'H7': 1.66,
'H9': 1.66,
'J7': 1.265,
'J9': 1.265,
'K7': 2.18,
'K9': 2.18,
'L9': 3.80
}

inv_microns = {
2.24: ['IB2.24'],
3.740: ['NB3.74'],
4.051: ['NB4.05'],
#1.66: ['H7', 'H9'],
1.66: ['H7'],
#1.265: ['J7', 'J9'],
1.265: ['J7'],
#2.18: ['K7', 'K9'],
2.18: ['K7'],
3.80: ['L9']
}

target_microns = []
target_pfluxes = []

df = pd.read_csv(path+'calibrator_flux.csv')

for micron in list(inv_microns.keys()):
    target_microns.append(micron)
    data_points = 0
    target_pflux = 0
    for wav in inv_microns[micron]:
        cal_pflux = float(df[cal_names[wav]][df['wav_microns']==micron])
        target_pflux += cal_pflux * np.sum(ratios[wav])
        data_points += len(ratios[wav])
    target_pfluxes.append(target_pflux/data_points)
    

plt.subplots(1)
plt.plot(target_microns, target_pfluxes, 'X')
plt.title('SED of Apep based on NACO')
plt.xlabel('Wavelengths (microns)')
plt.ylabel('Flux (Jy)')
if plot_log == 1:
    plt.xscale('log')
    plt.yscale('log')
plt.show()
