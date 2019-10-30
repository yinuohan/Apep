'''
Author: Yinuo Han

Performs aperture photometry on the 2019 epoch of 
NACO images of Apep. 
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

plt.close('all')

path = 'Where is the data?'

plt.rcParams.update({'font.size': 18})


## Configuration
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
#'L9'
]

radii = {
# central, northern, everything, background
'IB2.24': [17, 12, 150, 200],
'NB3.74': [15, 7, 150, 200],
'NB4.05': [15, 7, 150, 200],
'H7': [26, 26, 100, 120],
'H9': [28, 28, 100, 120],
'J7': [24, 28, 100, 120],
'J9': [24, 28, 100, 120],
'K7': [29, 27, 100, 120],
'K9': [29, 27, 100, 120],
'L9': [20, 7, 40, 55],
}

pixel = 13.221/1000
low = -13.221 * 128/1000
high = 13.221 * 128/1000
extent = [low, high, low, high]

## Test
wav = 'H9'
num = '_1'

# Read data
image = fits.getdata(path+'stacked_cube_'+wav+num+'.fits')

# Normalise image
image /= np.max(image)

# Get image dimensions
c0, c1 = image.shape[0], image.shape[1]
b0, b1 = int(c0/2), int(c1/2)

# Make primary windows
xx = np.tile(np.arange(-b1,c1-b1),[c0,1])
yy = np.tile(np.arange(-b0,c0-b0).reshape(c0,1),[1,c1])
rr = (xx**2 + yy**2)**(1/2)

# Primary component
engine_r = radii[wav][0]
engine = rr < engine_r
engine_c = plt.Circle(((b1-128)*pixel, (b0-128)*pixel), engine_r*pixel, color='w', linestyle='--', fill=False)

# Everything
everything_r = radii[wav][2]
everything = rr < everything_r
everything_c = plt.Circle(((b1-128)*pixel, (b0-128)*pixel), everything_r*pixel, color='w', fill=False, linestyle='--')

# Background
bg_r = radii[wav][3]
bg = (rr < bg_r) * (rr >= everything_r)
bg_c = plt.Circle(((b1-128)*pixel, (b0-128)*pixel), bg_r*pixel, color='w', fill=False, linestyle='-')

# Find secondary component
search = image * (1-engine)
nc = np.array(np.unravel_index(np.argmax(search),search.shape))
print('\nNorthern companion found at pixel: ', nc)

# Make secondary windows
xxnc = np.tile(np.arange(-nc[1],c1-nc[1]),[c0,1])
yync = np.tile(np.arange(-nc[0],c0-nc[0]).reshape(c0,1),[1,c1])
rrnc = (xxnc**2 + yync**2)**(1/2)

# Secondary component
northern_r = radii[wav][1]
northern = rrnc < northern_r
northern_c = plt.Circle((-(nc[0]-128)*pixel, -(nc[1]-128)*pixel), northern_r*pixel, color='w', linestyle='--', fill=False)

# Subtract background from image
background = np.sum(image*bg)/np.sum(bg)
image2 = np.array(image - background)
print("\nBackground of image: ", background)
print("Background of image2: ", np.sum(image2*bg)/np.sum(bg))

# Calculate fluxes
engine_flux = np.sum(image2 * engine)
northern_flux = np.sum(image2 * northern)
total_flux = np.sum(image2 * everything)
excess = total_flux - engine_flux - northern_flux

# Print results
print("\nThe ratios (%) using", wav,
"\nCentral: ", engine_flux/total_flux*100,
"\nNorthern: ", northern_flux/total_flux*100,
"\nExcess: ", excess/total_flux*100,
"\nSUM: ", (engine_flux + northern_flux + excess)/total_flux*100
)

# Plot
plt.figure(figsize=(8,7))
plt.imshow(np.flip(image.transpose()), origin='lower', extent=extent)
plt.gcf().gca().add_artist(engine_c)
plt.gcf().gca().add_artist(northern_c)
plt.gcf().gca().add_artist(everything_c)
plt.gcf().gca().add_artist(bg_c)
#plt.title('Image of ' + wav + num)
plt.xlabel('Relative RA (″)',fontsize=20)
plt.ylabel('Relative Dec (″)',fontsize=20)
plt.show()

## Test L
superimage = np.zeros([128,128])

for i in range(6):
    image = fits.getdata(path+'stacked_cube_L9_'+str(i+1)+'.fits')
    superimage += image

plt.figure(2)
plt.imshow(superimage, origin='lower')
plt.show()



## Get files
all_files = glob.glob(path+'stacked_cube'+'*'+'.fits')

image_files = {}

for wav in wavelengths:
    image_files[wav] = [file for file in all_files if (wav in file and ('_PSF' not in file))]

## Analysis
relative_fluxes = {}

for wav in wavelengths:
    print('\n----------', wav, '----------')
    output = np.zeros([len(image_files[wav]),3])
    
    for i in range(len(image_files[wav])):
        print('\n-----'+wav+'_'+str(i+1)+'-----')
        
        # Read data
        image = fits.getdata(image_files[wav][i])
        
        # Normalise image
        image /= np.max(image)
        
        # Get image dimensions
        c0, c1 = image.shape[0], image.shape[1]
        b0, b1 = int(c0/2), int(c1/2)
        
        # Make primary windows
        xx = np.tile(np.arange(-b1,c1-b1),[c0,1])
        yy = np.tile(np.arange(-b0,c0-b0).reshape(c0,1),[1,c1])
        rr = (xx**2 + yy**2)**(1/2)
        
        # Primary component
        engine_r = radii[wav][0]
        engine = rr < engine_r
        engine_c = plt.Circle((b1, b0), engine_r, color='r', fill=False)
        
        # Everything
        everything_r = radii[wav][2]
        everything = rr < everything_r
        everything_c = plt.Circle((b1, b0), everything_r, color='g', fill=False)
        
        # Background
        bg_r = radii[wav][3]
        bg = (rr < bg_r) * (rr >= everything_r)
        bg_c = plt.Circle((b1, b0), bg_r, color='g', fill=False)
        
        # Find secondary component
        search = image * (1-engine)
        nc = np.array(np.unravel_index(np.argmax(search),search.shape))
        print('\nNorthern companion found at pixel: ', nc)
        
        # Make secondary windows
        xxnc = np.tile(np.arange(-nc[1],c1-nc[1]),[c0,1])
        yync = np.tile(np.arange(-nc[0],c0-nc[0]).reshape(c0,1),[1,c1])
        rrnc = (xxnc**2 + yync**2)**(1/2)
        
        # Secondary component
        northern_r = radii[wav][1]
        northern = rrnc < northern_r
        northern_c = plt.Circle((nc[1], nc[0]), northern_r, color='r', fill=False)
        
        
        # Subtract background from image
        background = np.sum(image*bg)/np.sum(bg)
        image2 = np.array(image - background)
        #print("\nBackground of image: ", background)
        #print("Background of image2: ", np.sum(image2*bg)/np.sum(bg))
        
        # Calculate fluxes
        engine_flux = np.sum(image2 * engine)
        northern_flux = np.sum(image2 * northern)
        total_flux = np.sum(image2 * everything)
        excess = total_flux - engine_flux - northern_flux
        
        # Print results
        print(#"\nThe ratios (%) using", wav,
        "\nCentral: ", engine_flux/total_flux*100,
        "\nNorthern: ", northern_flux/total_flux*100,
        "\nExcess: ", excess/total_flux*100,
        #"\nSUM: ", (engine_flux + northern_flux + excess)/total_flux*100
        )
        
        print("\nOnly central and northern"
        "\nCentral: ", engine_flux/(engine_flux+northern_flux),
        "\nNorthern: ", northern_flux/(engine_flux+northern_flux),
        )
        
        # Output in fraction: central, northern, excess
        output[i,:] = [engine_flux/(engine_flux+northern_flux), northern_flux/(engine_flux+northern_flux), excess/(engine_flux+northern_flux)]
    
    relative_fluxes[wav] = output

## Summarise output
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
#'L9': 3.80
}

inv_microns = {
2.24: ['IB2.24'],
3.740: ['NB3.74'],
4.051: ['NB4.05'],
1.66: ['H7', 'H9'],
1.265: ['J7', 'J9'],
2.18: ['K7', 'K9'],
#3.80: ['L9']
}

df = pd.DataFrame.from_dict(relative_fluxes, orient='index')
df.to_csv(path+'Relative_NACO_flux.csv')

target_microns = []
target_ratios = []

for micron in list(inv_microns.keys()):
    target_microns.append(micron)
    sum_ratio = np.zeros([2])
    data_points = 0
    for wav in inv_microns[micron]:
        sum_ratio += np.sum(relative_fluxes[wav][:,0:2],0)
    sum_ratio /= np.sum(sum_ratio)
    target_ratios.append(sum_ratio)
    
