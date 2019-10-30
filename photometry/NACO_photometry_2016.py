'''
Author: Yinuo Han

Performs aperture photometry on the 2016 epoch of 
NACO images of Apep. 
'''


import os
import math
import time
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.ndimage as img
import scipy.interpolate as intp
from datetime import date


## Settings
calculate = '2016_4'

invert = 1

colourmap = 'viridis'

plt.close('all')


## Configuration
path = 'Where is the data?'

VISIR_IMAGES = {
'2016_1': 'stacked_cube_2.24_1.fits',
'2016_2': 'stacked_cube_2.24_2.fits',
'2016_3': 'stacked_cube_3.74.fits',
'2016_4': 'stacked_cube_4.05.fits',
}

cutoff = {
'2016_1': 0.02,
'2016_2': 0.02,
'2016_3': 0.005,
'2016_4': 0.005,
}

radii = {
'2016_1': [17, 12],
'2016_2': [17, 12],
'2016_3': [15, 7],
'2016_4': [15, 7],
}

## Read data
image = fits.getdata(path + VISIR_IMAGES[calculate])

# Normalise image
image /= np.max(image)

# Confirm brightest pixel at centre
c0, c1 = np.shape(image)[0], np.shape(image)[1]
print("Image dimensions: ", (c0, c1))

b0, b1 = np.unravel_index(image.argmax(), image.shape)[0], np.unravel_index(image.argmax(), image.shape)[1]
print("Brightest pixel at: ", (b0, b1))

'''
if (b0, b1) != (c0/2, c1/2):
    image = image[b0[0]-c0:b0[0]+c0, b1[1]-c1:b1[1]+c1]
    print("Adjusted brightest pixels at: " + 
str(np.unravel_index(image.argmax(), image.shape)))
'''

## Plot image
low = -13.221 * 128/1000
high = 13.221 * 128/1000
extent = [low, high, low, high]

fig, ax = plt.subplots(0)
#plt.imshow(image*invert, origin='lower', extent=[-lim,lim,-lim,lim],cmap=colourmap)
plt.imshow(np.clip(image,0,1)*invert, origin='lower', cmap=colourmap)
plt.title('Image of ' + calculate)
#plt.xlabel('Relative RA (″)')
#plt.ylabel('Relative Dec (″)')
plt.show()


## Analysis
# Make primary windows
xx = np.tile(np.arange(-b1,c1-b1),[c0,1])
yy = np.tile(np.arange(-b0,c0-b0).reshape(c0,1),[1,c1])
rr = (xx**2 + yy**2)**(1/2)

'''
[Name] = grid of region with 1s and 0s
[Name]_r = radius of region
[Name]_c = circle surrounding region for plotting
'''

# Primary component
engine_r = radii[calculate][0]
engine = rr < engine_r
engine_c = plt.Circle((b1, b0), engine_r, color='r', fill=False)

# Everything
everything_r = 100
everything = rr < everything_r
everything_c = plt.Circle((b1, b0), everything_r, color='g', fill=False)

# Background
bg_r = 170
bg = (rr < bg_r) * (rr >= everything_r)
bg_c = plt.Circle((b1, b0), bg_r, color='g', fill=False)

# Find secondary component
search = image * (1-engine)
nc = np.array(np.unravel_index(np.argmax(search),search.shape))
print('Northern companion found at pixel: ', nc)

# Make secondary windows
xxnc = np.tile(np.arange(-nc[1],c1-nc[1]),[c0,1])
yync = np.tile(np.arange(-nc[0],c0-nc[0]).reshape(c0,1),[1,c1])
rrnc = (xxnc**2 + yync**2)**(1/2)

# Secondary component
northern_r = radii[calculate][1]
northern = rrnc < northern_r
northern_c = plt.Circle((nc[1], nc[0]), northern_r, color='r', fill=False)


## Calculate fluxes
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


print("\nThe ratios (%) using", calculate,
"\nCentral: ", engine_flux/total_flux*100,
"\nNorthern: ", northern_flux/total_flux*100,
"\nExcess: ", excess/total_flux*100,
"\nSUM: ", (engine_flux + northern_flux + excess)/total_flux*100,
"\n\nOnly central and northern"
"\nCentral: ", engine_flux/(engine_flux+northern_flux),
"\nNorthern: ", northern_flux/(engine_flux+northern_flux),
)


## Plot things
# Plot circles
plt.imshow(engine, alpha = 0, origin='lower')
plt.gcf().gca().add_artist(engine_c)

plt.imshow(northern, alpha = 0, origin='lower')
plt.gcf().gca().add_artist(northern_c)

plt.imshow(bg, alpha = 0, origin='lower')
plt.gcf().gca().add_artist(bg_c)

plt.imshow(everything, alpha = 0, origin='lower')
plt.gcf().gca().add_artist(everything_c)

# Plot northern companion
plt.scatter(nc[1], nc[0], marker='x', color='r')
plt.show()


# Plot background
plt.subplots(0)
image3 = image * (1-engine-northern)
plt.imshow(np.clip(image3,0,cutoff[calculate]), origin='lower')
plt.show()
