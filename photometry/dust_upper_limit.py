'''
A quick estimate for an upper limit of the flux from
dust in the VISIR images of Apep. 
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
invert = 1

colourmap = 'viridis'

lowerlim = -0.00075
upperlim = 0.00075

down = 210
up = 300
left = 210
right = 300


## Configuration
path = 'Where is the data?'


## Read data
image = fits.getdata(path + 'stacked_cube.fits')

# Normalise image
image /= np.max(image)

# Confirm brightest pixel at centre
c0, c1 = np.shape(image)[0], np.shape(image)[1]
print("Image dimensions: ", (c0, c1))

b0, b1 = np.unravel_index(image.argmax(), image.shape)[0], np.unravel_index(image.argmax(), image.shape)[1]
print("Brightest pixel at: ", (b0, b1))

clipped = np.clip(image,lowerlim,upperlim)

## Plot image
fig, ax = plt.subplots(0)
plt.imshow(clipped*invert, origin='lower', cmap=colourmap)
plt.title('Stacked cube')
plt.show()


## Analysis
# Calculate centre flux
centre_flux = np.sum(image[down:up,left:right])
print("Centre (original): ", centre_flux)

# Get rid of centre
background = np.array(clipped)
background[down:up,left:right] = 0
plt.subplots(1)
plt.imshow(background*invert, origin='lower', cmap=colourmap)
plt.title('Background')
plt.show()

# Calculate background flux
background_flux = np.sum(background)
print("Background (clipped): ", background_flux)

'''
# Plot pixels
fig, ax = plt.subplots(0)
#plt.plot(np.mean(clipped,1))
plt.plot(background)
plt.plot([0,512],[0,0])
plt.title('Clipped background')
plt.show()
'''

xx = np.tile(np.arange(-b1,c1-b1),[c0,1])
yy = np.tile(np.arange(-b0,c0-b0).reshape(c0,1),[1,c1])
rr = (xx**2 + yy**2)**(1/2)

spiral = (rr < 500) * (rr >= 480)

#plt.plot()

# No more than 5% for spiral total flux in 4.05. 
# Using stacked cube. 

# No more than 10% for spiral total flux in 3.74. 
# Using stacked cube. 

# No more than 10% for spiral total flux in 2.24. 
# Using stacked cube2. 

