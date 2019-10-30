'''
Author: Yinuo Han

This is a watered down version of ridge-crawler. 
The program draws the skeleton of an input image of one 
epoch only without calculating radial displacements. 
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
compare = '2018'

new_step = 5
step_factor = new_step/3
#step_factor = 1

plot_sample_points = 0

plot_slices = 0

background_epoch = 1

invert = -1

colourmap = 'gray'

textcolour = 'k'


## Configuration
path = 'Where is the data?'

VISIR_IMAGES = {
'2016': 'APEP_VISIR_2016_J8.9_large.fits',
#'2017': 'second_epoch_visJ8.9.fits',
'2017': '2nd_epoch_coord1.fits',
'2018': 'APEP_VISIR_2018_J8.9_large.fits',
}

dates = {
'2016': date(2016,8,13),
'2017': date(2017,7,31),
'2018': date(2018,5,21),
}

cutoff = {
'2016': 0.02,
'2017': 0.02,
'2018': 0.008,
}

epoch1_filename = VISIR_IMAGES[compare]

# 'ridge name': image, startxy, starttheta, step, nstepmax, width
ridge_params_full = {
'a2016': [[119,145.5], 180, 3, 47, 14],
'a2017': [[115,144], 180, 3, 47, 14],
'a2018': [[116,145.5], 180, 3, 65, 16],#max 55

'b2016': [[130,220], 50, 3, 40, 14],#max 65
'b2017': [[127,218], 50, 3, 20, 14],#40
'b2018': [[117,209], 50, 3, 84, 14],

'c2016': [[251,251], -40, 3, 30, 16],
'c2017': [[251,251], -40, 3, 30, 16],
'c2018': [[160,179], -130, 3, 80, 16],

'd2016': [[305,143], -90, 3, 30, 14],
'd2017': [[304,141], -90, 3, 30, 14],
'd2018': [[236,45], 45, 3, 100, 16],

'e2016': [[132,77], -10, 3, 35, 16],#max 60
'e2017': [[138,70], -10, 3, 35, 14],
'e2018': [[138,70], -10, 3, 35, 14],

'f2016': [[42,157], -60, 3, 30, 14],
'f2017': [[41.5,157], -60, 3, 30, 14],
'f2018': [[41,157], -60, 3, 30, 14],

'g2018': [[237,282], -60, 3, 30, 14],
} 


## Read data
epoch1 = fits.getdata(path + epoch1_filename)

# Normalise images
epoch1 /= np.max(epoch1)

# Confirm brightest pixel at centre
print("Image dimensions: " + str(np.shape(epoch1)))

c0 = 180
c1 = 180

b1 = np.unravel_index(epoch1.argmax(), epoch1.shape)
print("Brightest pixels at: " + str(b1))

if b1 != (c0, c1):
    epoch1 = epoch1[b1[0]-c0:b1[0]+c0, b1[1]-c1:b1[1]+c1]
    
print("Adjusted brightest pixels at: " + 
str(np.unravel_index(epoch1.argmax(), epoch1.shape)))


## High pass filter
# Low pass filter
fwhm = 20
sd = fwhm/2.355
epoch1_blur = img.gaussian_filter(epoch1, sd)

# High pass filter
epoch1_sharp = epoch1 - epoch1_blur

# Clip
epoch1_sharp_chopped = np.clip(epoch1_sharp, 0, cutoff[compare])

# Blur
fwhm = 5
sd = fwhm/2.355
epoch_1_skeleton = img.gaussian_filter(epoch1_sharp_chopped, sd)


## Find ridges
def ridgecrawler(ridgename, image, startxy, starttheta, step, nstepmax, w, plot_sample_points=0, plot_slices=0, pixel=1):
    
    # Initialise parameters
    theta = math.radians(starttheta)
    point = startxy
    rx = [startxy[0]]
    ry = [startxy[1]]
    
    
    # Start the loop
    for s in range(nstepmax):
        
        # Plot current point
        if plot_sample_points == 1:
            plt.scatter(point[0], point[1], marker='x', c='r', s=50)
            
        # Draw line
        tpoint = point + step * np.array([math.cos(theta), math.sin(theta)])
        
        linex = (np.array(range(w))-w/2+0.5) * math.sin(theta) + tpoint[0]
        liney = -(np.array(range(w))-w/2+0.5) * math.cos(theta) + tpoint[1]
        
        
        # Find profile along line
        slice = image[np.round(liney).astype(int), np.round(linex).astype(int)]
        
        if plot_slices == 1:
            plt.plot(linex,liney,c='darkseagreen')
        #print(point)
        
        # Find max point location on line
        pmax = slice.argmax()
        if pmax in [0, len(slice)-1]:
            print('\nSlice not wide enough for ' + ridgename)
            print('Exited at step ' + str(s))
            break;
        
        # Turning point is 1D coordinate along the line
        # Calculated using quadratic fit to 3 highest points
        tp = -(slice[pmax+1] - slice[pmax-1]) / 2 / (slice[pmax-1] - 2*slice[pmax] + slice[pmax+1])
        
    
        # Calculate new point and new theta
        newpoint = [linex[pmax] + math.cos(theta) * tp, liney[pmax] + math.sin(theta) * tp]
        newtheta = np.arctan2(newpoint[1] - point[1], newpoint[0] - point[0])
        
        # Update point and theta
        point = newpoint
        theta = newtheta
        
        # Add to record
        rx = rx + [point[0]]
        ry = ry + [point[1]]
        
    ridgeline = np.array([rx, ry])
    return ridgeline


## Fit splines
def fit_spline(ridgename1, rl):
    # Points
    x = rl[ridgename1][0]
    y = rl[ridgename1][1]
    
    # Fit spline
    tck, u = intp.splprep([x, y], k=3, s=1)
    unew = np.arange(0, 1.01, 1/len(x)*0.01)
    out = intp.splev(unew, tck)
    
    # Plot fitted spline
    colours = ['darkseagreen', 'lightcoral']
    plt.plot(out[0], out[1], color=colours[1])
    
    # Return as array
    out = np.array(out)
    
    return out
    

## Plot original image
fig, ax = plt.subplots()
plt.imshow(epoch_1_skeleton*invert, origin='lower',cmap=colourmap)
plt.show()

## Apply ridgecrawler and fit_spline
rl = {}
rl_fitted = {}

ridges = ['a','b','c','d','e','f','g']

sketch = np.zeros(epoch_1_skeleton.shape)

for name in ridges:
    ridge_name = name+compare
    
    params = ridge_params_full[ridge_name]
    
    rl[name] = ridgecrawler(
    name, 
    epoch_1_skeleton, 
    params[0], 
    params[1], 
    np.round(params[2]*step_factor), 
    int(params[3]/step_factor), 
    params[4], 
    plot_sample_points, plot_slices)
    
    rl_fitted[name] = fit_spline(name, rl)
    
    ys = np.rint(rl_fitted[name][1]).astype(int)
    xs = np.rint(rl_fitted[name][0]).astype(int)
    
    sketch[ys,xs] = 1

#fig, ax = plt.subplots()
#plt.imshow(sketch, origin='lower',cmap=colourmap)
#plt.show()


