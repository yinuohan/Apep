'''
Author: Yinuo Han

Performs aperture photometry on any epoch of 
VISIR image of Apep. 
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
calculate = '2018J1'

invert = 1

colourmap = 'viridis'

subtract_northern_background = 1

fix_2017 = 1

plt.close('all')

plt.rcParams.update({'font.size': 18})


## Configuration
path = 'Where is the data?'

VISIR_IMAGES = {
'2016J1': 'APEP_VISIR_2016_J8.9_large.fits',
'2016J2': 'APEP_VISIR_2016_B11.7_large.fits',
'2017J1': 'second_epoch_visJ8.9.fits',
'2018J1': 'APEP_VISIR_2018_J8.9_large.fits',
'2018Q3': 'APEP_VISIR_2018_Q3_large.fits',
}

cutoff = {
'2016J1': 0.1,#0.02,
'2016J2': 0.1,#Something
'2017J1': 0.012,
'2018J1': 0.05,
'2018Q3': 1,#Something
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
pixel = 45 / 1000 # seconds / pixel

left = 180
right = 430-180
up = 180
down = 430-180

extent = [-left*pixel, right*pixel, -down*pixel, up*pixel]

fig, ax = plt.subplots(0,figsize=(8,7))
plt.imshow(np.clip(image,0,cutoff[calculate])*invert, origin='lower', cmap=colourmap, extent=extent)
#plt.title('Image of ' + calculate)
plt.xlabel('Relative RA (″)',fontsize=20)
plt.ylabel('Relative Dec (″)',fontsize=20)
plt.yticks([-10,-5,0,5])
plt.show()


## Analysis
# Make central windows
xx = np.tile(np.arange(-b1,c1-b1),[c0,1])
yy = np.tile(np.arange(-b0,c0-b0).reshape(c0,1),[1,c1])
rr = (xx**2 + yy**2)**(1/2)

'''
[Name] = grid of region with 1s and 0s
[Name]_r = radius of region
[Name]_c = circle surrounding region for plotting
'''

# Central engine
engine_r = 10
engine = rr < engine_r
engine_c = plt.Circle((b1, b0), engine_r, color='r', fill=False)

# Dust
dust_r = 150
dust = (rr < dust_r) * (rr >= engine_r)
dust_c = plt.Circle(((b1-left)*pixel, (b0-down)*pixel), dust_r*pixel, color='w', fill=False, linestyle='--')

# Dust background ( = everything background)
dustbg_r = 170
dustbg = (rr < dustbg_r) * (rr >= dust_r)
dustbg_c = plt.Circle(((b1-left)*pixel, (b0-down)*pixel), dustbg_r*pixel, color='w', fill=False)

# Everything
everything = rr < dust_r

# Shoot lines
nlines = 1000
length = dust_r

theta_bin = np.pi*2 / nlines
r_bin = 1

theta = np.arange(0, np.pi*2, theta_bin)
r = np.arange(0,length, r_bin)

lines = np.zeros([nlines,length])
linex = np.round((np.reshape(np.cos(theta),[len(theta),1]) * r + b1)).astype(int)
liney = np.round((np.reshape(np.sin(theta),[len(theta),1]) * r + b0)).astype(int)

radial = np.mean(image[liney,linex], axis = 0)

# Make azimuthally averaged image
radial2 = np.zeros(max(c0,c1))
radial2[0:len(radial)] = radial
imaverage = radial2[np.round(rr).astype(int)]
'''
imaverage = np.array(image)

for i in range(imaverage.shape[1]):
    for j in range(imaverage.shape[0]):
        index = int(np.round(rr[j,i]))
        if index < len(radial):
            imaverage[j,i] = radial[index]
'''

# Find northern companion
diff = image - imaverage

searchx = np.array([182,189]) - [180,180] + b1
searchy = np.array([193,200]) - [180,180] + b0
search = diff[searchy[0]:searchy[1],searchx[0]:searchx[1]]
nc = np.array(np.unravel_index(np.argmax(search),search.shape)) + [searchy[0], searchx[0]]
print('Northern companion found at pixel: ', nc)

# Make northern windows
xxnc = np.tile(np.arange(-nc[1],c1-nc[1]),[c0,1])
yync = np.tile(np.arange(-nc[0],c0-nc[0]).reshape(c0,1),[1,c1])
rrnc = (xxnc**2 + yync**2)**(1/2)

# Northern companion
northern_r = 5
northern = rrnc < northern_r
northern_c = plt.Circle(((nc[1]-left)*pixel, (nc[0]-down)*pixel), northern_r*pixel, color='w', fill=False, linestyle='--')

# Northern background
northernbg_r = 7
northernbg = (rrnc < northernbg_r) * (rrnc >= northern_r)
northernbg_c = plt.Circle(((nc[1]-left)*pixel, (nc[0]-down)*pixel), northernbg_r*pixel, color='w', fill=False)

# Find faraway star
searchx2 = np.array([195,205]) + b1 # 200 from centre
searchy2 = np.array([-123,-113]) + b0 # -118 from centre
search2 = image[searchy2[0]:searchy2[1],searchx2[0]:searchx2[1]]
search2 = img.gaussian_filter(search2,5/2.3)
fs = np.array(np.unravel_index(np.argmax(search2),search2.shape)) + [searchy2[0], searchx2[0]]
print('Faraway star found at pixel: ', fs)

# Make faraway windows
xxfs = np.tile(np.arange(-fs[1],c1-fs[1]),[c0,1])
yyfs = np.tile(np.arange(-fs[0],c0-fs[0]).reshape(c0,1),[1,c1])
rrfs = (xxfs**2 + yyfs**2)**(1/2)

# Faraway star
faraway_r = 7
faraway = rrfs < faraway_r
faraway_c = plt.Circle(((fs[1]-left)*pixel, (fs[0]-down)*pixel), faraway_r*pixel, color='w', fill=False, linestyle='--')
((b1-left)*pixel, (b0-down)*pixel), dustbg_r*pixel

# Faraway background
farawaybg_r = 20
farawaybg = (rrfs < farawaybg_r) * (rrfs >= faraway_r)
farawaybg_c = plt.Circle(((fs[1]-left)*pixel, (fs[0]-right)*pixel), farawaybg_r*pixel, color='w', fill=False)


## Calculate fluxes
# Subtract background from image
background = np.sum(image*dustbg)/np.sum(dustbg)
image2 = np.array(image - background)
print("\nBackground of image: ", background)
print("Background of image2: ", np.sum(image2*dustbg)/np.sum(dustbg))

# Do the same for nothern companion
northern_background = np.sum(diff*northernbg)/np.sum(northernbg)
diff2 = np.array(diff - northern_background)
print("\nBackground of diff: ", northern_background)
print("Background of diff2: ", np.sum(diff2*northernbg)/np.sum(northernbg))

# Do the same for faraway star
faraway_background = np.sum(image*farawaybg)/np.sum(farawaybg)
image3 = np.array(image - faraway_background)
print("\nBackground of image at faraway: ", faraway_background)
print("Background of image3: ", np.sum(image3*farawaybg)/np.sum(farawaybg))

# Calculate fluxes
if subtract_northern_background == 1: # diff or diff2?
    northern_flux = np.sum(diff2 * northern)
else:
    northern_flux = np.sum(diff * northern)
dust_flux = np.sum(image2 * dust) - northern_flux
engine_flux = np.sum(image2 * engine)
faraway_flux = np.sum(image3 * faraway)
if c1-b1 <= 200 and fix_2017 == 1:
    faraway_flux = faraway_flux * 2
total_flux = np.sum(image2 * everything) + faraway_flux


print("\nThe ratios (%) using", calculate,
"\nNorthern companion: ", northern_flux/total_flux*100,
"\nCentral engine: ", engine_flux/total_flux*100,
"\nDust: ", dust_flux/total_flux*100,
"\nDistant star: ", faraway_flux/total_flux*100,
"\nSUM: ", (northern_flux + engine_flux + dust_flux + faraway_flux) / total_flux * 100,
)

## Plot things
# Plot circles
#plt.imshow(dust, alpha = 0, origin='lower')
plt.gcf().gca().add_artist(dust_c)

#plt.imshow(dustbg, alpha = 0, origin='lower')
plt.gcf().gca().add_artist(dustbg_c)

#plt.imshow(engine, alpha = 0, origin='lower')
#plt.gcf().gca().add_artist(engine_c)

#plt.imshow(engine, alpha = 0, origin='lower')
#plt.gcf().gca().add_artist(northern_c)


# Plot lines
#plt.plot(linex.transpose(), liney.transpose())

# Plot northern companion
plt.scatter((nc[1]-left)*pixel, (nc[0]-down)*pixel, marker='+', color='crimson', s=50)
plt.show()

# Plot faraway search
#plt.scatter(searchx2, searchy2)
#plt.imshow(faraway, alpha = 0, origin='lower')
plt.gcf().gca().add_artist(faraway_c)
#plt.imshow(farawaybg, alpha = 0, origin='lower')
plt.gcf().gca().add_artist(farawaybg_c)
plt.show()


'''
# Plot image
plt.subplots(1)
plt.imshow(image, origin='lower')
plt.imshow(engine, alpha = 0, origin='lower')
plt.gcf().gca().add_artist(engine_c)
plt.show()

# Plot azimuthally averaged image
plt.subplots(1)
plt.imshow(imaverage, origin='lower')
plt.show()
'''

# Plot difference
plt.subplots(1,figsize=(8,7))
plt.imshow(np.clip(diff,-0.02,0.02), origin='lower', extent=extent)

# Plot northern search
#plt.scatter(searchx, searchy)
#plt.imshow(northern, alpha = 0, origin='lower')
plt.gcf().gca().add_artist(northern_c)
#plt.imshow(northernbg, alpha = 0, origin='lower')
plt.gcf().gca().add_artist(northernbg_c)

plt.xlabel('Relative RA (″)',fontsize=20)
plt.ylabel('Relative Dec (″)',fontsize=20)

plt.xlim([-1.2, 1.2])
plt.ylim([-0.8, 1.6])


plt.show()

'''
Clearly elongated in the east-west direction! Is the same seen in all other images? Does it move across the 3 epochs?
'''

'''
# Plot radial distribution
plt.subplots(1)
plt.plot(radial)
plt.show()
'''
