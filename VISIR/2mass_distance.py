'''
Author: Yinuo Han

Detects the location of 2MASS J16004953-5142506 visible in the VISIR frames of Apep and calculats displacement
between two epochs for plate scale calibration. 
'''

import os
import math
import time
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.ndimage as img
import scipy.interpolate as intp


## Import data
path = 'Where is the data?'

VISIR_IMAGES = {
'2016': 'APEP_VISIR_2016_J8.9_large.fits',
'2017': '2nd_epoch_coord1.fits',
'2018': 'APEP_VISIR_2018_J8.9_large.fits',
}

ep1 = '2017'
ep2 = '2018'

img_1 = fits.getdata(path+VISIR_IMAGES[ep1])
img_2 = fits.getdata(path+VISIR_IMAGES[ep2])

img_1 = img_1/img_1.max()
img_2 = img_2/img_2.max()

## Find brightest pixels
# Gaussian blur
fwhm = 1
sd = fwhm/2.355

img_1_blur = img.gaussian_filter(img_1, sd)
img_2_blur = img.gaussian_filter(img_2, sd)

# Shift
cent_1 = np.unravel_index(img_1_blur.argmax(), img_1_blur.shape)
cent_2 = np.unravel_index(img_2_blur.argmax(), img_2_blur.shape)

if ep1 == '2017':
    dim = img_2.shape
    img_1 = img_1[cent_1[0]-cent_2[0]:cent_1[0]+(dim[0]-cent_2[0]), cent_1[1]-cent_2[1]:cent_1[1]+(dim[1]-cent_2[1])]
    img_1_blur = img.gaussian_filter(img_1, sd)

if ep2 == '2017':
    dim = img_1.shape
    img_2 = img_2[cent_2[0]-cent_1[0]:cent_2[0]+(dim[0]-cent_1[0]), cent_2[1]-cent_1[1]:cent_2[1]+(dim[1]-cent_1[1])]
    img_2_blur = img.gaussian_filter(img_2, sd)

# Search range
cent_xl = 160
cent_yl = 220
cent_range = 60

star_xl = 350
star_yl = 110
star_range = 50


# Epoch 1
print(ep1)

ind_cent_1 = img_1_blur[cent_yl:cent_yl+cent_range, cent_xl:cent_xl+cent_range].argmax()

cent_1 = np.unravel_index(ind_cent_1, img_1_blur[cent_yl:cent_yl+cent_range, cent_xl:cent_xl+cent_range].shape) + np.array([cent_yl, cent_xl])

print('Centre', cent_1)

ind_star_1 = img_1_blur[star_yl:star_yl+star_range, star_xl:star_xl+star_range].argmax()

star_1 = np.unravel_index(ind_star_1, img_1_blur[star_yl:star_yl+star_range, star_xl:star_xl+star_range].shape) + np.array([star_yl, star_xl])

print('2mass', star_1)


# Epoch 2
print(ep2)
ind_cent_2 = img_2_blur[cent_yl:cent_yl+cent_range, cent_xl:cent_xl+cent_range].argmax()

cent_2 = np.unravel_index(ind_cent_2, img_2_blur[cent_yl:cent_yl+cent_range, cent_xl:cent_xl+cent_range].shape) + np.array([cent_yl, cent_xl])

print('Centre', cent_2)

ind_star_2 = img_2_blur[star_yl:star_yl+star_range, star_xl:star_xl+star_range].argmax()

star_2 = np.unravel_index(ind_star_2, img_2_blur[star_yl:star_yl+star_range, star_xl:star_xl+star_range].shape) + np.array([star_yl, star_xl])

print('2mass', star_2)


## Plot
scale = 0.045

plt.figure(1)
plt.imshow(np.clip(img_1_blur,0,0.05), origin='lower', cmap='hot', extent=[-(430-250)*scale, (430-180)*scale, -(430-180)*scale, (430-250)*scale])
#plt.scatter((cent_1[1]-180)*scale,(cent_1[0]-250)*scale, c='b', marker='+', s=100)
#plt.scatter((star_1[1]-180)*scale,(star_1[0]-250)*scale, c='g', marker='+', s=100)
plt.title("Epoch 1")
plt.xlabel('Relative RA (″)')
plt.ylabel('Relative Dec (″)')


plt.figure(2)
plt.imshow(np.clip(img_2_blur,0,0.03), origin='lower', cmap='hot', extent=[-(430-250)*scale, (430-180)*scale, -(430-180)*scale, (430-250)*scale])
#plt.scatter((cent_2[1]-250)*scale,(cent_2[0]-180)*scale, c='b', marker='+', s=100)
#plt.scatter((star_2[1]-250)*scale,(star_2[0]-180)*scale, c='g', marker='+', s=100)
plt.title("Epoch 2")
plt.xlabel('Relative RA (″)')
plt.ylabel('Relative Dec (″)')

'''
plt.figure(3)
plt.imshow(np.clip(img.gaussian_filter(img_2-img_2_blur-img_1+img_1_blur,5),-0.005,0.005), origin='lower')
plt.title('2018 - 2016 low pass filtered + blurred')
plt.scatter(cent_2[1],cent_2[0], c='r', marker='+', s=100)
plt.scatter(star_2[1],star_2[0], c='g', marker='+', s=100)

plt.figure(4)
plt.imshow(img.gaussian_filter(img_2-img_2_blur,5), origin='lower')
plt.title('2018 low pass filtered blurred')
plt.scatter(cent_2[1],cent_2[0], c='r', marker='+', s=100)
plt.scatter(star_2[1],star_2[0], c='g', marker='+', s=100)

plt.figure(5)
plt.imshow(img.gaussian_filter(img_1-img_1_blur,5), origin='lower')
plt.title('2016 low pass filtered blurred')
plt.scatter(cent_2[1],cent_2[0], c='r', marker='+', s=100)
plt.scatter(star_2[1],star_2[0], c='g', marker='+', s=100)

plt.figure(6)
plt.imshow(img_1_blur - img_2_blur, origin='lower')
plt.title('Difference of blurred image')
'''


plt.show()

dist_1 = np.linalg.norm(cent_1 - star_1)
dist_2 = np.linalg.norm(cent_2 - star_2)

print('The percentage change is (%)')
print((dist_2 - dist_1) / dist_1 * 100 )

