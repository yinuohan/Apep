#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 13:06:27 2019

@author: Yinuo Han
"""

"""This file takes all the raw images from the VISIR instrument
 and spits out two final 360x360 images at wavelengths 8.9 micrometres,
 and 11.7 micrometres, using the chopping and nodding technique, 
 for images of dimensions 900x1024. There were 48 images taken at 11.7 micrometres
 and 40 images taken at 8.9 micrometres"""

import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# Settings
#wavelength = 'J8.9'
wavelength = 'Q3'

add_Q3_first = 1

#S pecify path to the directory where the VISIR files are stored
path = "/import/silo4/snert/TINTAGEL_DATA/VISIR/2018_05/raw/"
savedir = "/import/silo4/snert/TINTAGEL_DATA/VISIR/2018_05/Analysis/"

prefixes = {'J8.9': 'VISIR_IMG_OBS141',
            'Q3': 'VISIR_IMG_OBS156',
            }

# Open and store all the VISIR image files
files = []

for filename in os.listdir(path):
    if filename.startswith(prefixes[wavelength]):
        files.append(filename)
    
files.sort()
#print (files)

'''
#opening the flat file that was made in Making_flats.py
flat_file = fits.open(path+"Making_flat.fits")
f = flat_file[0].data

flat_file.close()
'''

# Chopping and nodding the images
images = []

for i in range(int(len(files)/2)):
    a_file = fits.open(path + files[2 * i])
    b_file = fits.open(path + files[2 * i + 1])
    
    print(a_file[0].header['HIERARCH ESO INS FILT1 NAME'])
    print(b_file[0].header['HIERARCH ESO INS FILT1 NAME'])
    
    nodA = a_file[3].data
    nodB = b_file[3].data
    
    a_file.close()
    b_file.close()
    
    image = (nodA - nodB)#/f
    images.append(image)

# Individual frames low signal to noise
if wavelength == 'Q3' and add_Q3_first == 1:
    new_images = []
    new_images.append(np.sum(images[0::2],0))
    new_images.append(np.sum(images[1::2],0))
    images = list(new_images)

# After chopping and nodding the image needs to be extracted and re-centred
z = len(images)
y = len(images[0])
#x = len(images[0][0])

# Dimensions of the final image
#dim = 430
dim = 450
#d_left = 180
d_left = 225
d_right = dim - d_left
#d_top = 180
d_top = 225
d_bottom = dim - d_top

# Array for the positive and negatives of the chopped and nodded images
pos = np.zeros((2*z,dim,dim))
neg = np.zeros((2*z,dim,dim))

"""Finding the brightest and darkest spots and cutting a 360x360 image out
   and storing them in pos and neg"""
   
for k in range(z):
    print(k)
    image = images[k]
   
    # Consider only the area with the image and it's negative on the left half
    if wavelength == 'Q3':
        left = 200
        right = 300
        up = 700
        down = 200
    else:
        left = 150
        right = 400
        up = y
        down = 0
    
    left_slice = image[down:up,left:right]
    
    indmax = np.argmax(left_slice)
    c_max = np.unravel_index(indmax,left_slice.shape)
    indmin = np.argmin(left_slice)
    c_min = np.unravel_index(indmin,left_slice.shape)
    
    print(c_max)
    print(c_min)
    
    pos_image = image[(down+c_max[0])-d_bottom:(down+c_max[0])+d_top,(left+c_max[1])-d_left:(left+c_max[1])+d_right]
    neg_image = image[(down+c_min[0])-d_bottom:(down+c_min[0])+d_top,(left+c_min[1])-d_left:(left+c_min[1])+d_right]
    
    pos[2*k] = pos_image
    neg[2*k] = neg_image
    
    # Consider only the area with the images and it's negative on the right half
    if wavelength == 'Q3':
        left = 550
        right = 650
        up = 700
        down = 200
    else:
        left = 500
        right = 750
        up = y
        down = 0
    
    right_slice = image[down:up,left:right]
    
    indmax = np.argmax(right_slice)
    c_max = np.unravel_index(indmax,right_slice.shape)
    indmin = np.argmin(right_slice)
    c_min = np.unravel_index(indmin,right_slice.shape)
    
    print(c_max)
    print(c_min)
    
    pos_image = image[(down+c_max[0])-d_bottom:(down+c_max[0])+d_top,(left+c_max[1])-d_left:(left+c_max[1])+d_right]
    neg_image = image[(down+c_min[0])-d_bottom:(down+c_min[0])+d_top,(left+c_min[1])-d_left:(left+c_min[1])+d_right]
    
    pos[2*k+1] = pos_image
    neg[2*k+1] = neg_image
    
# Final image = positive image - negative image
final_image = np.mean(pos,0) - np.mean(neg,0)

# Plot final image
plt.imshow(final_image, origin='lower')

# Save final image
hdu1 = fits.PrimaryHDU(final_image)
hdu1.writeto(savedir+'APEP_VISIR_2018_'+wavelength+'_huge2.fits')

