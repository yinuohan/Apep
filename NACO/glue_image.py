#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 11:49:10 2019

@author: Yinuo Han

Pieces all NACO frames together into a large glued image.
"""

import os
import csv
import pandas as pd
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.ndimage as ndimage

# Settings
wavelength = 'NB3.74'
#wavelength = 'NB4.05'

# Configurations
path = "/import/silo4/snert/TINTAGEL_DATA/NACO/2016_04/raw/"
flatpath = "/import/silo4/snert/TINTAGEL_DATA/NACO/2016_04/Testing/Flats/"

prefixes = {'NB3.74': 'NACO_IMG_SCI119*',
            'NB4.05': 'NACO_IMG_SCI119*',
            }

suffixes = {'NB3.74': ['00' + str(i+13) for i in range(12)],
            'NB4.05': ['00' + str(i+1) for i in range(12)],
            }
            
badpixelnames = {'NB3.74': 'bad_pixels_3.74.csv',
                'NB4.05': 'bad_pixels_4.05.csv',
               }
               
savedirs = {'NB3.74': '/import/silo4/snert/TINTAGEL_DATA/NACO/2016_04/Analysis/Apep/NB_3.74/',
            'NB4.05': '/import/silo4/snert/TINTAGEL_DATA/NACO/2016_04/Analysis/Apep/NB_4.05/',
            }

# Find all image files
files = [glob.glob(path + prefixes[wavelength] + str(suffix) + '.fits')[0] for suffix in suffixes[wavelength]]

# Get sky and bad pixels
savedir = savedirs[wavelength]
supersky = fits.getdata(savedir+'supersky.fits')
c0, c1 = supersky.shape[0], supersky.shape[1]

badarray = np.array(pd.read_csv(flatpath+badpixelnames[wavelength])).flatten()
badarray = np.unravel_index(badarray,[c0,c1])
bad_pixels = np.zeros([c0,c1])
bad_pixels[badarray] = 1

# Glue onto superimage
superimage = np.zeros([2*c0,2*c1])
glue = np.zeros([2*c0,2*c1])

fwhm = 5

up = int(3*c0/4)
down = int(c0/4)
left = int(c1/4)
right = int(3*c1/4)

for file in files:
    image = fits.getdata(file)
    
    for i in range(len(image)):
        frame = image[i]
        cleanframe = frame - supersky
        #cleanframe[0:int(c0/2),0:int(c1/2)] = np.median(cleanframe)
        cleanframe[badarray] = np.median(cleanframe)
        '''change this to interpolate around'''
        
        blurframe = ndimage.gaussian_filter(cleanframe,fwhm/2.355)
        cmax = np.unravel_index(np.argmax(blurframe[down:up,left:right]),blurframe[down:up,left:right].shape)
        cmax = np.array(cmax) + [down, left]
        
        superimage[c0-cmax[0]:2*c0-cmax[0],c1-cmax[1]:2*c1-cmax[1]] += cleanframe
        glue[c0-cmax[0]:2*c0-cmax[0],c1-cmax[1]:2*c1-cmax[1]] += 1
        
        print(i)
    
final_image = superimage/glue
plt.imshow(final_image,origin='lower')

fits.PrimaryHDU(final_image).writeto(savedir+'glued_image.fits')
fits.PrimaryHDU(glue).writeto(savedir+'glue.fits')

clipped = np.clip(np.nan_to_num(final_image),0,5000)
fits.PrimaryHDU(clipped).writeto(savedir+'clipped.fits')


