#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 13:09:20 2019

@author: Yinuo Han
"""

"""This code makes a flat file from images taken from VISIR, and saves it in the same directory. The file made here will be used in the VISIR_CALIBRATION.py file """

from astropy.io import fits
import numpy as np
import os
import matplotlib.pyplot as plt
#path to the directory where all the VISIR files are stored
path = "where is the data?"
#Store the files in a list
files = []
for filename in os.listdir(path):
    if filename.startswith('VISIR_IMG_OBS'):
        files.append(filename)
files.sort()
#print (files)
print (len(files))
#Taking the average of the sky frames
skies = []
for i in range(len(files)):
    hdulist = fits.open(path+files[i])
    upper = hdulist[1].data
    lower = hdulist[2].data
    sky = upper + lower
    y = upper.shape[0]
    x = upper.shape[1]
    print (str(y) + ' ' + str(x))
    if i % 2 == 0:
        sky[:,450:] = np.zeros((y,x-450))
    if i % 2 == 1:
        sky[:,:450] = np.zeros((y,450))
    skies.append(sky)

final_sky = sum(skies) / len(files)
#hdu = fits.PrimaryHDU(final_sky)
#hdu.writeto('Making_sky.fits')

#Getting calibrator images to make the dark file
darkfiles = []
for filename in os.listdir(path):
    if filename.startswith('VISIR_IMG_CAL'):
        darkfiles.append(filename)
darkfiles.sort()
#print (darkfiles)
print (len(darkfiles))

#Getting the dark file by taking the average of the calibrator files
darks = []
for i in range(len(files)):
    hdulist = fits.open(path+files[i])
    dark = hdulist[3].data
    darks.append(dark)

#Subtracting the darks from the skies to get the flat file
final_dark = sum(darks) / float(len(darkfiles))
final_flat = final_sky - final_dark
plt.imshow(final_dark)
#hdu = fits.PrimaryHDU(final_flat)
#hdu.writeto(path+'Making_flat.fits')
