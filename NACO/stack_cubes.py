#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 11:49:10 2019

@author: Yinuo Han

Stacks all frames in each NACO cube into one frame. 
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
#wavelength = 'H7'
#wavelength = 'H9'
#wavelength = 'J7'
#wavelength = 'J9'
#wavelength = 'K7'
#wavelength = 'K9'
#wavelength = 'L9'

#stack_PSF = 0 # 0 to stack target, 1 to stack PSF

def stack(wavelength, stack_PSF):
    # Configurations
    xrmin, xrmax = 0, 105
    yrmin, yrmax = 0, np.inf
    r = 25
    
    suffixes_target = {'H7': [111, 112, '0005', '0006'],
                       'H9': [214, 215],
                       'J7': [113, 114, '0007', '0008'],
                       'J9': [216, 217],
                       'K7': [109, 110, '0003', '0004'],
                       'K9': [212, 213],
                       'L9': [65, 66, 67, 68, 69, 70],
                       'M9': [],
                      }
    
    suffixes_PSF = {'H7': [117, 118, '0011', '0012'],
                    'H9': [220, 221],
                    'J7': [119, 120, '0013', '0014'],
                    'J9': [222, 223],
                    'K7': [115, 116, '0009', '0010'],
                    'K9': [218, 219],
                    'L9': [77, 78, 79, 80, 81, 82],
                    'M9': [],
                   }
    
    rootdir = '/import/silo4/snert/TINTAGEL_DATA/NACO/2019_03/Analysis/Apep/'
    savedirs = {'H7': rootdir + 'H_7_holes/',
                'H9': rootdir + 'H_9_holes/',
                'J7': rootdir + 'J_7_holes/',
                'J9': rootdir + 'J_9_holes/',
                'K7': rootdir + 'Ks_7_holes/',
                'K9': rootdir + 'Ks_9_holes/',
                'L9': rootdir + 'L_prime_9_holes/',
                'M9': rootdir + 'M_prime_9_holes/',
               }
    
    # Initialise
    suffixes = 0
    filename_suffix = 0
    
    if stack_PSF == 0:
        suffixes = suffixes_target
        filename_suffix = 'stacked_cube_' + str(wavelength)
    elif stack_PSF == 1:
        suffixes = suffixes_PSF
        filename_suffix = 'stacked_cube_PSF_'+  str(wavelength)
    
    # Find all image files
    files = [glob.glob(savedirs[wavelength] + 'cube*' + str(suffix) + '.fits')[0] for suffix in suffixes[wavelength]]
    
    # Get dimensions
    sample = fits.getdata(files[0])[0]
    c0, c1 = sample.shape[0], sample.shape[1]
    
    i = 1
    
    # Stack
    for file in files:
        image = fits.getdata(file)
        
        needshift = 0
        fine = 0
        nframes = 0
    
        # Glue onto superimage
        superimage = np.zeros([c0,c1])
    
        for frame in image:
            # Shift brightest pixel to centre
            blurframe = ndimage.gaussian_filter(frame, sigma=1)
            cmax = np.unravel_index(np.argmax(blurframe),blurframe.shape)
            
            newframe = np.roll(frame,int(c0/2-cmax[0]),axis=0)
            newframe = np.roll(newframe,int(c1/2-cmax[1]),axis=1)
            
            blurframe = ndimage.gaussian_filter(newframe, sigma=1)
            
            if  stack_PSF == 0:
                # Find brightest pixel not near centre
                pframe = np.array(blurframe)
                pframe[int(c0/2)-r:int(c0/2)+r,int(c0/2)-r:int(c0/2)+r] = np.min(pframe)
                pmax = np.array(np.unravel_index(pframe.argmax(),pframe.shape))
                
                # Shift if frame flipped
                if not (yrmin < pmax[0] < yrmax and xrmin < pmax[1] < xrmax):
                    # Shift brightest peripheral pixel to centre
                    newframe = np.roll(newframe,int(c0/2)-pmax[0],axis=0)
                    newframe = np.roll(newframe,int(c1/2)-pmax[1],axis=1)
                    
                    #print('Frame requiring shift: ' + str(nframes))
                    needshift += 1
                else:
                    #print('Frame FINE!')
                    fine += 1
            
            #plt.figure(nframes)
            #plt.imshow(newframe)
            #plt.title(str(nframes) + str(cmax) + str(pmax))
            
            superimage += newframe
            nframes += 1
        
        superimage /= nframes
        print('Shifted ' + str(needshift) + '/' + str(needshift + fine) + ' frames')
        plt.imshow(superimage,origin='lower')
        #fits.PrimaryHDU(superimage).writeto(rootdir+filename_suffix+'_'+str(i)+'.fits')
        i += 1

wavelengths = ['H7', 'H9', 'J7', 'J9', 'K7', 'K9', 'L9']
for wavelength in wavelengths:
    print(wavelength)
    for stack_PSF in [0, 1]:
        print(stack_PSF)
        stack(wavelength, stack_PSF)

        
        
        