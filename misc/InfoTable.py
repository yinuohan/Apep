#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 14:29:49 2019

@author: yhan9672

Extracts useful header information from .FITS files
in a given directory and outputs a .CSV file. 
"""

import os
from astropy.io import fits
import numpy as np
import csv

owd = os.getcwd()
print ('\n' + 'Orignal working directory: ' + owd)

#Comment this line when run on Linux

cwd = os.getcwd()
print ('New working directory: ' + cwd + '\n')

#par = ['ORIGFILE', 'TELESCOP', 'INSTRUME', 'DATE-OBS', 'RA', 'DEC', 'NAXIS1', 'NAXIS2', 'HIERARCH ESO INS RETA2 ROT', 'HIERARCH ESO INS OPTI1 ID', 'HIERARCH ESO INS OPTI2 ID', 'HIERARCH ESO INS OPTI3 ID', 'HIERARCH ESO INS OPTI4 ID', 'HIERARCH ESO INS OPTI5 ID', 'HIERARCH ESO INS OPTI6 ID', 'HIERARCH ESO INS OPTI7 ID', 'PARANG END', 'HIERARCH ESO DET DIT', 'HIERARCH ESO TEL ALT']

# OOBSNAME changed to TELESCOP and INSTRUME
# HIERARCH ESO INS RETA2 ROT = Angle relative to Wollaston
# PARANG END = Parallactic angle at end
# HIERARCH ESO DET DIT = Integration time
# HIERARCH ESO TEL ALT = Alt angle at start
# 
# File#ID  ObsBlk  TIME  RA  DEC  Filter  OPTI4ID  Mask  AX1  AX3  T_int  ALT  PsAng  ROT


Filename = ['ORIGFILE']
Type = ['HIERARCH ESO DPR TYPE']
Telescope = ['TELESCOP']
Instrument = ['INSTRUME']
Date = ['DATE-OBS']
RA = ['RA']
Dec = ['DEC']
IntTime = ['HIERARCH ESO DET DIT']
X_axis = ['NAXIS1']
Y_axis = ['NAXIS2']
Frames = ['NAXIS3']
AngWol = ['HIERARCH ESO INS RETA2 ROT']
OneID = ['HIERARCH ESO INS OPTI1 ID']
TwoID = ['HIERARCH ESO INS OPTI2 ID']
ThreeID = ['HIERARCH ESO INS OPTI3 ID']
FourID = ['HIERARCH ESO INS OPTI4 ID']
FiveID = ['HIERARCH ESO INS OPTI5 ID']
SixID = ['HIERARCH ESO INS OPTI6 ID']
SevenID = ['HIERARCH ESO INS OPTI7 ID']
ParAngEnd = ['PARANG END']
AltAngStart = ['HIERARCH ESO TEL ALT']

par = [Filename, Type, Telescope, Instrument, Date, RA, Dec, IntTime, X_axis, Y_axis, Frames, AngWol, OneID, TwoID, ThreeID, FourID, FiveID, SixID, SevenID, ParAngEnd, AltAngStart]


for filename in os.listdir():
    if (filename.startswith('NACO_') or filename.startswith('VISIR_')) and filename.endswith('.fits'):
        
        hdulist = fits.open(filename)
        
        for parname in par:
            while True:
                try:
                    parname.append(hdulist[0].header[parname[0]])
                    break
                except KeyError:
                    parname.append('')
                    break
        
        hdulist.close()



with open('Info.csv', 'w', newline='') as info:
    writer = csv.writer(info)
    for i in range(0,len(par[0])):
        row = []
        for parname in par:
            row.append(parname[i])
        writer.writerow(row)
