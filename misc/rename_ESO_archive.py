#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 12:51:29 2019

@author: Yinuo Han

This program renames ESO's archived name of each 
individual .FITS file to their original names. 
"""

import os
from astropy.io import fits

owd = os.getcwd()
print ('\n' + 'Orignal working directory: ' + owd)

os.chdir('Where is the data?')

cwd = os.getcwd()
print ('New working directory: ' + cwd + '\n')

for oldname in os.listdir():
    if oldname.startswith('NACO.') and oldname.endswith('.fits'):
        hdulist = fits.open(oldname)
        for i in range(0,100):
            newname = str(hdulist[0].header[i])
            if newname.startswith('NACO_'):
                break
        hdulist.close()
        
        print (oldname)
        print (newname)
        print()
            
        namefile = open('Names.txt', 'a')
        namefile.write(oldname + '\n' + newname + '\n\n')
        namefile.close()
        
        os.rename (oldname, newname)
    
    if oldname.startswith('VISIR.') and oldname.endswith('.fits'):
        hdulist = fits.open(oldname)
        for i in range(0,100):
            newname = str(hdulist[0].header[i])
            if newname.startswith('VISIR_'):
                break
        hdulist.close()
        
        print (oldname)
        print (newname)
        print()
        
        namefile = open('Names.txt', 'a')
        namefile.write(oldname + '\n' + newname + '\n\n')
        namefile.close()
        
        os.rename (oldname, newname)