'''
Author: Yinuo Han

Corrects for shifts of each NACO frame, 
windows components in each frame with a super-Gaussian,
stores files according to calibration pairs. 

See settings below for details. 

Warning: may alter directories in file system. 

Warning: new .FITS files are named to hack pipeline. 
Please keep track of your components!
'''


import numpy as np
import os
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import ndimage
from statistics import mode 
from scipy.signal import general_gaussian
from shutil import copyfile

# Define working directories
workingdir = './'

# Define files
#apep_files = [1]
apep_files = ['0' + str(i+1) for i in range(6)] + [i+49 for i in range(6)]
calibrator_files = ['0' + str(i+7) for i in range(6)] + [i+55 for i in range(6)]

# Define binary and companion star diameters
bdim = 45
cdim = 35

# Define roughly where to expect companion star
xrmin, xrmax = 0, np.inf
yrmin, yrmax = 0, np.inf

# Define which window to use
use_pos_window = 0

# Define window aggressiveness
super_gaussian_p = 3

# Define whether to window calibrator star
window_calibrator = 0


# Auto from here --------------------

# Initialise
print('Current directory: ' + os.getcwd())
os.chdir(workingdir)
print('\ncurrent directory: ' + os.getcwd())

apep_files = [glob.glob('./cube*'+str(file)+'.fits')[0] for file in apep_files]
calibrator_files = [glob.glob('./cube*'+str(file)+'.fits')[0] for file in calibrator_files]
                    
apep_files = np.sort(apep_files)
calibrator_files = np.sort(calibrator_files)

print('\nApep files: ')
print(apep_files)
print('\nCalibrator files: ')
print(calibrator_files)

bfiles = []
cfiles = []
bposfiles = []
cposfiles = []

cal_bfiles = []
cal_cfiles = []
cal_bposfiles = []
cal_cposfiles = []

# Go through each file
for n in range(len(apep_files)):
    file = apep_files[n]
    
    print('\n\nFile name: ' + file)
    #print('\nList of peripheral star centres')
    
    img = fits.getdata(file)
    z0, y0, x0 = img.shape
    
    bmin, bmax = int(y0/2 - bdim/2), int(y0/2 + bdim/2)
    
    # Find centre of binary and companion star
    central_bright = np.empty([z0,2],dtype='int')
    periph_bright = np.empty([z0,2],dtype='int')
    
    for z in range(z0):
        blur_img = ndimage.gaussian_filter(img[z], sigma=1)
        central_bright[z,:] = np.array(np.unravel_index(blur_img[bmin:bmax,bmin:bmax].argmax(),blur_img[bmin:bmax,bmin:bmax].shape)) + [bmin, bmin]        
        
        periph_img = np.array(blur_img)
        periph_img[bmin:bmax,bmin:bmax] = np.min(periph_img)
        periph_bright[z,:] = np.array(np.unravel_index(periph_img.argmax(),periph_img.shape))
        
        #print(central_bright[z,:])
        #print('Frame ' + str(z) + ': ' + str(periph_bright[z,:]))
        
        # Shift if frame is flipped
        if not (yrmin < periph_bright[z,0] < yrmax and xrmin < periph_bright[z,1] < xrmax):
            print('Frame requiring shift: ' + str(z))
            
            # Shift brightest peripheral pixel to centre
            img[z] = np.roll(img[z],int(y0/2)-periph_bright[z,0],axis=0)
            img[z] = np.roll(img[z],int(x0/2)-periph_bright[z,1],axis=1)
            
            periph_bright[z,0] = central_bright[z,0] + int(y0/2)-periph_bright[z,0]
            periph_bright[z,1] = central_bright[z,1] + int(x0/2)-periph_bright[z,1]
            
            central_bright[z,0] = int(y0/2)
            central_bright[z,1] = int(x0/2)
            
            print('Frame ' + str(z) + ' new: ' + str(periph_bright[z,:]))
            
            #plt.figure(z)
            #plt.imshow(img[z])
    
    
    mpy = np.median(periph_bright[:,0])
    mpx = np.median(periph_bright[:,1])
    print('\nMedian peripheral brightest pixel: ' + str([mpy, mpx]))
    
    
    # Delete bad frames
    periph_norm = np.linalg.norm(np.array([mpy,mpx])-periph_bright,axis=1)
    bad_frames = np.where(periph_norm > 5)
    img = np.delete(img, bad_frames, 0)
    print('\nBad frames removed: ' + str(bad_frames))
    
    # Prepare to window
    bimg = np.array(img) # Only binary left
    cimg = np.array(img) # Only companion left
    
    bimgpos = np.array(img)
    cimgpos = np.array(img)
    
    # Make windows
    x = np.arange(0, x0, 1)
    x = x - np.median(x)
    y = np.arange(0, y0, 1)
    y = y - np.median(y)
    xx, yy = np.meshgrid(x, y)
    
    rr2 = (xx**2 + yy**2)
    P = super_gaussian_p
    
    bsigma = bdim/2.2
    csigma = cdim/2.2
    
    bwinpos_centre = np.exp(-(rr2/2/bsigma**2)**P) # Centred at centre of image
    cwinpos_centre = np.exp(-(rr2/2/csigma**2)**P)
    
    # Shift windows to right location - same shift across all frames
    bwinpos = np.array(bwinpos_centre) # To keep binary
    cwinpos = np.array(cwinpos_centre) # To keep companion star
    cwinpos = np.roll(cwinpos,int(mpx-x0/2),axis=1) 
    cwinpos = np.roll(cwinpos,int(mpy-y0/2),axis=0)
    
    bwin = 1-cwinpos # To keep binary
    cwin = 1-bwinpos # To keep companion star
    
    
    
    # Apply 2 types of windows
    bimg = bwin * img
    cimg = cwin * img
        
    bimgpos = bwinpos * img
    cimgpos = cwinpos * img
    
    cimg = np.roll(cimg,int(x0/2-mpx),axis=2)
    cimg = np.roll(cimg,int(y0/2-mpy),axis=1)
    
    cimgpos = np.roll(cimgpos,int(x0/2-mpx),axis=2)
    cimgpos = np.roll(cimgpos,int(y0/2-mpy),axis=1)
    
    # Store Apep images
    if use_pos_window == 0:
        bfiles.append(bimg)
        cfiles.append(cimg)
    elif use_pos_window == 1:
        bposfiles.append(bimgpos)
        cposfiles.append(cimgpos)
    
        
    # Apply same windows to calibrator star
    if window_calibrator == 1:
        calfile = calibrator_files[n]
        print('\n\nCalibration file name: ' + file)
        
        cal = fits.getdata(calfile)
        
        # Positive windows can just reuse the centered windows for img of Apep
        bcalpos = cal * bwinpos_centre
        ccalpos = cal * cwinpos_centre
        
        # Negative window for calibrating companion star needs shifting
        bcalwin = np.array(bwin)
        ccalwin = np.array(cwin) # To keep companion star
        ccalwin = np.roll(ccalwin,int(-mpx+x0/2),axis=1) 
        ccalwin = np.roll(ccalwin,int(-mpy+y0/2),axis=0)
        
        bcal = cal * bcalwin
        ccal = cal * ccalwin
        
        # Store calibrator images
        if use_pos_window == 0:
            cal_bfiles.append(bcal)
            cal_cfiles.append(ccal)
        elif use_pos_window == 1:
            cal_bposfiles.append(bcalpos)
            cal_cposfiles.append(ccalpos)
             

# Diagnostic plots
plt.figure(1)
plt.imshow(bwin*cwin)

plt.figure(2)
plt.imshow(bwin*cwin*img[0])

plt.figure(3)
plt.imshow((1-bwin*cwin)*img[0])


# Make directories and save .fits with disguised names
maindir = os.getcwd()
    
binary_vs_companion = maindir + '/binary_vs_companion' # No.0
companion_vs_calibrator = maindir + '/companion_vs_calibrator' # No.1
binary_vs_calibrator = maindir + '/binary_vs_calibrator' # No.2

subdirs = [binary_vs_companion, companion_vs_calibrator, binary_vs_calibrator]


# Copy dcubes, .idlvar and .script
idlvar = glob.glob(maindir + '/cubeinfo*.idlvar')[0]
script = glob.glob(maindir + '/*_cube.script')[0]
dcubes = glob.glob(maindir + '/dcube*.fits')

for n in range(len(subdirs)):
    if use_pos_window == 0:
        subdirs[n] += '_neg'
    elif use_pos_window == 1:
        subdirs[n] += '_pos'
    
    try:
        os.mkdir(subdirs[n])
    except FileExistsError:
        pass
    
    copyfile(idlvar, subdirs[n] + '/' + idlvar.split('/')[-1])
    copyfile(script, subdirs[n] + '/' + script.split('/')[-1].split('_cube')[0] + '_bis_cal.script')
    for dcube in dcubes:
        copyfile(dcube, subdirs[n] + '/' + dcube.split('/')[-1])

print('\nDirectories ready')

# Initialise saving pramaeters
if use_pos_window == 0:
    finalbfiles = bfiles
    finalcfiles = cfiles
    cal_finalbfiles = cal_bfiles
    cal_finalcfiles = cal_cfiles
elif use_pos_window == 1:
    finalbfiles = bposfiles
    finalcfiles = cposfiles
    cal_finalbfiles = cal_bposfiles
    cal_finalcfiles = cal_cposfiles
else:
    print('Variable use_pos_window must be 0 or 1!')
    
    
# Save calibrator star 
if window_calibrator == 0:
    for i in range(len(calibrator_files)):
        copyfile(maindir + calibrator_files[i][1:], subdirs[1] + calibrator_files[i][1:])
        copyfile(maindir + calibrator_files[i][1:], subdirs[2] + calibrator_files[i][1:])
elif window_calibrator == 1:
    for i in range(len(calibrator_files)):
        hdu = fits.PrimaryHDU(cal_finalbfiles[i])
        hdu.writeto(subdirs[2] + calibrator_files[i][1:])
    
        hdu = fits.PrimaryHDU(cal_finalcfiles[i])
        hdu.writeto(subdirs[1] + calibrator_files[i][1:])
    
print('\nCalibrators done')

# Save Apep
for i in range(len(apep_files)):
    hdu = fits.PrimaryHDU(finalbfiles[i])
    hdu.writeto(subdirs[0] + apep_files[i][1:])
    hdu.writeto(subdirs[2] + apep_files[i][1:])
    
    hdu = fits.PrimaryHDU(finalcfiles[i])
    hdu.writeto(subdirs[0] + calibrator_files[i][1:])
    hdu.writeto(subdirs[1] + apep_files[i][1:])

print('\nApep done')

