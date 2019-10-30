'''
A quick plot to compare two VISIR images. 
'''

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.ndimage as img

## Settings
scale = 0.045

plt.rcParams.update({'font.size': 14})


## Import data
path = 'where is the data?'

VISIR_IMAGES = {
'2016': 'APEP_VISIR_2016_J8.9_large.fits',
'2017': '2nd_epoch_coord1.fits',
'2018J8.9': 'APEP_VISIR_2018_J8.9_huge_460.fits',
'2018Q3': 'APEP_VISIR_2018_Q3_huge.fits',
}

ep1 = '2018J8.9'
ep2 = '2018Q3'

img_1 = fits.getdata(path+VISIR_IMAGES[ep1])
img_2 = fits.getdata(path+VISIR_IMAGES[ep2])

img_1 = img_1/img_1.max()
img_2 = img_2/img_2.max()

img_1[372,100] = 0
img_1 = img_1[5:-5,5:-5]

img_1_blur = img.gaussian_filter(img_1,0)
img_2_blur = img.gaussian_filter(img_2,0)

c1 = int(img_1.shape[0]/2)
c2 = int(img_2.shape[0]/2)


## Plot
plt.subplots(figsize=(17,7))
plt.subplot(121)
plt.imshow(np.clip(img_1_blur,0,0.03), origin='lower', cmap='hot', extent=[-c1*scale, c1*scale, -c1*scale, c1*scale])
#plt.title("J8.9 (2018)")
plt.xlabel('Relative RA (″)',fontsize=20)
plt.ylabel('Relative Dec (″)',fontsize=20)

plt.subplot(122)
plt.imshow(np.clip(img_2_blur,0,0.5), origin='lower', cmap='hot', extent=[-c2*scale, c2*scale, -c2*scale, c2*scale])
#plt.title("Q3 (2018)")
plt.xlabel('Relative RA (″)',fontsize=20)
plt.ylabel('Relative Dec (″)',fontsize=20)

plt.show()

