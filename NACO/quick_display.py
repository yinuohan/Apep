'''
A quick display of NACO images. 
'''


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

path = 'Where is the data?'

filenames = {
'stacked_cube_J7_1': 'J 7-hole',
'stacked_cube_H7_1': 'H 7-hole',
'stacked_cube_K7_1': 'Ks 7-hole' ,
'stacked_cube_J9_1': 'J 9-hole',
'stacked_cube_H9_1': 'H 9-hole',
'stacked_cube_K9_1': 'Ks 9-hole',
'stacked_cube_L9_1': 'L 9-hole',
'stacked_cube_IB2.24_2': 'K full frame',
'stacked_cube_NB3.74': 'L full frame',
}

clip_lower = 0
clip_upper = 1

nrows = 3
ncols = 3

low = -13.221 * 128/1000
high = 13.221 * 128/1000

fig, axs = plt.subplots(nrows,ncols,figsize=(7,7.5))
for i in range(len(filenames)):
    filename = list(filenames.keys())[i]
    name = filenames[filename]
    
    image = fits.getdata(path + filename + '.fits')
    
    if 'IB' in filename:
        lowdim = int(512/4*1)
        highdim = int(512/4*3)
        image = image[lowdim:highdim,lowdim:highdim]
        image /= np.max(image)
        image = np.clip(image, clip_lower, 0.5)
    elif 'NB' in filename:
        lowdim = int(512/8*3)
        highdim = int(512/8*5)
        image = image[lowdim:highdim,lowdim:highdim]
        image /= np.max(image)
        image = np.clip(image, clip_lower, 0.5)
    #elif 'L' in filename:
    #    new_image = np.zeros((256,256))
    #    new_image[64:192,64:192] = image
    #    image = np.array(new_image)
    #    image /= np.max(image)
    #    image = np.clip(image, clip_lower, clip_upper)
    else:
        image /= np.max(image)
        image = np.clip(image, clip_lower, clip_upper)
    
    
    axs[int(i/ncols), i%ncols].imshow(image, origin='lower', cmap='viridis', extent=[low,high,low,high])
    axs[int(i/ncols), i%ncols].set_title(name)
    #axs[int(i/ncols), i%ncols].

for ax in axs.flat:
    ax.label_outer()

plt.show()
