'''
Makes NACO cube into a gif. 
'''


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import imageio

path = 'Where is the data?'

filename = 'cube0001'

cube = fits.getdata(path+filename+'.fits')

plt.imshow(np.sum(cube,1))

imageio.mimsave(path+filename+'.gif', cube)

plt.show()

