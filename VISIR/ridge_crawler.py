'''
Author: Yinuo Han

This program takes 2 images as input and outputs
the radial displacements of ridges. 
Output displacements are fitted sub-pixel. 

Users must specify the initial location and direction
of each ridge to be analysed. 

A range of settings are offered. Check settings and
configurations for details. 

'''


import os
import math
import time
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.ndimage as img
import scipy.interpolate as intp
from datetime import date

## Settings
compare = ['2016', '2018']

use_sem = 2

sem_factor = 4 # How many independent data points are there in each ridge?

new_step = 5
step_factor = new_step/3
#step_factor = 1

plot_sample_points = 0

plot_slices = 0

background_epoch = 1

invert = -1

colourmap = 'gray'

textcolour = 'k'

distance_in_kpc = 2.4

plt.rcParams.update({'font.size': 18})

apply_lim = 0

save_rc = 0

save_diff = 0

figsuffix = '_n'

## Configuration
path = 'Where is the data?'

VISIR_IMAGES = {
'2016': 'APEP_VISIR_2016_J8.9_large.fits',
'2017': '2nd_epoch_coord1.fits',
'2018': 'APEP_VISIR_2018_J8.9_large.fits',
}

dates = {
'2016': date(2016,8,13),
'2017': date(2017,7,31),
'2018': date(2018,5,21),
}

'''
cutoff = {
'2016': 0.02,
'2017': 0.02,
'2018': 0.008,
}

'''
# For difference image
cutoff = {
'2016': 0.01,
'2017': 0.012,
'2018': 0.006,
}

epoch1_filename = VISIR_IMAGES[compare[0]]
epoch2_filename = VISIR_IMAGES[compare[1]]

# 'ridge name': image, startxy, starttheta, step, nstepmax, width
ridge_params_full = {
'a2016': [[120,145.5], 180, 3, 47, 16],
'a2017': [[115,144], 180, 3, 47, 14],
'a2018': [[116,145.5], 180, 3, 47, 14],#max 55

'b2016': [[130,220], 50, 3, 40, 14],#max 65
'b2017': [[127,218], 50, 3, 20, 14],#40
'b2018': [[127,218], 50, 3, 40, 14],

'c2016': [[251,251], -40, 3, 30, 16],
'c2017': [[251,251], -40, 3, 30, 16],
'c2018': [[251,251], -40, 3, 30, 16],

'd2016': [[305,143], -90, 3, 30, 14],
'd2017': [[304,141], -90, 3, 30, 14],
'd2018': [[307,138], -90, 3, 30, 14],

'e2016': [[132,77], -10, 3, 35, 16],#max 60
'e2017': [[138,70], -10, 3, 35, 14],
'e2018': [[138,70], -10, 3, 35, 14],

'f2016': [[42,157], -60, 3, 30, 14],
'f2017': [[41.5,157], -60, 3, 30, 14],
'f2018': [[41,157], -60, 3, 30, 14],
} 


## Read data
epoch1 = fits.getdata(path + epoch1_filename)
epoch2 = fits.getdata(path + epoch2_filename)

# Normalise images
epoch1 /= np.max(epoch1)
epoch2 /= np.max(epoch2)

# Confirm brightest pixel at centre
print("Image dimensions: " + str(np.shape(epoch1)) + str(np.shape(epoch2)))

#c0 = int(np.shape(epoch1)[0] / 2)
#c1 = int(np.shape(epoch1)[1] / 2)
c0 = 180
c1 = 180

b1 = np.unravel_index(epoch1.argmax(), epoch1.shape)
b2 = np.unravel_index(epoch2.argmax(), epoch2.shape)

print("Brightest pixels at: " + str(b1) + str(b2))

if b1 != (c0, c1):
    epoch1 = epoch1[b1[0]-c0:b1[0]+c0, b1[1]-c1:b1[1]+c1]
if b2 != (c0, c1):
    epoch2 = epoch2[b2[0]-c0:b2[0]+c0, b2[1]-c1:b2[1]+c1]
    
print("Adjusted brightest pixels at: " + 
str(np.unravel_index(epoch1.argmax(), epoch1.shape)) + str(np.unravel_index(epoch2.argmax(), epoch2.shape)))

# Plot images
'''
plt.figure()
fig1 = plt.imshow(epoch1, origin='lower')
plt.show()

plt.figure()
fig2 = plt.imshow(epoch2, origin='lower')
plt.show()
'''


## High pass filter for dislay
fwhm = 20
sd = fwhm/2.355

epoch1_blur = img.gaussian_filter(epoch1, sd)
epoch1_sharp = epoch1 - epoch1_blur

epoch2_blur = img.gaussian_filter(epoch2, sd)
epoch2_sharp = epoch2 - epoch2_blur

# Plot filtered images for display
'''
plt.figure()
plt.imshow(epoch1_blur, origin='lower')
plt.show()

plt.figure()
plt.imshow(epoch1_sharp, origin='lower')
plt.show()

plt.figure()
plt.imshow(epoch2_blur, origin='lower')
plt.show()

plt.figure()
plt.imshow(epoch2_sharp, origin='lower')
plt.show()
'''

# High pass filter for analysis
fwhm = 5
sd = fwhm/2.355

epoch1_sharp_chopped = np.clip(epoch1_sharp, 0, cutoff[compare[0]]) #change this max cutoff
epoch_1_skeleton = img.gaussian_filter(epoch1_sharp_chopped, sd)

epoch2_sharp_chopped = np.clip(epoch2_sharp, 0, cutoff[compare[1]]) #change this max cutoff
epoch_2_skeleton = img.gaussian_filter(epoch2_sharp_chopped, sd)

# Plot filtered images for analysis
'''
plt.figure()
plt.imshow(epoch_1_skeleton, origin='lower')
plt.show()

plt.figure()
plt.imshow(epoch_2_skeleton, origin='lower')
plt.show()
# '''





## Find ridges
def ridgecrawler(ridgename, image, startxy, starttheta, step, nstepmax, w, plot_sample_points=0, plot_slices=0, pixel=1):
    
    # Initialise parameters
    theta = math.radians(starttheta)
    point = startxy
    rx = [startxy[0]]
    ry = [startxy[1]]
    
    
    # Start the loop
    for s in range(nstepmax):
        
        # Plot current point
        if plot_sample_points == 1:
            plt.scatter(-(point[0]*pixel-lim), point[1]*pixel-lim, marker='x', c='r', s=50)
            
        # Draw line
        tpoint = point + step * np.array([math.cos(theta), math.sin(theta)])
        
        linex = (np.array(range(w))-w/2+0.5) * math.sin(theta) + tpoint[0]
        liney = -(np.array(range(w))-w/2+0.5) * math.cos(theta) + tpoint[1]
        
        
        # Find profile along line
        slice = image[np.round(liney).astype(int), np.round(linex).astype(int)]
        
        if plot_slices == 1:
            plt.plot(-(linex*pixel-lim),liney*pixel-lim,c='darkseagreen')
        #print(point)
        
        # Find max point location on line
        pmax = slice.argmax()
        if pmax in [0, len(slice)-1]:
            print('\nSlice not wide enough for ' + ridgename)
            print('Exited at step ' + str(s))
            break;
        
        # Turning point is 1D coordinate along the line
        # Calculated using quadratic fit to 3 highest points
        tp = -(slice[pmax+1] - slice[pmax-1]) / 2 / (slice[pmax-1] - 2*slice[pmax] + slice[pmax+1])
        
    
        # Calculate new point and new theta
        newpoint = [linex[pmax] + math.cos(theta) * tp, liney[pmax] + math.sin(theta) * tp]
        newtheta = np.arctan2(newpoint[1] - point[1], newpoint[0] - point[0])
        
        #print(newtheta*180/np.pi)
        #print()
        
        # What does this do?
        #if (point[1] - newpoint[1]) == 0.:
        #    newtheta = newtheta + np.pi()
        #    print("Adjusted")
        
        # Update point and theta
        point = newpoint
        theta = newtheta
        
        # Add to record
        rx = rx + [point[0]]
        ry = ry + [point[1]]
        
    ridgeline = np.array([rx, ry])
    return ridgeline
        


## Fit splines
def fit_spline(ridgename1, rl):
    # Points
    x = rl[ridgename1][0]
    y = rl[ridgename1][1]
    
    # Fit spline
    tck, u = intp.splprep([x, y], k=3, s=1)
    unew = np.arange(0, 1.01, 1/len(x)*0.01)
    out = intp.splev(unew, tck)
    
    # Plot fitted spline
    colours = ['darkseagreen', 'lightcoral']
    plt.plot(-(out[0]*pixel-lim), out[1]*pixel-lim, color=colours[int(ridgename1[-1])-1])
    
    for i in range(0, len(out[0])):
        #plt.plot([180, out[0][i]], [180, out[1][i]])
        pass
    
    # Return as array
    out = np.array(out)
    
    return out
    
    
    
## Find displacement
def ridge_displacement(ridgename1, ridgename2, rl_fitted, pixel=1, centre=[c0,c1]):
    # Convert to polar coordinates
    cx = centre[1]
    cy = centre[0]
    
    rf1_x, rf1_y = rl_fitted[ridgename1][0], rl_fitted[ridgename1][1]
    rf2_x, rf2_y = rl_fitted[ridgename2][0], rl_fitted[ridgename2][1]
    
    rf1_r = ((rf1_x-cx)**2 + (rf1_y-cy)**2)**0.5
    rf1_theta = np.mod(np.arctan2(rf1_y-cy, rf1_x-cx), 2*np.pi)
    
    rf2_r = ((rf2_x-cx)**2 + (rf2_y-cy)**2)**0.5
    rf2_theta = np.mod(np.arctan2(rf2_y-cy, rf2_x-cx), 2*np.pi)
    
    # Scanning parameters
    #scan_r_max = max(np.max(rf1_r), np.max(rf2_r))*1.2
    #bin_r = 0.1
    #scan_r = np.arange(0, scan_r_max, bin_r)
    
    min_theta = max(np.min(rf1_theta), np.min(rf2_theta))
    max_theta = min(np.max(rf1_theta), np.max(rf2_theta))
    bin_theta = 1*np.pi/180
    scan_theta = np.arange(min_theta, max_theta, bin_theta)
    
    # Scan for intersections
    index_1 = []
    index_2 = []
    for theta in scan_theta:
        i1 = np.argmin(np.abs(rf1_theta-theta))
        i2 = np.argmin(np.abs(rf2_theta-theta))
        
        index_1 += [i1]
        index_2 += [i2]
        
        # Plot distance
        plt.plot([-(rf1_x[i1]*pixel-lim), -(rf2_x[i2]*pixel-lim)], [rf1_y[i1]*pixel-lim, rf2_y[i2]*pixel-lim],color='k',alpha=0.2)
        #plt.plot([180, rf2_x[i2]], [180, rf2_y[i2]])
    
    # Add text to plot
    if not apply_lim:
        plt.text(-(np.median(rl[ridgename1], axis=1)[0]*pixel-lim-0.2),np.median(rl[ridgename1], axis=1)[1]*pixel-lim+0.7,ridgename1[0],color=textcolour)
        
    dist = np.array(rf2_r[index_2] - rf1_r[index_1])
    sel1 = np.array([rf1_x[index_1], rf1_y[index_1]])
    sel2 = np.array([rf2_x[index_2], rf2_y[index_2]])
    
    return dist, sel1, sel2, scan_theta



## Plot
def plot_ridges(epoch_1_name, epoch_2_name, epoch_1_skeleton, epoch_2_skeleton, rl, rl_fitted, rl_fitted_selected):
    #plt.figure()
    #plt.imshow(epoch_1_skeleton, origin='lower')
    
    plt.plot(rl_fitted[epoch_1_name][0]*pixel-lim, rl_fitted[epoch_1_name][1]*pixel-lim, 'b-')
    plt.plot(rl_fitted[epoch_2_name][0]*pixel-lim, rl_fitted[epoch_2_name][1]*pixel-lim, 'g--')
    
    plt.plot(rl_fitted_selected[epoch_1_name][0]*pixel-lim, rl_fitted_selected[epoch_1_name][1]*pixel-lim, 'b*')
    plt.plot(rl_fitted_selected[epoch_2_name][0]*pixel-lim, rl_fitted_selected[epoch_2_name][1]*pixel-lim, 'g*')
    
    plt.title(epoch_1_name)
    
    plt.show()





## Parameters
all_ridge_names = ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'e1', 'e2', 'f1', 'f2']
epoch_1_ridge_names = all_ridge_names[0:len(all_ridge_names):2]
epoch_2_ridge_names = all_ridge_names[1:len(all_ridge_names):2]

paired = {
'a': ['a1', 'a2'],
'b': ['b1', 'b2'],
'c': ['c1', 'c2'],
'd': ['d1', 'd2'],
'e': ['e1', 'e2'],
'f': ['f1', 'f2'],
}

all_pairs = ['a', 'b', 'c', 'd', 'e', 'f']

ridge_params = {}
for letter in paired.keys():
    ridge_params[letter+'1'] = [epoch_1_skeleton] + ridge_params_full[letter+compare[0]]
    ridge_params[letter+'2'] = [epoch_2_skeleton] + ridge_params_full[letter+compare[1]]


## Plot original image
pixel = 45 / 1000 # seconds / pixel

axis = np.array(range(len(epoch_1_skeleton)))
axis = (axis - np.median(axis)) * pixel

lim = len(epoch_1_skeleton) / 2 * pixel

plt.subplots(figsize=(8,7))

if background_epoch == 1:
    plt.imshow(epoch_1_skeleton*invert, origin='lower', extent=[lim,-lim,-lim,lim],cmap=colourmap)
elif background_epoch == 2:
    plt.imshow(epoch_2_skeleton*invert, origin='lower', extent=[lim,-lim,-lim,lim],cmap=colourmap)

#plt.title('Displacement from '+compare[0]+' to '+compare[1])
plt.xlabel('Relative RA (″)',fontsize=20)
plt.ylabel('Relative Dec (″)',fontsize=20)
ticks = [-7.5, -5, -2.5, 0, 2.5, 5, 7.5]
xxticks = [1.5, 1, 0.5, 0, -0.5]
yyticks = [2, 2.5, 3, 3.5, 4]
if apply_lim:
    plt.xticks(xxticks)
    plt.yticks(yyticks)
else:
    plt.xticks(ticks)
    plt.yticks(ticks)
plt.show()

## Apply ridgecrawler and fit_spline
rl = {}
rl_fitted = {}

for ridge_name in all_ridge_names:
    params = ridge_params[ridge_name]
    rl[ridge_name] = ridgecrawler(ridge_name, params[0], params[1], params[2], np.round(params[3]*step_factor), int(params[4]/step_factor), params[5], plot_sample_points, plot_slices, pixel)
    rl_fitted[ridge_name] = fit_spline(ridge_name, rl)
    
    

## Apply ridge_displacement
ridge_dist = {}
mean_dist = {}
sem_dist = {}
arc_len = []
rl_fitted_selected = {} # for plotting
scan_theta = {}

for pair_name in all_pairs:
    epoch_1_name = paired[pair_name][0]
    epoch_2_name = paired[pair_name][1]
    
    ridge_dist[pair_name], rl_fitted_selected[epoch_1_name], rl_fitted_selected[epoch_2_name], scan_theta[pair_name]  = ridge_displacement(epoch_1_name, epoch_2_name, rl_fitted, pixel, [c0,c1])
    mean_dist[pair_name] = np.mean(ridge_dist[pair_name])
    if use_sem == 1:
        sem_dist[pair_name] = np.std(ridge_dist[pair_name])/np.sqrt(len(ridge_dist[pair_name]))
    elif use_sem == 2:
        sem_dist[pair_name] = np.std(ridge_dist[pair_name])/np.sqrt(sem_factor)
    else:
        sem_dist[pair_name] = np.std(ridge_dist[pair_name])
    arc_len.append(len(ridge_dist[pair_name]))
    
    #plot_ridges(epoch_1_name, epoch_2_name, epoch_1_skeleton, epoch_2_skeleton, rl, rl_fitted, rl_fitted_selected)




## Calculate real displacement
print()

pixelscale = 45 * 1/3600000 / 180 * np.pi # radians per pixel
distance = distance_in_kpc * 3.0857e16 # km
time = (dates[compare[1]]-dates[compare[0]]).days * 24 * 3600 # seconds
speed = np.array(list(mean_dist.values())) * pixelscale * distance / time # km/s

time_year = (dates[compare[1]]-dates[compare[0]]).days/365.25 # years
speed_mas = np.array(list(mean_dist.values())) * 45 / time_year # mas/year
sem_mas = np.array(list(sem_dist.values())) * 45 / time_year

# Speed in km/s
print("Wind speed from",compare[0],'to',compare[1],'(km/s)')
for i in range(len(speed)):
    print(list(paired.keys())[i], " ", speed[i])

print("\nMean: " + str(np.sum(speed*np.array(arc_len))/np.sum(arc_len)))
print("SD: " + str(np.std(speed)/np.sqrt(len(all_pairs))))

# Speed in mas/year
print("\n\nWind speed from",compare[0],'to',compare[1],'(mas/year)')
for i in range(len(speed_mas)):
    print(list(paired.keys())[i], " ", speed_mas[i],'+-',sem_mas[i])

print("\nMean: " + str(np.sum(speed_mas*np.array(arc_len))/np.sum(arc_len)))
print("SEM: " + str(np.std(speed_mas)/np.sqrt(len(all_pairs))))

#plt.text(np.median(rl['a1'], axis=1)[0]*pixel-lim+0.2,np.median(rl['a1'], axis=1)[1]*pixel-lim+0.2,'hey',color='w')

if apply_lim:
    plt.xlim([1.7, -0.9])
    plt.ylim([1.8, 4.4])

if save_rc:
    plt.savefig('C://users//Yinuo//Desktop//High_rc'+figsuffix+'.png', dpi = 300)

## Plot difference
#ticks = [-8, -6, -4, -2, 0, 2, 4, 6, 8]
ticks = [-7.5, -5, -2.5, 0, 2.5, 5, 7.5]

plt.subplots(figsize=(8,7))
diff = epoch1_sharp_chopped/np.max(epoch1_sharp_chopped) - epoch2_sharp_chopped/np.max(epoch2_sharp_chopped)
plt.imshow(np.clip(img.gaussian_filter(diff,1),-0.15,0.15),origin='lower',extent=[lim,-lim,-lim,lim],cmap = colourmap)
#plt.title('Displacement from '+compare[0]+' to '+compare[1])
plt.xlabel('Relative RA (″)',fontsize=20)
plt.ylabel('Relative Dec (″)',fontsize=20)
plt.xticks(ticks)
plt.yticks(ticks)
plt.plot(0,0) # Somehow makes axes look right
plt.show()

if save_diff:
    plt.savefig('C://users//Yinuo//Desktop//High_diff'+figsuffix+'.png', dpi = 300)


## Plot angular dependance
'''
def sliding_mean(data_array, window=5):
    data_array = np.array(data_array)
    new_list = []
    for i in range(len(data_array)):
        indices = range(max(i - window + 1, 0),  
                        min(i + window + 1, len(data_array)))
        avg = 0
        for j in indices:
            avg += data_array[j]
        avg /= float(len(indices))
        new_list.append(avg)
    return np.array(new_list)


angles = np.concatenate([scan_theta[x] for x in scan_theta], 0)
distances = np.concatenate([ridge_dist[x] for x in ridge_dist], 0)

#x = []
#y = []
plt.figure(3)
for pair_name in paired:
    plt.plot(scan_theta[pair_name]/np.pi*180,sliding_mean(ridge_dist[pair_name],2),c='cornflowerblue')
    #x += list(scan_theta[pair_name]/np.pi*180)
    #y += list(sliding_mean(ridge_dist[pair_name],2))
plt.xlabel('Angle (degrees)')
plt.ylabel('Displacement (pixels)')

#plt.plot(angles, distances)

# Model
theta = 30
phi = 0
rmax = 5

rev = np.arange(phi,phi+2*np.pi,0.01)
dis = rmax * (np.cos(rev)**2 + np.sin(rev)**2*np.cos(theta)**2)**0.5
#plt.plot(x,y)
#plt.plot(mx, my)
plt.plot(rev/np.pi*180,dis,'--',c='steelblue')
#plt.legend(['a', 'b', 'c', 'd', 'e', 'f', 'planar disk model'])
plt.show()


plt.figure()
plt.imshow(-np.clip(epoch_2_skeleton,0,0.005),origin='lower', cmap='gray', extent=[lim,-lim,-lim,lim])
plt.xlabel('Relative RA (″)')
plt.ylabel('Relative Dec (″)')
plt.show()
'''
