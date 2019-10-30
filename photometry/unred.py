'''
Performs interstellar dereddening based on Mathis 1990.
'''


def unred(fluxes, wavelengths, A_wav, wav, out_cloud=0):
    import numpy as np
    import scipy.interpolate as intp
    
    lam=[2, 4, 23, 41, 73, 91, 120, 130, 150, 180, 200, 218, 240, 260,  280, 330, 365, 440, 550, 700, 900, 1250, 1650, 2200, 3400, 5000, 7000,  9000, 9700, 10000, 12000, 15000, 18000, 20000, 6e4, 1e5, 2e5]
    lam=np.array(lam)/1000.0 # now in microns
    red=[1.35, 3.39, 7.31, 9.15, 19.1, 17.2, 12.71, 11.09, 9.44, 8.93, 10.08, 11.29, 9.03, 7.63, 6.90, 5.87, 5.53, 4.70, 3.55, 2.66, 1.70, 1.00, 0.624, 0.382, 0.182, 0.095, 0.070, 0.157, 0.208, 0.192, 0.098, 0.053, 0.083, 0.075, 0, 0, 0]
    
    if out_cloud:
        lam=[120, 130, 150, 180, 200, 218, 240, 260, 280, 330, 365, 440, 550, 700, 900, 1250, 1650, 2200, 3400, 5000, 7000, 9000, 9700, 10000,  12000, 15000, 18000, 20000]
        lam=np.array(lam)/1000.0 # now in microns
        red=[5.32, 4.89, 4.57, 4.66, 5.32, 6.03, 5.13, 4.59, 4.34, 4.12, 4.07, 3.67, 3.06, 2.43, 1.70, 1.00, 0.624, 0.382, 0.182, 0.095, 0.070, 0.157, 0.208, 0.192, 0.098, 0.053, 0.083, 0.075]
    
    if(max(wavelengths) > max(lam) or min(wavelengths) < min(lam)):
        print('### ERROR - lam out of interpolation range ###')
    
    tck = intp.splrep(lam,red)
    
    red_wav0 = intp.splev(wav,tck)
    reddening = intp.splev(wavelengths,tck) * A_wav / red_wav0
    
    star_flx = fluxes/10**(reddening/(-2.5))
    
    return star_flx