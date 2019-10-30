'''
Converts from AB mag to Jy. 
'''

def abmag2jy(mag_AB):
    flux_Jy = 3631 * 10**(mag_AB/(-2.5))
    return flux_Jy