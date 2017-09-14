'''

code for convenience 

'''
import os
import sys
import numpy as np 
import scipy as sp
import cosmolopy as cosmos


def t_mod(file):   
    ''' return the modified time of `file`. if `file` does 
    not exist, returns 0 
    '''
    if not os.path.isfile(file): 
        return 0 
    else: 
        return os.path.getmtime(file)


def data_dir(type, catname):
    ''' Return directory location of data/fft/spec type
    data given catalog name . 
    '''
    if type not in ('data', 'fft', 'spec'): 
        raise ValueError()
    # needs to be updated to make more sense  
    if catname == 'lasdamasgeo': 
        catdir = 'LasDamas/Geo/'
    elif catname == 'tilingmock': 
        catdir = 'tiling_mocks/'
    elif catname == 'qpm': 
        catdir = 'qpm/dr12d/'
    elif catname == 'nseries': 
        catdir = 'nseries/'
    elif catname == 'patchy': 
        catdir = 'PATCHY/dr12/v6c/'
    elif 'bigmd' in catname:
        catdir = 'BigMD/'
    elif 'cmass' in catname: 
        if catname == 'cmass': 
            catdir = 'CMASS/'
        elif 'cmasslowz' in catname: 
            catdir = 'CMASS/dr12v5/'
    else: 
        raise NotImplementedError()
    return  ''.join([dat_dir(), type, '/', catdir])


def code_dir(): 
    ''' Directory where all the code is located (the directory that this file is in!)
    '''
    return os.path.dirname(os.path.realpath(__file__))+'/'


def fortran_dir(): 
    return code_dir()+'fortran/'


def dat_dir(): 
    ''' dat directory is symlinked to a local path where the data files are located
    '''
    return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))+'/dat/'


def fig_dir(): 
    ''' 
    '''
    return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))+'/figs/'


def comdis2z(comdis, **cosmo): 
    ''' Given comoving distance and cosmology, determine z 
    using cubic spline

    Notes
    -----
    * Comoving distance *has* to be in Mpc/h
    '''
    z_arr = np.array([0.0+0.05*np.float(i) for i in range(21)]) 
    dm_arr = cosmos.distance.comoving_distance(z_arr, **cosmo)*cosmo['h']

    dmz_spline = sp.interpolate.interp1d(dm_arr, z_arr, kind='cubic') 
    
    z = dmz_spline(comdis)

    return z 


def find_nearest(array, value, index=False):
    ''' Find nearest element in array to value. If index is True then return index
    '''
    idx = (np.abs(array-value)).argmin()
    if index == False: 
        return array[idx]
    else: 
        return idx


def radecz_to_xyz(ra, dec, z, **cosmo):
    ''' Given RA, Dec, redshift AND cosmology, calculate x,y,z in Mpc/h
    '''
    phi = ra 
    theta = 90.0 - dec 
    r = cosmos.distance.comoving_distance(z, **cosmo)*cosmo['h']    # Mpc/h

    x = r * np.cos(np.deg2rad(phi)) * np.sin(np.deg2rad(theta)) 
    y = r * np.sin(np.deg2rad(phi)) * np.sin(np.deg2rad(theta))
    z = r * np.cos(np.deg2rad(theta))
    return (x,y,z)


def ang_sep(ra1, dec1, ra2, dec2): 
    ''' Given a pair of ra and decs in DEGREES gives angular separation in DEGREES
    '''
    # convert to radians 
    ra1 = ra1*np.pi/180.
    dec1 = dec1*np.pi/180.
    ra2 = ra2*np.pi/180.
    dec2 = dec2*np.pi/180.

    x = np.cos(ra1)*np.cos(dec1)*np.cos(ra2)*np.cos(dec2) 
    y = np.sin(ra1)*np.cos(dec1)*np.sin(ra2)*np.cos(dec2) 
    z = np.sin(dec1)*np.sin(dec2)

    rad = np.arccos(x+y+z)
    
    sep = rad
    #sep = np.choose( rad<0.000004848 , ( np.sqrt( (np.cos(dec1)*(ra1-ra2))**2+(dec1-dec2)**2), rad))

    sep = sep*180./np.pi
    return sep
