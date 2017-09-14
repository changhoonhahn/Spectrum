'''

Spectrum class of FiberCollisions project

'''
import numpy as np
import os.path
import subprocess

import util as UT
from corrdata import CorrData 

from fft import Fft
from fort import Fcode


class Spec(object): 
    def __init__(self, catinfo, **kwargs):
        ''' Class describing power/bispectrum measurements specify 
        catalog, version, mock file number, file specifications (e.g. Nrandom), 
        fiber collision correction method, correction specifications (e.g. sigma, fpeak)
        '''
        self.catinfo = catinfo.copy() 

        # make sure all the spectrum parameters are appropriately specified 
        if 'spec' not in self.catinfo.keys(): 
            raise ValueError
        if 'P0' not in self.catinfo['spec'].keys():
            self.catinfo['spec']['P0'] = 20000
        if 'Lbox' not in self.catinfo['spec'].keys():
            self.catinfo['spec']['Lbox'] = 3600
        if 'Ngrid' not in self.catinfo['spec'].keys(): 
            raise ValueError
        if 'ell' not in self.catinfo['spec'].keys(): 
            raise ValueError
        self.ell = self.catinfo['spec']['ell']
        if 'type' not in self.catinfo['spec'].keys(): 
            raise ValueError
        if self.catinfo['spec']['type'] not in ['pk', 'bk']: 
            raise ValueError
        self.type = self.catinfo['spec']['type']
        
        if self.type == 'bk':
            self.scale = np.float(self.catinfo['spec']['Lbox'])
            k_fund = (2.0*np.pi)/self.scale     # k fundamental 
            self.k_fund = k_fund 
       
        data_catinfo = catinfo.copy() 
        data_catinfo['catalog']['type'] = 'data'
        self.fft_data = Fft(data_catinfo)

        rand_catinfo = catinfo.copy() 
        rand_catinfo['catalog']['type'] = 'random'
        self.fft_rand = Fft(rand_catinfo)
        
        self.file_name = self.file()

    def read(self): 
        ''' Read power/bispectrum of simulated/observed data catalog
        '''
        if self.type == 'pk':   # power spectrum
            if self.ell == 0:   # monopole
                col_index = [0, 1, -1]
                data_cols = ['k', 'p0k', 'count']
            elif self.ell == 2:     # quadrupoel
                col_index = [0, 2, 1, 3, -2]
                data_cols = ['k', 'p2k', 'p0k', 'p4k', 'count']
            elif self.ell == 4:     # hexadecapole
                col_index = [0, 3, 1, 2, -2]
                data_cols = ['k', 'p4k', 'p0k', 'p2k', 'count']
            else: 
                raise NotImplementedError()
        else:                   # bispectrum
            # bispectrum v5 columns k1, k2, k3, P0(k1), P0(k2), P0(k3), B0, Q0, P2(k1), P2(k2), P2(k3), B2, Q2, dum, dum 
            col_index = [0, 1, 2, 3, 4, 5, 6, 7, 11, 12]
            data_cols = ['k1', 'k2', 'k3', 'p0k1', 'p0k2', 'p0k3', 'bk', 'qk', 'b2k', 'q2k']

        spec_data = np.loadtxt(self.file_name, unpack = True, usecols = col_index)
        for i_col, col in enumerate(data_cols): 
            setattr(self, col, spec_data[i_col])

        if self.type == 'bk': 
            self.k1 *= self.k_fund          # k1 * k_fundamental to get h/Mpc
            self.k2 *= self.k_fund 
            self.k3 *= self.k_fund 
            
            # some extra useful values
            self.i_triangle = range(len(self.k1))               # triangle ID
            self.avgk = (self.k1 + self.k2 + self.k3)/3.0       # average k
            self.kmax = np.amax(np.vstack([self.k1, self.k2, self.k3]), axis=0) # max(k1,k2,k3)
        return None 

    def file(self):
        ''' power/bispectrum file name
        '''
        # powerspectrum or bispectrum 
        if self.type  == 'pk': 
            spec_str = 'POWER_'
        elif self.type == 'bk':
            spec_str = 'BISP_'
        
        # string specific the catalog
        str_cat = self.fft_data.file_name.split('/')[-1].replace('FFT_', spec_str)

        if self.type == 'pk': 
            str_extra = ''
        elif self.type == 'bk': 
            str_extra = '.nmax40.ncut3.s3', 
    
        return ''.join([UT.data_dir('spec', self.catinfo['catalog']['name']), str_cat, str_extra])

    def build(self, clobber=False): 
        ''' Calculate power/bispectrum of simulated/observed data catalog 
        '''
        fcode = Fcode(self.type, self.catinfo) 
        spec_code, spec_exe = fcode.code, fcode.exe
        
        # code and exe modification time 
        t_code, t_exe = UT.t_mod(spec_code), UT.t_mod(spec_exe)  

        if t_exe <= t_code: 
            print('FFT code was modified after .exe file was compiled') 
            #codeclass.compile()

        # fft files 
        if not os.path.isfile(self.fft_data.file_name) or clobber:
            self.fft_data.build()

        if not os.path.isfile(self.fft_rand.file_name): 
            self.fft_rand.build()
        
        spec_cmd = fcode.commandline_call(
                datafft = self.fft_data.file_name, 
                randfft = self.fft_rand.file_name, 
                powerfile = self.file_name)
        print spec_cmd

        if any([not os.path.isfile(self.file_name), clobber]):
            print ''
            print '-----------------------'
            print 'Constructing '
            print self.file_name  
            print '-----------------------'
            print ''
            print spec_cmd
            print '-----------------------'

            subprocess.call(spec_cmd.split())
        else: 
            print ''
            print '-----------------------'
            print self.file_name  
            print 'Already Exists'
            print '-----------------------'
            print ''
        return None


if __name__=='__main__':
    cat_corr = {'catalog': {'name': 'qpm', 'n_mock': 1}}
    spectrum = Spec('pk', cat_corr, ell=0, Ngrid=360)
    #spectrum = Spec('bk', cat_corr, Ngrid=360)
    print spectrum.file()
    print spectrum.build()
