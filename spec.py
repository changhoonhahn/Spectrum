'''

Spectrum class of FiberCollisions project

'''

import numpy as np
import os.path
import subprocess

from rand import Random
from corrdata import CorrData 

from fft import Fft
from util.direc import direc
from util.fortran import Fcode

# Classes ------------------------------------------------------------
class Spec(object): 

    def __init__(self, spectype, cat_corr, ell=None, **kwargs):
        """ Class that describes power/bispectrum measurements 
        
        specify catalog, version, mock file number, file specifications (e.g. Nrandom), 
        fiber collision correction method, correction specifications (e.g. sigma, fpeak)

        Parameters 
        ----------
        spectype : 'power' or 'bispec'
        cat_corr : catalog and correction dictionary 


        """

        if spectype not in ['pk', 'bk']: 
            raise ValueError()
        else: 
            self.type = spectype

        if 'spec' not in cat_corr.keys(): 
            if self.type == 'bk': 
                ell = 2         # this is a hack so that the FFT is from the quadrupole FFT

            if ell is None: 
                raise ValueError("Specify ell (monopole: 0, quadrupole: 2, hexadecapole: 4)")

            # default spectrum parameters
            cat_corr['spec'] = {
                    'P0': 20000, #P0 
                    'Lbox': 3600, 
                    'Ngrid':360, 
                    'ell': ell 
                    }
            self.ell = ell

            if 'Ngrid' in kwargs.keys(): 
                cat_corr['spec']['Ngrid'] = kwargs.pop('Ngrid')
        else: 

            if self.type == 'bk': 
                cat_corr['spec']['ell'] = 2
            elif 'ell' not in cat_corr['spec'].keys():
                raise ValueError("Specify ell (monopole: 0, quadrupole: 2, hexadecapole: 4) in catcorr dictionary")

            if ell is not None: 
                if ell != cat_corr['spec']['ell']: 
                    raise ValueError

            self.ell = cat_corr['spec']['ell']
        
        self.cat_corr = cat_corr.copy()
        self.kwargs = kwargs

        self.file_name = self.file()
            
        if self.type == 'bk':
            self.scale = np.float(self.cat_corr['spec']['Lbox'])
            k_fund = (2.0*np.pi)/self.scale        # k fundamental 
            self.k_fund = k_fund 
    
    def read(self): 
        """ Read power/bispectrum of simulated/observed data catalog
        """
        if self.type == 'pk':   # power spectrum
            if self.ell == 0:   # monopole
                    
                col_index = [0, 1]
                data_cols = ['k', 'p0k']

            elif self.ell == 2:     # quadrupoel
                
                col_index = [0, 2, 1, 3]
                data_cols = ['k', 'p2k', 'p0k', 'p4k']

            elif self.ell == 4:     # hexadecapole
                
                col_index = [0, 3, 1, 2]
                data_cols = ['k', 'p4k', 'p0k', 'p2k']

            else: 
                raise NotImplementedError()
        else:                   # bispectrum
            # bispectrum v5 columns k1, k2, k3, P0(k1), P0(k2), P0(k3), B0, Q0, P2(k1), P2(k2), P2(k3), B2, Q2, dum, dum 
            col_index = [0, 1, 2, 3, 4, 5, 6, 7]
            data_cols = ['k1', 'k2', 'k3', 'p0k1', 'p0k2', 'p0k3', 'bk', 'qk']

        spec_data = np.loadtxt(
                    self.file_name, 
                    unpack = True, 
                    usecols = col_index
                    )

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
        """ power/bispectrum file 
        """

        specdict = (self.cat_corr)['spec']

        # powerspectrum or bispectrum 
        if self.type in ('pk'): 
            spec_str = 'POWER_'
        elif self.type == 'bk':
            spec_str = 'BISP_'

        #if 'quad' not in specdict.keys(): 
        #    specdict['quad'] = False
        
        #if specdict['quad']:          
        #    spec_str += 'Q_'
        if (self.type == 'pk') and (self.ell != 0): 
            spec_str += 'Q_'

        gal_data = CorrData(self.cat_corr, **self.kwargs)
        self.data_file = gal_data.file_name
        gal_file = (gal_data.file_name).split('/')[-1]

        rand_data = Random(self.cat_corr, **self.kwargs)
        self.random_file = rand_data.file_name

        spec_dir = direc('spec', self.cat_corr)

        if self.type == 'pk': 
            specparam_str = ''.join([
                '.grid', str(specdict['Ngrid']), 
                '.P0', str(specdict['P0']), 
                '.box', str(specdict['Lbox'])
                ])
        elif self.type == 'bk': 
            specparam_str = ''.join([
                '.grid', str(specdict['Ngrid']), 
                '.nmax40.ncut3.s3', 
                '.P0', str(specdict['P0']), 
                '.box', str(specdict['Lbox'])
                ])
        else: 
            raise NotImplementedError
    
        file_name = ''.join([
            spec_dir, 
            spec_str,
            gal_file, 
            specparam_str
            ])

        return file_name

    def build(self): 
        """ Calculate power/bispectrum of simulated/observed data catalog 
        """
        
        if 'clobber' not in (self.kwargs).keys(): 
            bool_clobber = False
        else: 
            bool_clobber = self.kwargs['clobber']

        #if corrdict['name'] == 'fourier_tophat':
        #    if self.ell != 2: 
        #        raise ValueError

        #    true_cat_corr = {
        #            'catalog': catdict, 
        #            'correction': {'name': 'true'}
        #            }
        #    tr_gal = Data('data', true_cat_corr)

        #    fourier_tophat_Pk(self.cat_corr, self.file_name, tr_gal.file_name)
        #    return None

        codeclass = Fcode(self.type, self.cat_corr) 
        spec_code = codeclass.code
        spec_exe = codeclass.fexe()
        
        # code and exe modification time 
        speccode_t_mod, specexe_t_mod = codeclass.mod_time()

        if specexe_t_mod <= speccode_t_mod: 
            codeclass.compile()

        # fft files 
        datafft = Fft('data', self.cat_corr, **self.kwargs)
        if not os.path.isfile(datafft.file_name+'_0') or bool_clobber:
            datafft.build()

        randfft = Fft('random', self.cat_corr, **self.kwargs)
        if not os.path.isfile(randfft.file_name+'_0'): 
            randfft.build()
        
        spec_cmd = codeclass.commandline_call(
                datafft = datafft.file_name, 
                randfft = randfft.file_name, 
                powerfile = self.file_name
                )
        print spec_cmd

        if any([not os.path.isfile(self.file_name), bool_clobber]):
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
    cat_corr = {'catalog': {'name': 'nseries', 'n_mock': 1}}
    spectrum = Spec('bk', cat_corr, Ngrid=360)
    print spectrum.file()
    print spectrum.build()
