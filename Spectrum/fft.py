"""

FFT of galaxy catalogs. 

"""
import numpy as np 
import os.path
import time 
import subprocess
# -- local -- 
import util as UT 
from data import Data
from fortran import Fcode


class Fft(object): 
    def __init__(self, catinfo): 
        ''' class object for the FFT of a given catinfo dictionary 
        '''
        self.catinfo = catinfo.copy() 
        if 'spec' not in self.catinfo.keys(): 
            # default pk/bk parameters
            self.catinfo['spec'] = {
                    'P0': 20000, #P0 
                    'Lbox': 3600, 
                    'Ngrid':360, 
                    'ell': 0}

        self.data_obj = Data(catinfo)
        self.data_file = self.data_obj.file_name # data file name 
        self.file_name = self.file()

    def file(self): 
        ''' FFT data file name 
        '''
        # spectrum specifiers 
        specdict = self.catinfo['spec'] 
        str_spec = ''.join(['.grid', str(specdict['Ngrid']), '.P0', str(specdict['P0']), '.box', str(int(specdict['Lbox']))])

        # FFT label 
        str_fft = 'FFT_'
        if specdict['ell'] != 0: 
            str_fft += 'Q_'

        # FFTs from data file 
        fft_file = ''.join([UT.data_dir('fft', self.data_obj.catalog), 
            str_fft, (self.data_file).rsplit('/')[-1], str_spec])
        return fft_file  

    def build(self, clobber=False): 
        ''' Run FFT FORTRAN code to calculate FFT of data
        '''
        specdict = self.catinfo['spec'] 

        if not os.path.isfile(self.data_file):
            # if the data doesn't exist, build data first        
            print('data file does not exist!')
            self.data_obj.build()

        fcode = Fcode('fft', self.catinfo) 
        fftcode, fftexe = fcode.code, fcode.exe
        
        # code and exe modification time 
        fftcode_t_mod, fftexe_t_mod = codeclass.mod_time()
        if fftexe_t_mod < fftcode_t_mod: 
            print('FFT code was modified after .exe file was compiled') 
            #codeclass.compile()

        # command line call 
        FFTcmd = fcode.commandline_call(
                    DorR = self.type, 
                    datafile = self.data_file,
                    fftfile = self.file_name
                    ) 

        if any([not os.path.isfile(self.file_name), clobber]):
            print ''
            print '-----------------------'
            print 'Constructing '
            print self.file_name  
            print '-----------------------'
            print ''
            print FFTcmd
            print '-----------------------'
            subprocess.call(FFTcmd.split())
        else: 
            print ''
            print '-----------------------'
            print self.file_name  
            print 'Already Exists'
            print '-----------------------'
            print ''
        return None 


if __name__=='__main__': 
    cat_corr = {'catalog': {'name': 'qpm', 'n_mock': 1}, 'correction': {'name': 'true'}}
    for DorR in ['data', 'random']:
        fftee = Fft(DorR, cat_corr, clobber=True)
        print fftee.build()

    """
        # Quad FFT argument sequence (SUBJECT TO CHANGE) 
            # determine "idata"
            if catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz'): 
                idata = 2
                ifc = 0 
            elif catalog['name'].lower() == 'qpm': 
                idata = 3 
                ifc = 0 
            elif catalog['name'].lower() == 'tilingmock':
                idata = 9 
                ifc = 0 
            elif catalog['name'].lower() == 'nseries': 
                idata = 10 
                ifc = 0 
            elif catalog['name'].lower() == 'ldgdownnz':  
                idata = 11 
                ifc = 0 
            elif 'bigmd' in catalog['name'].lower(): 
                idata = 12
                ifc = 0 
            elif catalog['name'].lower() == 'cmass': 
                idata = 13 
                ifc = 0 
            else: 
                raise NameError('not included in Quadrupole code') 
                
            # bash commend 
            # is of the form 
            # FFT_FKP_BOSS_cic_il4_v3.exe idata box Ngrid interpol iflag P0  ifc icomp input_file output_file
            # icomp is hardcoded 0 so that it takes into account completeness!
            FFT_cmd = ' '.join([
                FFT_exe, str(idata), 
                str(spec['box']), str(spec['grid']), 
                "4", str(DorR_number), str(spec['P0']), 
                str(ifc), "0", data_file, fft_file]) 
            print FFT_cmd

            if DorR.lower() == 'data':  # don't bother checking if the file exists for mocks and run the damn thing 
                subprocess.call(FFT_cmd.split()) 

            elif DorR.lower() == 'random':      # random takes longer so check to see if it exists first
                # call FFT randomc ommand 
                if os.path.isfile(fft_file) == False: 
                    print "Building ", fft_file 
                    subprocess.call(FFT_cmd.split())
                else: 
                    print fft_file, " already exists" 

            print 'Constructing ', 
    """
