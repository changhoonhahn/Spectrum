'''

Interface with FORTRAN code for galaxy powerspectrum and bispectrum calculations 

'''
import subprocess
import os.path

import util as UT 

class Fcode: 

    def __init__(self, type, catinfo): 
        ''' Class for all the FORTRAN stuff for P(k)/B(k) calculations 
        '''
        self.catinfo = catinfo.copy()
        self.type = type

        specdict = catinfo['spec']
        
        self.code = self.code_name()
        self.exe = self.exe_name()  

    def code_name(self): 
        ''' code name 
        '''
        if self.type == 'fft': # fft code
            multipole_str = ''
            if specdict['ell'] != 0: 
                multipole_str = '_quad'
            corr_str = ''
            if 'correction' in self.catinfo.keys():
                if self.catinfo['correction']['name'] == 'noweight': 
                    corr_str = '_noweight'
                elif self.catinfo['correction']['name'] == 'floriansn': 
                    corr_str = '_florian'
                elif self.catinfo['correction']['name'] == 'hectorsn': 
                    corr_str = '_hector'
            f_name = ''.join(['FFT_fkp', multipole_str, corr_str, '.f'])
        elif type == 'pk': # P(k) code
            if specdict['ell'] == 0: 
                f_name = 'power-Igal-Irand.f' 
            else:
                f_name = 'power_quad.f'
        elif type == 'bk': # B(k1,k2,k3) code
            f_name = 'bisp_fast_bin_fftw2_quad.f' 
        else: 
            raise NotImplementedError
        
        f_code = ''.join([UT.fortran_dir(), f_name])
        return f_code

    def exe_name(self): 
        ''' name of compiled fortran executable.
        '''
        code_file = (self.code).split('/')[-1]
        return ''.join([UT.fortran_dir(), 'exe/', '.'.join(code_file.rsplit('.')[:-1]), '.exe'])

    def compile(self):
        ''' Compile fortran code
        '''
        specdict = self.catinfo['spec']

        # compile command for fortran code. Quadruple codes have more
        # complex compile commands specified by Roman 
        if (self.type == 'fft') and (specdict['ell'] != 0): 
            compile_cmd = ' '.join([
                'ifort -fast -o', 
                self.exe, 
                self.code, 
                '-L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw -lm'
                ])
        elif (self.type == 'pk') and (specdict['ell'] != 0): 
            compile_cmd = ' '.join([
                'ifort -fast -o', 
                self.exe, 
                self.code
                ])
        elif (self.type == 'bk'): 
            compile_cmd = ' '.join([
                'ifort -fast -o', 
                self.exe, 
                self.code, 
                '-L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw'
                ])
        else: 
            compile_cmd = ' '.join([
                'ifort -O3 -o', 
                self.exe, 
                self.code, 
                '-L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw'
                ])

        print ' ' 
        print 'Compiling -----'
        print compile_cmd
        print '----------------'
        print ' ' 

        # call compile command 
        subprocess.call(compile_cmd.split())
        return None 

    def commandline_call(self, **kwargs): 
        ''' Command line call for Fortran code
        '''
        specdict = self.catinfo['spec']
        
        fcode_t_mod, fexe_t_mod = self.mod_time()
        if fcode_t_mod < fcode_t_mod: 
            raise ValueError("Compile failed")

        fort_exe = self.fexe() 

        if self.type == 'fft': 

            if specdict['ell'] == 0:    # monopole

                if not all([kwarg in kwargs.keys() for kwarg in ['DorR', 'datafile', 'fftfile']]):

                    err_msg = ''.join([
                        'For type = ', self.type, 
                        " 'DorR', 'datafile', 'fftfile' must be specified in kwargs"
                        ])
                    raise KeyError(err_msg)
                
                catname = ((self.catinfo)['catalog'])['name']
                specdict = (self.catinfo)['spec']

                if catname == 'nseries': 
                    n_cat = 7
                elif catname == 'qpm': 
                    n_cat = 3 
                elif catname == 'bigmd':
                    n_cat = 9 
                elif catname == 'cmass': 
                    n_cat = 1 
                elif catname == 'qso_bigmd':
                    n_cat = 10
                else: 
                    raise NotImplementedError()
                
                try: 
                    if kwargs['cosmology'] == 'survey': 
                        n_cosmo = 1
                    else: 
                        n_cosmo = 0 
                except KeyError: 
                    n_cosmo = 0 

                if kwargs['DorR'] == 'data': 
                    n_DorR = 0 
                elif kwargs['DorR'] == 'random': 
                    n_DorR = 1 

                cmdline_call = ' '.join([
                    self.exe, 
                    str(n_cat), 
                    str(n_cosmo), 
                    str(specdict['Lbox']), 
                    str(specdict['Ngrid']), 
                    str(n_DorR),
                    str(specdict['P0']), 
                    kwargs['datafile'], 
                    kwargs['fftfile']
                    ])
            else: 
                if not all([kwarg in kwargs.keys() for kwarg in ['DorR', 'datafile', 'fftfile']]):

                    err_msg = ''.join([
                        'For type = ', self.type, 
                        " 'DorR', 'datafile', 'fftfile' must be specified in kwargs"
                        ])
                    raise KeyError(err_msg)
                
                catname = ((self.catinfo)['catalog'])['name']
                specdict = (self.catinfo)['spec']

                if catname == 'nseries': 
                    n_cat = 7
                elif catname == 'qpm': 
                    n_cat = 3
                elif catname == 'bigmd': 
                    n_cat = 9 
                elif catname == 'cmass': 
                    n_cat = 1
                elif catname == 'qso_bigmd':
                    n_cat = 10
                else: 
                    raise NotImplementedError()
                
                try: 
                    if kwargs['cosmology'] == 'survey': 
                        n_cosmo = 1
                    else: 
                        n_cosmo = 0 
                except KeyError: 
                    n_cosmo = 0 

                if kwargs['DorR'] == 'data': 
                    n_DorR = 0 
                elif kwargs['DorR'] == 'random': 
                    n_DorR = 1 

                n_interp = 4    # 4th order interlaced CIC

                cmdline_call = ' '.join([
                    self.exe, 
                    str(n_cat), 
                    str(n_cosmo), 
                    str(specdict['Lbox']), 
                    str(specdict['Ngrid']), 
                    str(n_interp), 
                    str(n_DorR),
                    str(specdict['P0']), 
                    kwargs['datafile'], 
                    kwargs['fftfile']
                    ])

        elif self.type == 'pk': 

            if specdict['ell'] == 0:    # monopole

                if not all([kwarg in kwargs.keys() for kwarg in ['datafft', 'randfft', 'powerfile']]):
                    err_msg = ''.join([
                        'For type = ', self.type, 
                        " must be specified in kwargs"
                        ])
                    raise KeyError(err_msg)
                
                catname = ((self.catinfo)['catalog'])['name']
                specdict = (self.catinfo)['spec']

                nbins = int(specdict['Ngrid']/2) # number of bins

                cmdline_call = ' '.join([
                    self.exe, 
                    kwargs['datafft'], 
                    kwargs['randfft'],
                    kwargs['powerfile'], 
                    str(specdict['Lbox']), 
                    str(nbins),
                    ])

            else:   # quadrupole, hexadecapole, etc

                if not all([kwarg in kwargs.keys() for kwarg in ['datafft', 'randfft', 'powerfile']]):
                    err_msg = ''.join([
                        "For type = ", self.type, "datafft, randfft, and powerfile must be specified in kwargs"
                        ])
                    raise KeyError(err_msg)
                
                specdict = (self.catinfo)['spec']

                nbins = int(specdict['Ngrid']/2) # number of bins hardcoded to be Ngrid/2

                cmdline_call = ' '.join([
                    self.exe, 
                    kwargs['datafft'], 
                    kwargs['randfft'],
                    kwargs['powerfile'], 
                    str(specdict['Lbox']), 
                    str(nbins),
                    ])
        elif self.type == 'bk':
            
            if specdict['Ngrid'] != 360: 
                raise ValueError('Ngrid has to be equal to 360 due to count file')

            count_file = '/home/users/rs123/Code/Fortran/counts2quad_n360_nmax40_ncut3_s3'            
            if not os.path.isfile(count_file): 
                raise NotImplementedError('Count File does not exist') 
            
            if not all([kwarg in kwargs.keys() for kwarg in ['datafft', 'randfft', 'powerfile']]):
                err_msg = ''.join([
                    'For type = ', self.type, 
                    " must be specified in kwargs"
                    ])
                raise KeyError(err_msg)
            
            # bispectrum code input: period/data, random fft file, data fft file, bispectrum file  
            cmdline_call = ' '.join([
                self.exe, 
                '2', 
                kwargs['randfft'], 
                kwargs['datafft'], 
                kwargs['powerfile']
                ]) 
        else: 
            raise NotImplementError()

        return cmdline_call
