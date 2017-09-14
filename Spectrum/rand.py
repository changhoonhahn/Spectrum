'''

Random data for calculating power/bispectrum

'''
import numpy as np
import time 
# --- Local ---
from data import Data
from util import data_dir as direc
from ChangTools.fitstables import mrdfits


class Random(Data): 
    def __init__(self, cat_corr, **kwargs): 
        ''' 
        Child class of Data class to describe the random data of 
        simulation/data catalogs. Altogether this class does not 
        inherit much from the Data class aside from cosmology and
        such. However, child class is designated for consistency. 
        ''' 

        super(Random, self).__init__(cat_corr, **kwargs)

    def file(self): 
        ''' 
        Override parent class file() method in order to return
        file name of random catalog file. Everything is quite
        hardcoded.
        ''' 

        cat_name = ((self.cat_corr)['catalog'])['name'].lower()
        cat_dict = (self.cat_corr)['catalog']
        corr_name = ((self.cat_corr)['correction'])['name'].lower()
    
        data_dir = direc('data', self.cat_corr)

        if 'cmass' in cat_name: # CMASS

            if cat_name == 'cmass': 
                # CMASS random catalog 
                file_name = 'cmass-dr12v4-N-Reid.ran.dat'

            elif 'cmasslowz' in cat_name:  

                # CMASS LOWZ combined random catalog
                if 'e2' in cat_name: 
                    cmasslowz_str = 'e2'
                elif 'e3' in cat_name: 
                    cmasslowz_str = 'e3'
                else: 
                    cmasslowz_str = ''
                        
                if '_low' in cat_name: 
                    zbin_str = '_LOW' 
                elif 'high' in cat_name: 
                    zbin_str = '_HIGH'
                else: 
                    raise NameError("Must specify redshift bin of CMASS LOWZ sample") 

                file_name = ''.join([
                    'random0_DR12v5_CMASSLOWZ', 
                    cmasslowz_str.upper(), 
                    zbin_str, 
                    '_North.ran.dat'
                    ])
            else: 
                raise NotImplementedError()

        elif cat_name == 'nseries':     # Nseries

            file_name = 'Nseries_cutsky_randoms_50x_redshifts_comp.dat'

        elif cat_name == 'qpm':         # QPM 

            file_name = 'a0.6452_rand50x.dr12d_cmass_ngc.vetoed.dat'

        elif cat_name == 'tilingmock':  # Tiling Mock 

            file_name = 'randoms-boss5003-icoll012-vetoed.zlim.dat'

        elif cat_name == 'bigmd':       # Big MultiDark

            file_name = 'BigMD-cmass-dr12v4-RST-standHAM-Vpeak_vetoed.ran'

        elif cat_name == 'qso_bigmd': 
            
            if cat_dict['version'] == 'evo': 
                file_name = 'QSO-BigMD-1481.43deg2-bin0_v1.0-evo.ran'

            elif cat_dict['version'] == 'noevo': 
                file_name = 'QSO-BigMD-1481.43deg2-bin0_v1.0-noevo.ran'

            elif cat_dict['version'] == 'eboss': 
                file_name = 'eboss_v1.0-QSO-NS-eboss_v1.0-bin0.ran'

            elif cat_dict['version'] == 'ebossv1.5': 
                file_name = 'eboss_v1.5-QSO.ran'

            elif cat_dict['version'] == 'ebossnew': 
                file_name = 'rand20_y1_comp_cut_0.5_double_weight_col_Z.ran'

            elif 'jackknife' in cat_dict['version']: 
                n_jack = cat_dict['version'].split('jackknife')[-1]
                file_name = ''.join(['QSO-bin0-jackknife', str(n_jack), '.ran'])
                
            elif 'v2' in cat_dict['version']: 
                if 'z' in cat_dict['version']: 
                    file_name = 'BigMDPLv2-QSOZ.ran' 
                elif 'nsat' in cat_dict['version'] : 
                    file_name = 'BigMDPLv2-QSO-NSAT.ran'
                else: 
                    file_name = 'BigMDPLv2-QSO.ran'
        else: 
            raise NotImplementedError()

        if corr_name == 'noweight': 
            file_name = ''.join(['.'.join(file_name.rsplit('.')[:-1]), '.noweight.ran.dat'])
        
        return  ''.join([data_dir, file_name])



if __name__=="__main__":
    for cat in ['bigmd', 'nseries', 'qpm']: 
        cat_corr = {'catalog': {'name': cat}, 'correction': {'name': 'noweight'}}
        NoweightRandom(cat_corr)
        #rand_class = Random(cat_corr)
        #print rand_class.file()
        #rand_class.build()
