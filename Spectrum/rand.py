'''

Random data for calculating power/bispectrum

'''
import numpy as np
import time 
# --- Local ---
from defutility.fitstables import mrdfits
from data import Data
from util.direc import direc

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

        else: 
            raise NotImplementedError()
        
        return  ''.join([data_dir, file_name])

    def build(self): 
        ''' 
        Completely replace parent class build attribute in order
        to ruild the random catalogs from original survey data. 
        Everything is hardcoded. 
        '''
        cat_name = ((self.cat_corr)['catalog'])['name'].lower()

        if cat_name == 'nseries':     # Nseries ----------------------------

            # original random catalog 
            data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
            orig_rand_file = ''.join([data_dir, 'Nseries_cutsky_randoms_50x_redshifts.dat']) 
            ra, dec, z = np.loadtxt(orig_rand_file, unpack=True, usecols=[0,1,2]) # RA, Decl, Redhsift
        
            # sector completeness catalog
            orig_comp_file = ''.join([data_dir, 'Nseries_cutsky_randoms_50x_maskinfo.dat'])  
            comp = np.loadtxt(orig_comp_file, unpack=True, usecols=[0])
            
            header_str = 'Columns : ra, dec, z, comp'
            data_list = [ra, dec, z, comp]
            data_fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f']

        elif cat_name == 'qpm':             # QPM -----------------------------------

            data_dir = '/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/randoms/'
            ra, dec, z = np.loadtxt(
                    data_dir+'a0.6452_rand50x.dr12d_cmass_ngc.rdz', unpack=True, usecols=[0,1,2])   # ra, dec, z, wfkp
            comp = np.loadtxt(
                    data_dir+'a0.6452_rand50x.dr12d_cmass_ngc.rdz.info', unpack=True, skiprows=3, usecols=[1])   # galid, comp?
            veto = np.loadtxt(
                    data_dir+'a0.6452_rand50x.dr12d_cmass_ngc.veto')       # veto  

            vetomask = np.where(veto == 0)

            header_str = 'Columns : ra, dec, z, comp'
            data_list = [ra[vetomask], dec[vetomask], z[vetomask], comp[vetomask]]
            data_fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'] 

        elif 'cmass' in cat_name:          # CMASS -------------------------------- 

            data_dir = direc('data', self.cat_corr) 

            if cat_name == 'cmass': 
                # random data fits file
                data_file = ''.join([data_dir, 'cmass-dr12v4-N-Reid.ran.fits']) 
                cmass = mrdfits(data_file) 
            
                # mask file 
                mask_file = ''.join([data_dir, 'mask-cmass-dr12v4-N-Reid.fits']) 
                mask = mrdfits(mask_file) 
                ipoly = cmass.ipoly # polygon index
                comp = mask.weight[ipoly]
            
                # redshift limit 
                zlimit = np.where((cmass.z >= 0.43) & (cmass.z <= 0.7))

            elif 'cmasslowz' in cat_name:   
                # CMASS LOWZ combined data
                
                # three different CMASS LOWZ  
                if 'e2' in cat_name: 
                    cmasslowz_str = 'E2' 
                elif 'e3' in cat_name: 
                    cmasslowz_str = 'E3'
                else: 
                    cmasslowz_str = ''

                if 'high' in cat_name: 
                    zmin, zmax = 0.5, 0.75
                elif '_low' in cat_name:
                    zmin, zmax = 0.2, 0.5
                else: 
                    raise NameError("CMASSLOWZ Catalog must specify high or lowr edshift bin") 
                
                # mask file 
                mask_file = ''.join([
                    data_dir, 
                    'mask_DR12v5_CMASSLOWZ', 
                    cmasslowz_str, 
                    '_North.fits.gz'
                    ])
                start_time = time.time()
                print 'Reading ', mask_file 
                mask = mrdfits(mask_file) 
                print 'took ', (time.time() - start_time)/60.0, ' minutes'
                
                # random data fits file
                data_file = ''.join([
                    data_dir, 
                    'random0_DR12v5_CMASSLOWZ', 
                    cmasslowz_str, 
                    '_North.fits.gz'
                    ])
                start_time = time.time()
                print 'Reading ', data_file 
                cmass = mrdfits(data_file) 
                print 'took ', (time.time() - start_time)/60.0, ' minutes'
            
                ipoly = cmass.ipoly # polygon index
                comp = mask.weight[ipoly]
            
                # redshift limit 
                zlimit = np.where((cmass.z >= zmin) & (cmass.z < zmax))

            else: 
                raise NotImplementedError("Only CMASS and CMASS+LOWZ combined sample implemented") 
        
            header_str = 'columns : ra, dec, z, nbar, comp'  #ra, dec, z, nz, comp 
            data_list = [(cmass.ra)[zlimit], (cmass.dec)[zlimit], (cmass.z)[zlimit], (cmass.nz)[zlimit], comp[zlimit]]
            data_fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f']
        
        else:
            raise NotImplementedError()

        # write to corrected file 
        output_file = self.file()
        np.savetxt(
                output_file, 
                (np.vstack(np.array(data_list))).T, 
                fmt=data_fmt, 
                delimiter='\t', 
                header=header_str
                ) 

        return None 

    def read(self): 
        ''' 
        Read random catalog data
        '''
        raise ValueError("You're trying to read the random catalog -- don't do it.")

    def datacolumns(self): 
        ''' 
        Data columns for given catalog and correction
        '''
        cat_name = ((self.cat_corr)['catalog'])['name'].lower()

        if cat_name == 'nseries':
            data_cols = ['ra', 'dec', 'z', 'comp']
        elif cat_name == 'qpm': 
            data_cols = ['ra', 'dec', 'z', 'comp']
        elif cat_name == 'cmass': 
            data_cols = ['ra', 'dec', 'z', 'nbar', 'comp']

        return data_cols 

    def datacols_fmt(self): 
        ''' 
        Data format of columns of catalog data
        '''
        cat_name = ((self.cat_corr)['catalog'])['name'].lower()

        if cat_name == 'nseries':
            data_fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f']
        elif cat_name == 'qpm':
            data_fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f']
        elif cat_name == 'cmass': 
            data_fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f']

        return data_fmt 

    def datacols_header(self): 
        ''' 
        Header string that describes data columsn
        '''
        cat_name = ((self.cat_corr)['catalog'])['name'].lower()

        if cat_name == 'nseries':
            hdr_str = 'Columns : ra, dec, z, comp'
        elif cat_name == 'qpm':
            hdr_str = 'Columns : ra, dec, z, comp'
        elif cat_name == 'cmass': 
            hdr_str = 'Columns : ra, dec, z, nbar, comp'

        return hdr_str

if __name__=="__main__":
    cat_corr = {'catalog': {'name': 'qpm'}, 'correction': {'name': 'true'}}
    rand_class = Random(cat_corr)
    print rand_class.file()
    rand_class.build()
