'''

Random data for calculating power/bispectrum

'''
import numpy as np
import time 
# --- Local ---
from ChangTools.fitstables import mrdfits
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

    def build(self): 
        ''' 
        Completely replace parent class build attribute in order
        to ruild the random catalogs from original survey data. 
        Everything is hardcoded. 
        '''
        catdict = (self.cat_corr)['catalog']
        corrdict = (self.cat_corr)['correction']
        cat_name = catdict['name'].lower()

        if corrdict['name'].lower() == 'noweight': 
            NoweightRandoms(self.cat_corr) 
            return None

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

        elif cat_name == 'bigmd':
            # Big MultiDark 
            P0 = 20000.0
            # original random catalog 
            data_dir = '/mount/riachuelo1/hahn/data/BigMD/'
            if 'version' in catdict.keys():
                orig_rand_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-wcp-veto.ran']) 
                orig_rand_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-standardHAM-veto.ran']) 
                orig_rand_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-quadru-veto.ran']) 
            else: 
                orig_rand_file = ''.join([
                    data_dir, 
                    'BigMD-cmass-dr12v4-RST-standHAM-Vpeak-veto.ran'
                    ])

            # RA, Decl, Redhsift, veto  
            ra, dec, z, wfkp, veto = np.loadtxt(orig_rand_file, unpack=True, usecols=[0,1,2,3,4]) 
            #nbar = (1.0 / P0) * (1.0/wfkp - 1.0)    # nbar(z) 

            vetomask = np.where(veto == 1)  # impose vetomask 
        
            header_str = 'Columns : ra, dec, z'
            data_list = [ra[vetomask], dec[vetomask], z[vetomask]]
            data_fmt=['%10.5f', '%10.5f', '%10.5f']
        
        elif cat_name == 'qso_bigmd':
            # original random catalog 
            P0 = 20000.
            data_dir = '/mount/riachuelo1/hahn/data/BigMD/'
            if catdict['version'] == 'evo':
                orig_rand_file = ''.join([data_dir, 'QSO-BigMD-1481.43deg2-bin0_v1.0_evo.ran'])
            elif catdict['version'] == 'noevo': 
                orig_rand_file = ''.join([data_dir, 'QSO-BigMD-1481.43deg2-bin0_v1.0_noevo.ran'])
            elif catdict['version'] == 'eboss': 
                orig_rand_file = ''.join([data_dir, 'eboss_v1.0-QSO-NS-eboss_v1.0_bin0.ran'])
            elif catdict['version'] == 'ebossv1.5': 
                P0 = 6000.
                orig_rand_file = ''.join([data_dir, 'eboss_v1.5-QSO-eboss_v1.5.ran'])
            elif catdict['version'] == 'ebossnew': 
                orig_rand_file = ''.join([data_dir, 
                    'rand20_y1_comp_cut_0.5_double_weight_col_Z.dat'])
            elif 'jackknife' in catdict['version']: 
                n_jack = catdict['version'].split('jackknife')[-1]
                orig_rand_file = ''.join([data_dir, 'QSO-bin0_', str(n_jack), '.ran']) 
            elif 'v2' in catdict['version'] : 
                P0 = 6000.
                if 'z' in catdict['version']: 
                    orig_rand_file = ''.join([data_dir, 'BigMDPL-QSOZ.ran']) 
                elif 'nsat' in catdict['version'] : 
                    orig_rand_file = ''.join([data_dir, 'BigMDPL-QSO-NSAT.ran']) 
                else: 
                    orig_rand_file = ''.join([data_dir, 'BigMDPL-QSO.ran']) 
            else: 
                raise NotImplementedError

            # RA, Decl, Redhsift, wfkp 
            ra, dec, z, wfkp = np.loadtxt(orig_rand_file, unpack=True, usecols=[0,1,2,3]) 
            nbar = (1.0 / P0) * (1.0/wfkp - 1.0)    # nbar(z) 

            header_str = 'Columns : ra, dec, z, nbar'
            data_list = [ra, dec, z, nbar]
            data_fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e']

        elif cat_name == 'tilingmock':  # Tiling Mock  
            
            orig_file = ''.join([
                '/mount/riachuelo1/hahn/data/tiling_mocks/', 
                'randoms-boss5003-icoll012-vetoed.dat'
                ])
            orig_ra, orig_dec, orig_z, orig_w = np.loadtxt(
                    orig_file, unpack=True, usecols=[0,1,2,3])

            zlim = np.where((orig_z > 0.43) & (orig_z < 0.7))

            header_str = 'Columns : ra, dec, z, wfc'
            data_list = [orig_ra[zlim], orig_dec[zlim], orig_z[zlim], orig_w[zlim]]
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
        elif cat_name == 'bigmd': 
            data_cols = ['ra', 'dec', 'z']
        elif cat_name == 'qso_bigmd': 
            data_cols = ['ra', 'dec', 'z', 'nbar']
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
        elif cat_name == 'bigmd': 
            data_fmt=['%10.5f', '%10.5f', '%10.5f']
        elif cat_name == 'qso_bigmd': 
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
        elif cat_name == 'bigmd':
            hdr_str = 'Columns : ra, dec, z'
        elif cat_name == 'qso_bigmd':
            hdr_str = 'Columns : ra, dec, z, comp'
        elif cat_name == 'cmass': 
            hdr_str = 'Columns : ra, dec, z, nbar, comp'

        return hdr_str

def NoweightRandom(cat_corr): 
    ''' Hack to generate downsampled no weight randoms
    '''
    # original randoms
    true_cat_corr = cat_corr.copy()
    true_cat_corr['correction'] = {'name': 'true'}

    orig_rand = Random(true_cat_corr)
    orig_file = orig_rand.file()
    orig_cols = orig_rand.datacolumns()
    orig_data = np.loadtxt(orig_file, skiprows=1, unpack=True)

    z_col = orig_cols.index('z')
    orig_z = orig_data[z_col]

    # import f_noweight file 
    if cat_corr['catalog']['name'] == 'bigmd':  
        dir = '/mount/riachuelo1/hahn/data/BigMD/' 
    elif cat_corr['catalog']['name'] == 'nseries':
        dir = '/mount/riachuelo1/hahn/data/Nseries/' 
    elif cat_corr['catalog']['name'] == 'qpm': 
        dir ='/mount/riachuelo1/hahn/data/QPM/dr12d/'
    fnow_file = ''.join([dir, 
        'fnoweight_nbarz.', cat_corr['catalog']['name'], '.noweight', 
        '.dat'])
    zmid, zlow, zhigh, f_now = np.loadtxt(fnow_file, unpack=True, usecols=[0,1,2,3]) 

    downsample_i = [] 
    for iz in range(len(zmid)): 
        zbin = np.where(
                (orig_z > zlow[iz]) & 
                (orig_z <= zhigh[iz])
                )

        ngal_bin0 = len(zbin[0])
        ngal_bin_now = round(np.float(ngal_bin0) * f_now[iz]) 

        ngal_downsample = ngal_bin0 - ngal_bin_now 
        print zmid[iz], ngal_downsample
        if ngal_downsample > 0: 
            downsample_i += list(np.random.choice(zbin[0], int(ngal_downsample)))
        else: 
            continue
    keep = list(set(range(len(orig_z))) - set(downsample_i))
    #keep = [ i for i in range(len(orig_z)) if i not in downsample_i ]
    print 'Ngal = ', len(orig_z), 'Ngal downsampled = ', len(keep)

    data_list = [] 
    for i_col in range(len(orig_cols)):  
        data_list.append((orig_data[i_col])[keep]) 
    
    now_rand = Random(cat_corr)
    output_file = now_rand.file()
    print output_file 
    np.savetxt(
            output_file, 
            (np.vstack(np.array(data_list))).T, 
            fmt=orig_rand.datacols_fmt(), 
            delimiter='\t', 
            header=orig_rand.datacols_header()
            ) 
    return None




if __name__=="__main__":
    for cat in ['bigmd', 'nseries', 'qpm']: 
        cat_corr = {'catalog': {'name': cat}, 'correction': {'name': 'noweight'}}
        NoweightRandom(cat_corr)
        #rand_class = Random(cat_corr)
        #print rand_class.file()
        #rand_class.build()
