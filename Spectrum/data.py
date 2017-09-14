'''

Code to handle galaxy data for calculating power/bispectrum 

Author(s): ChangHoon Hahn 


'''
import numpy as np
import scipy as sp 
import os.path
import cosmolopy as cosmos
# -- local -- 
import util as UT


# Classes ------------------------------------------------------------
class Data(object): 
    def __init__(self, catinfo): 
        ''' Data class object describes a simulated/observed galaxy catalog. 

        Parameters
        ----------
        catinfo : Dictionary that specifies the details of the sim./obv. catalog.
        ''' 
        self.catinfo = catinfo.copy()    
        self.catalog = (catinfo['catalog'])['name'].lower() 
        
        # cosmology of catalog
        try: 
            self.cosmo = catinfo['cosmology']# type (fid or survey) 
        except KeyError: 
            self.cosmo = 'fid' # type (fid or survey) 

        # data or random 
        try: 
            self.type = catinfo['catalog']['type']
        except KeyError:  # by default it assumes data 
            self.type = 'data' 

        if self.type == 'data': 
            if 'n_mock' in catinfo['catalog'].keys(): 
                self.n_mock = (catinfo['catalog'])['n_mock']
            self.file_name = self.file()
        elif self.type == 'random': 
            self.file_name = self.rfile()
        
        # galaxy properties
        self.ra = None
        self.dec = None
        self.z = None
        self.wfc = None
        self.weight = None 
        self.comp = None

    def read(self): 
        ''' Read in catalog data 
        '''
        if self.type == 'random': 
            raise ValueError("You're trying to read the random catalog -- don't do it.")
        cols = self.datacolumns()
        datah = np.loadtxt(
                self.file_name, 
                skiprows=1, # for the header 
                unpack=True, 
                usecols=range(len(cols[0])))
        # now save the data into the class
        for i_col, col in enumerate(cols[0]): 
            setattr(self, col, datah[i_col])
        return None
    
    def file(self): 
        ''' Given 'catinfo' dictionary specifing the catalog, return 
        name of ASCII catalog file
        '''
        name_comp = self._file_comp() 
        # cosmology string 
        if self.cosmo == 'fid': 
            cosmos_str = '.fidcosmo'
        elif self.cosmo == 'survey': 
            cosmos_str = '.sureycosmo'
        else: 
            raise NotImplementedError() 
        name_comp.insert(-1, cosmos_str)

        return ''.join(name_comp)

    def rfile(self): 
        ''' file name of random catalog file. Everything is quite hardcoded.
        ''' 
        cat_name = self.catalog 
        cat_dict = (self.catinfo)['catalog']
        #corr_name = ((self.catinfo)['correction'])['name'].lower()

        if 'cmass' in cat_name: # CMASS
            if cat_name == 'cmass': # CMASS random catalog 
                file_name = 'cmass-dr12v4-N-Reid.ran.dat'
            elif 'cmasslowz' in cat_name:  # CMASS LOWZ combined random catalog
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
                    '_North.ran.dat'])
            else: 
                raise NotImplementedError()
        
        elif cat_name == 'nseries': # Nseries
            file_name = 'Nseries_cutsky_randoms_50x_redshifts_comp.dat'
        elif cat_name == 'qpm': # QPM 
            file_name = 'a0.6452_rand50x.dr12d_cmass_ngc.vetoed.dat'
        elif cat_name == 'tilingmock': # Tiling Mock 
            file_name = 'randoms-boss5003-icoll012-vetoed.zlim.dat'
        elif cat_name == 'bigmd': # Big MultiDark
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
        return  ''.join([UT.data_dir('data', self.catalog), file_name])

    def cosmology(self): 
        ''' Return `cosmolopy` cosmology dictionary
        '''
        if self.cosmo == 'fid': 
            omega_m = 0.31 
        elif self.cosmo == 'survey': 
            if self.catalog == 'nseries':
                omega_m =  0.286
            elif self.catalog == 'qpm': 
                omega_m = 0.31
            else:
                raise NotImplementedError()
        else: 
            raise NotImplementedError()

        # survey cosmology 
        cosmo = {} 
        cosmo['omega_M_0'] = omega_m 
        cosmo['omega_lambda_0'] = 1.0 - omega_m 
        cosmo['h'] = 0.676
        cosmo = cosmos.distance.set_omega_k_0(cosmo) 
        return cosmo
    
    def zlimits(self): 
        ''' return redshift limit of survey 
        '''
        # redshift limits 
        if self.catalog in ('lasdamasgeo', 'ldgdownnz'):     
            zlim = [0.16, 0.44]
        elif self.catalog in ('tilingmock', 'qpm', 'patchy', 'nseries', 'bigmd'): 
            zlim = [0.43, 0.7]
        elif self.catalog in ('qso_bigmd'): 
            zlim = [0.9, 2.2]
        elif 'bigmd' in self.catalog:             
            zlim = [0.43, 0.7]
        elif 'cmass' in self.catalog:
            if self.catalog == 'cmass': 
                zlim = [0.43, 0.7] 
            elif 'cmasslowz' in self.catalog: 
                if '_high' in self.catalog: 
                    zlim = [0.5, 0.75]
                elif '_low' in self.catalog: 
                    zlim = [0.2, 0.5] 
            else: 
                raise NotImplementedError('CMASS or CMASSLOWZ combined sample')
        else: 
            raise NotImplementedError()
        return zlim 
    
    def datacolumns(self): 
        ''' Data columns for the catalog. 
        Returns [column name, column format, column header]
        '''
        if self.type == 'data': 
            catcol = self._columns((self.catinfo['catalog'])['name'].lower()) 
        elif self.type == 'random': 
            catcol = self._columns((self.catinfo['catalog'])['name'].lower()) 
        return [catcol['cols'], catcol['fmts'], catcol['hdrs']]
   
    def _build(self): 
        ''' Construct the 'defaut' galaxy catalog by reading in 
        a variety of data.

        Notes
        -----
        * Tiling Mock : Reads in true mock catalogs and imposes predetermined redshift 
        limits on them and attaches a flag mostly hardcoded since it's a simple procedure 
        * QPM: Handles messy weights 
        * PATCHY: Returns mocks with only necessary columns and w_fc = 1
        * Everything is very hardcoded
        '''
        if self.catalog == 'nseries': 
            # read in original files from Jeremy Tinker and adjust them to make them
            # easier to use for fiber collisions
            # read rdzw file 
            data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
            orig_file = ''.join([
                data_dir, 
                'CutskyN', str(catalog['n_mock']), '.rdzwc'
                ]) 
            orig_ra, orig_dec, orig_z, orig_wfc, orig_zupw, orig_upw_index = np.loadtxt(
                    orig_file, 
                    unpack = True, 
                    usecols=[0,1,2,4,5,6]
                    )
        
            # file with completeness
            mask_file = ''.join([
                data_dir, 
                'CutskyN', str(catalog['n_mock']), '.mask_info'
                ]) 

            orig_wcomp = np.loadtxt(
                    mask_file, 
                    unpack = True, 
                    usecols = [0]
                    ) 

            # true wfc 
            true_wfc = np.array([ 1.0 for i in range(len(orig_wfc)) ]) 
            
            # write to file 
            data_list = [orig_ra, orig_dec, orig_z, true_wfc, orig_wcomp, orig_zupw, orig_upw_index]
        elif self.catalog == 'qpm': # QPM 
            # import original true data 
            orig_file = ''.join([UT.data_dir('data', self.catalog),
                'a0.6452_', str("%04d" % self.catinfo['catalog']['n_mock']), '.dr12d_cmass_ngc.rdz']) 
            orig_data = np.loadtxt(orig_file) 

            orig_info_file = orig_file+'.info'
            orig_info = np.loadtxt(orig_info_file)    # gal_id, comp, z_real, z_red, mass_halo, flag_sta, id_halo
        
            # consistency issue with #46
            if catalog['n_mock'] in (46, 52, 53, 54, 56, 61, 707, 756, 794, 819, 831, 835, 838): 
                orig_veto_file = ''.join(['/mount/riachuelo1/hahn/data/QPM/dr12d/', 
                    'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.veto']) 
            else:
                orig_veto_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
                    'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.veto']) 
            orig_veto = np.loadtxt(orig_veto_file) 
            n_gal = len(orig_veto)
       
            # assign RA, Dec, z, and w_veto 
            true_ra     = orig_data[:,0]
            true_dec    = orig_data[:,1]
            true_z      = orig_data[:,2]
            true_wfkp   = orig_data[:,3]
            true_wfc = np.repeat(1.0, n_gal)    # fiber collisions weights are all 1 for true

            # check to make sure that the redshifts correspond btw rdz file and rdz.info file 
            if np.max(np.abs(orig_info[:,3]/true_z-1.0)) > 10**-5: 
                raise ValueError('redshifts between the data file and info file dont match') 

            true_comp = orig_info[:,1]         # completness weights

            # remove veto mask 
            # Only keep galaxies with veto = 0 (for veto values in .veto file) 
            vetomask = np.where(orig_veto == 0)            
            data_list = [true_ra[vetomask], true_dec[vetomask], true_z[vetomask], true_wfc[vetomask], true_comp[vetomask]]
        
        elif catalog['name'].lower() == 'bigmd':                
            # Big MultiDark  
            P0 = 20000.0
            # read original file
            data_dir = '/mount/riachuelo1/hahn/data/BigMD/'
            if 'version' in catalog.keys():
                raise NotImplementedError
                if catalog['version'] == 'blah':
                    pass
                orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-wcp-veto.dat'])  # hardcoded
                orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-standardHAM-veto.dat'])
                orig_file = ''.join([data_dir, 'bigMD-cmass-dr12v4-RST-quadru-veto.dat'])
            else:       # default 
                orig_file = ''.join([data_dir, 'BigMD-cmass-dr12v4-RST-standHAM-Vpeak-veto.dat'])

            orig_ra, orig_dec, orig_z, orig_wfkp, orig_veto, orig_wfc = np.loadtxt(
                    orig_file, 
                    unpack=True, 
                    usecols=[0,1,2,3,4,5])
        
            # true wfc = 1 for all galaxies 
            true_wfc = np.array([ 1.0 for i in range(len(orig_wfc)) ]) 
            nbar = (1.0/P0) * (1.0/orig_wfkp - 1.0) 

            vetomask = np.where(orig_veto == 1)     # if veto = 1 then keep; otherwise discard 
            
            # write to file 
            data_list = [
                    orig_ra[vetomask], 
                    orig_dec[vetomask], 
                    orig_z[vetomask], 
                    nbar[vetomask], 
                    true_wfc[vetomask]
                    ]
    
        elif catalog['name'].lower() == 'tilingmock': 
            # read in original file and impose redshift limits 
            orig_ra, orig_dec, orig_z, orig_w = np.loadtxt(
                    '/mount/riachuelo1/hahn/data/tiling_mocks/', 
                    'cmass-boss5003sector-icoll012.dat', 
                    unpack=True, 
                    usecols=[0,1,2,3]) 
            
            zlow, zhigh = self.survey_zlimits()
            zlim = np.where((orig_z > zlow) & (orig_z < zhigh))
            
            true_wfc = np.repeat(1.0, len(orig_ra))

            data_list = [orig_ra[zlim], orig_dec[zlim], orig_z[zlim], true_wfc[zlim]]

        elif catalog['name'].lower() == 'lasdamasgeo':          
            orig_true_file = ''.join(['/mount/chichipio2/rs123/MOCKS/LRGFull_zm_geo/gaussian/zspace/', 
                'sdssmock_gamma_lrgFull_zm_oriana', str("%02d" % catalog['n_mock']), catalog['letter'], '_no.rdcz.dat']) 
            orig_true_data = np.loadtxt(orig_true_file, unpack=True, usecols=[0,1,2])         # ra, dec, ***CZ***

            true_ra = orig_true_data[0]
            true_dec = orig_true_data[1]
            true_z = orig_true_data[2]/299800.0         # convert cz to z

            true_weight = np.array([1.0 for j in range(len(true_z))])   # no weights for true (all 1) 

            header_str = "Columns : ra, dec, z, weight"
            data_list = [true_ra, true_dec, true_z, true_weight]
            data_fmt = ['%10.5f', '%10.5f', '%10.5f', '%10.5f']

        elif catalog['name'].lower() == 'patchy':       # PATCHY mocks ------------------------
            # read original mock data 
            orig_file = ''.join(['/mount/riachuelo1/hahn/data/PATCHY/dr12/v6c/', 
                'Patchy-Mocks-DR12CMASS-N-V6C-Portsmouth-mass_', 
                str("%04d" % catalog['n_mock']), '.dat']) 

            # ra, dec, z, nbar, wfc, veto 
            orig_ra, orig_dec, orig_z, orig_nbar, orig_veto = np.genfromtxt(orig_file, 
                    unpack=True, usecols=[0, 1, 2, 4, 6]) 
            n_gal = len(orig_ra) 
            
            new_wfc = np.array([1.0 for i in range(n_gal)])     # w_fc = 1.0 for true 

            vetomask = (orig_veto == 1)            # only keep galaxies with w_veto = 1
           
            header_str = "Columns : ra, dec, z, nz, wfc" 
            data_list = [orig_ra[vetomask], orig_dec[vetomask], orig_z[vetomask], 
                        orig_nbar[vetomask], new_wfc[vetomask]] 
            data_fmt=['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f']
        

        else: 
            raise NameError('not yet coded') 
            
        # data columns, format and header are determined in util.catalog
        # and are catalog specific
        data_fmt = self.datacols_fmt()
        header_str = self.datacols_header()
        
        # write to file 
        np.savetxt(
                self.file(), 
                (np.vstack(np.array(data_list))).T, 
                fmt=data_fmt, 
                delimiter='\t', 
                header=header_str) 
        return None 

    def _rbuild(self): 
        ''' build the random catalogs from original survey data. 
        Everything is hardcoded. 
        '''
        catdict = (self.catinfo)['catalog']
        corrdict = (self.catinfo)['correction']
        cat_name = self.catalog 

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

    def _file_comp(self): 
        ''' Given 'catinfo' dictionary specifing the catalog, return 
        naming components of the ASCII catalog file. Messy and therefore 
        hidden away!
        '''
        cat = self.catalog 
        dir = UT.data_dir('data', cat)
        
        if cat == 'nseries': # Nseries
            file_beg = ''.join(['CutskyN', str(self.n_mock)])
            file_end = '.dat'
        elif 'cmass' in cat: # CMASS
            if cat == 'cmass': 
                file_beg = 'cmass-dr12v4-N-Reid'
                file_end = '.dat'
            elif 'cmasslowz' in cat:
                file_beg = ''.join(['galaxy_DR12v5_', self.catalog_name.upper(), '_North'])
                file_end = '.dat'
        elif cat == 'qpm':            # QPM
            file_beg = ''.join(['a0.6452_', str("%04d" % self.n_mock), '.dr12d_cmass_ngc.vetoed'])
            file_end = '.dat'
        elif cat == 'tilingmock':     # Tiling Mock 
            file_beg = 'cmass-boss5003sector-icoll012'
            file_end = '.dat'
        elif cat == 'bigmd':          # Big MultiDark
            if 'version' in catinfo['catalog'].keys(): 
                if catinfo['catalog']['version'] == 'nowiggles': 
                    file_beg = 'BigMD-cmass-dr12v4-nowiggle-veto'
                    file_end = '.dat'
                else: 
                    raise NotImplementedError
            else: 
                file_beg = 'BigMD-cmass-dr12v4-RST-standHAM-Vpeak_vetoed'
                file_end = '.dat'
        elif cat == 'qso_bigmd': 
            if 'version' in catinfo['catalog'].keys(): 
                if catinfo['catalog']['version'] == 'evo': 
                    file_beg = 'QSO-BigMD-1481.43deg2-bin0_v1.0-evo'
                    file_end = '.dat'
                elif catinfo['catalog']['version'] == 'noevo': 
                    file_beg = 'QSO-BigMD-1481.43deg2-bin0_v1.0-noevo'
                    file_end = '.dat'
                elif catinfo['catalog']['version'] == 'eboss': 
                    file_beg = 'eboss_v1.0-QSO-NS-eboss_v1.0-bin0'
                    file_end = '.dat'
                elif catinfo['catalog']['version'] == 'ebossv1.5': 
                    file_beg = 'eboss_v1.5-QSO'
                    file_end = '.dat'
                elif catinfo['catalog']['version'] == 'ebossnew': 
                    file_beg = 'ebossQSOs_y1_comp_cut_0.5_double_weight_col_Z'
                    file_end = '.dat'
                elif 'jackknife' in catinfo['catalog']['version']: 
                    n_jack = catinfo['catalog']['version'].split('jackknife')[-1]
                    file_beg = ''.join(['QSO-bin0-jackknife', str(n_jack)])
                    file_end = '.dat'
                elif 'v2' in catinfo['catalog']['version']: 
                    if 'z' in catinfo['catalog']['version']: 
                        file_beg = 'BigMDPL-QSOZ'
                    elif 'nsat' in catinfo['catalog']['version']: 
                        file_beg = 'BigMDPL-QSO-NSAT'
                    else: 
                        file_beg = 'BigMDPL-QSO'
                    file_end = '.dat'
        else: 
            raise NotImplementedError
        return [dir, file_beg, file_end]
    
    def _columns(self, catname):
        ''' Given catalog name, return details on the data columns
        '''
        # hardcoded dictionary of catlaog columsn 
        catcol = {
                'nseries': {
                    'cols': ['ra', 'dec', 'z', 'wfc', 'comp', 'zupw', 'upw_index'], 
                    'fmts': ['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%i'], 
                    'hdrs': "Column : ra, dec, z, wfc, comp, z_upw, upw_index"
                    }, 
                'qpm': {
                    'cols': ['ra', 'dec', 'z', 'wfc', 'comp'], 
                    'fmts': ['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'],
                    'hdrs': "Column : ra, dec, z, wfc, comp"
                    }, 
                'tilingmock': {
                    'cols': ['ra', 'dec', 'z', 'wfc'], 
                    'fmts': ['%10.5f', '%10.5f', '%10.5f', '%10.5f'],
                    'hdrs': "Column : ra, dec, z, wfc"
                    }, 
                'bigmd': {
                    'cols': ['ra', 'dec', 'z', 'wfc'], 
                    'fmts': ['%10.5f', '%10.5f', '%10.5f', '%10.5f'],
                    'hdrs': "Columns : ra, dec, z, wfc"
                    },
                'qso_bigmd': {
                    'cols': ['ra', 'dec', 'z', 'nbar', 'wfc'], 
                    'fmts': ['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f'],
                    'hdrs': "Columns : ra, dec, z, 'nbar', wfc"
                    },
                'cmass': {
                    'cols': ['ra', 'dec', 'z', 'nbar', 'wsys', 'wnoz', 'wfc', 'comp'], 
                    'fmts': ['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f', '%10.5f', '%10.5f', '%10.5f'], 
                    'hdrs': 'Columns : ra, dec, z, nbar, w_systot, w_noz, w_cp, comp'
                    }
                }
        if catname not in catcol.keys(): 
            raise ValueError("catalog name is not in catalog_column dictionary!") 
        else: 
            return catcol[catname] 
    
    def _rcolumns(self, catname):
        ''' Given catalog name, return details on the data columns of randoms
        '''
        # hardcoded dictionary of catlaog columsn 
        catcol = {
                'nseries': {
                    'cols': ['ra', 'dec', 'z', 'comp'], 
                    'fmts': ['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f', '%i'], 
                    'hdrs': "Column : ra, dec, z, comp"
                    }, 
                'qpm': {
                    'cols': ['ra', 'dec', 'z', 'comp'],
                    'fmts': ['%10.5f', '%10.5f', '%10.5f', '%10.5f', '%10.5f'],
                    'hdrs': "Column : ra, dec, z, comp"
                    }, 
                'bigmd': {
                    'cols': ['ra', 'dec', 'z'],
                    'fmts': ['%10.5f', '%10.5f', '%10.5f'],
                    'hdrs': "Columns : ra, dec, z"
                    },
                'qso_bigmd': {
                    'cols': ['ra', 'dec', 'z', 'nbar'], 
                    'fmts': ['%10.5f', '%10.5f', '%10.5f', '%10.5f'],
                    'hdrs': "Columns : ra, dec, z, 'nbar'"
                    },
                'cmass': {
                    'cols': ['ra', 'dec', 'z', 'nbar', 'comp'], 
                    'fmts': ['%10.5f', '%10.5f', '%10.5f', '%.5e', '%10.5f'],
                    'hdrs': 'Columns : ra, dec, z, nbar, comp'
                    }
                }
        if catname not in catcol.keys(): 
            raise ValueError("catalog name is not in catalog_column dictionary!") 
        else: 
            return catcol[catname] 


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
