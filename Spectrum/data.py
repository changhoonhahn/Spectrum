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
    def __init__(self, catinfo, cosmo='fid', **kwargs): 
        ''' Data class object describes a simulated/observed galaxy catalog. 

        Parameters
        ----------
        catinfo : Dictionary that specifies the details of the sim./obv. catalog.
        ''' 
        self.catinfo = catinfo.copy()    
        self.catalog = (catinfo['catalog'])['name'].lower() 
        self.kwargs = kwargs.copy() 
        
        # cosmology of catalog
        self.cosmo = cosmo  # type (fid or survey) 

        # galaxy properties
        self.ra = None
        self.dec = None
        self.z = None
        self.wfc = None
        self.weight = None 
        self.comp = None

        self.file_name = self.file()

    def read(self): 
        ''' Read in catalog data 
        '''
        cols = self.datacolumns()
        datah = np.loadtxt(
                self.file_name, 
                skiprows=1, # for the header 
                unpack=True, 
                usecols=range(len(cols[0])))
        # now save the data into the class
        for i_col, col in enumerate(col[0]): 
            setattr(self, col, datah[i_col])
        return None
    
    def file(self): 
        ''' Given 'catinfo' dictionary specifing the catalog, return 
        name of ASCII catalog file
        '''
        name_comp = self._file_comp(self.catinfo) 
        # cosmology string 
        if self.cosmo == 'fid': 
            cosmos_str = '.fidcosmo'
        elif self.cosmo == 'survey': 
            cosmos_str = '.sureycosmo'
        else: 
            raise NotImplementedError() 
        name_comp.insert(-1, cosmos_str)

        return ''.join(file_list)

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
        catcol = self._catalog_columns((self.catinfo['catalog'])['name'].lower()) 
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

        elif catalog['name'].lower() == 'qpm':      # QPM 
            # import original true data 
            orig_file = ''.join([
                '/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
                'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.rdz']) 
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

    def _file_comp(self): 
        ''' Given 'catinfo' dictionary specifing the catalog, return 
        naming components of the ASCII catalog file. Messy and therefore 
        hidden away!
        '''
        cat = self.catalog 
        dir = UT.data_dir('data', cat)
        
        if cat == 'nseries': # Nseries
            file_beg = ''.join(['CutskyN', str(cat['n_mock'])])
            file_end = '.dat'
        elif 'cmass' in cat: # CMASS
            if cat == 'cmass': 
                file_beg = 'cmass-dr12v4-N-Reid'
                file_end = '.dat'
            elif 'cmasslowz' in cat:
                file_beg = ''.join(['galaxy_DR12v5_', self.catalog_name.upper(), '_North'])
                file_end = '.dat'
        elif cat == 'qpm':            # QPM
            file_beg = ''.join(['a0.6452_', str("%04d" % cat['n_mock']), '.dr12d_cmass_ngc.vetoed'])
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
    
    def _catalog_columns(self, catname):
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
        if catname not in catcol.key(): 
            raise ValueError("catalog name is not in catalog_column dictionary!") 
        else: 
            return catcol[catname] 
