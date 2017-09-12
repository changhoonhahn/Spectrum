'''



'''

class Catalog(object): 
    def __init__(self, cat_corr): 
        ''' 
        Class describing simulated/data catalogs 
        '''
        self.cat_corr = cat_corr.copy()
        self.catalog_name = (cat_corr['catalog'])['name'].lower()

        if self.catalog_name in ('nseries', 'qpm', 'bigmd', 'qso_bigmd'): 
            self.cat_col_dict_key = self.catalog_name
        elif 'cmass' in self.catalog_name: 
            self.cat_col_dict_key = 'cmass' 
        else: 
            raise NotImplementedError()
    
        # hardcoded dictionary of catlaog columsn 
        self.cat_col_dict = {
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

    def file(self): 
        ''' File elements that pertain to the simulated/data catalog
        '''
        cat = self.cat_corr['catalog']
        
        if self.catalog_name == 'nseries':          # Nseries

            data_dir = '/mount/riachuelo1/hahn/data/Nseries/'
            file_beg = ''.join(['CutskyN', str(cat['n_mock'])])
            file_end = '.dat'

        elif 'cmass' in self.catalog_name:          # CMASS

            if self.catalog_name == 'cmass': 
                data_dir = '/mount/riachuelo1/hahn/data/CMASS/'
                file_beg = 'cmass-dr12v4-N-Reid'
                file_end = '.dat'

            elif 'cmasslowz' in self.catalog_name:
                data_dir = '/mount/riachuelo1/hahn/data/CMASS/dr12v5/'
                file_beg = ''.join(['galaxy_DR12v5_', self.catalog_name.upper(), '_North'])
                file_end = '.dat'

        elif self.catalog_name == 'qpm':            # QPM

            data_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/'
            file_beg = ''.join(['a0.6452_', str("%04d" % cat['n_mock']), '.dr12d_cmass_ngc.vetoed'])
            file_end = '.dat'

        elif self.catalog_name == 'tilingmock':     # Tiling Mock 

            data_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'
            file_beg = 'cmass-boss5003sector-icoll012'
            file_end = '.dat'

        elif self.catalog_name == 'bigmd':          # Big MultiDark

            if 'version' in self.cat_corr['catalog'].keys(): 
                if self.cat_corr['catalog']['version'] == 'nowiggles': 
                    data_dir = '/mount/riachuelo1/hahn/data/BigMD/nowiggles/'
                    file_beg = 'BigMD-cmass-dr12v4-nowiggle-veto'
                    file_end = '.dat'
                else: 
                    raise NotImplementedError
            else: 
                data_dir = '/mount/riachuelo1/hahn/data/BigMD/'
                file_beg = 'BigMD-cmass-dr12v4-RST-standHAM-Vpeak_vetoed'
                file_end = '.dat'
    
        elif self.catalog_name == 'qso_bigmd': 
            data_dir = '/mount/riachuelo1/hahn/data/BigMD/'
            if 'version' in self.cat_corr['catalog'].keys(): 
                if self.cat_corr['catalog']['version'] == 'evo': 
                    file_beg = 'QSO-BigMD-1481.43deg2-bin0_v1.0-evo'
                    file_end = '.dat'
                elif self.cat_corr['catalog']['version'] == 'noevo': 
                    file_beg = 'QSO-BigMD-1481.43deg2-bin0_v1.0-noevo'
                    file_end = '.dat'
                elif self.cat_corr['catalog']['version'] == 'eboss': 
                    file_beg = 'eboss_v1.0-QSO-NS-eboss_v1.0-bin0'
                    file_end = '.dat'
                elif self.cat_corr['catalog']['version'] == 'ebossv1.5': 
                    file_beg = 'eboss_v1.5-QSO'
                    file_end = '.dat'
                elif self.cat_corr['catalog']['version'] == 'ebossnew': 
                    file_beg = 'ebossQSOs_y1_comp_cut_0.5_double_weight_col_Z'
                    file_end = '.dat'
                elif 'jackknife' in self.cat_corr['catalog']['version']: 
                    n_jack = self.cat_corr['catalog']['version'].split('jackknife')[-1]
                    file_beg = ''.join(['QSO-bin0-jackknife', str(n_jack)])
                    file_end = '.dat'
                elif 'v2' in self.cat_corr['catalog']['version']: 
                    if 'z' in self.cat_corr['catalog']['version']: 
                        file_beg = 'BigMDPL-QSOZ'
                    elif 'nsat' in self.cat_corr['catalog']['version']: 
                        file_beg = 'BigMDPL-QSO-NSAT'
                    else: 
                        file_beg = 'BigMDPL-QSO'
                    file_end = '.dat'
        else: 
            raise NotImplementedError

        return [data_dir, file_beg, file_end]

    def cosmo(self): 
        """ Survey cosmology. Note this is different than the fiducial cosmology
        used to compute powerspectrum 
        """

        if self.catalog_name == 'nseries': 
            omega_m =  0.286
        elif self.catalog_name == 'qpm': 
            omega_m = 0.31
        else: 
            raise NotImplementedError()
        
        cosmo = {} 
        cosmo['omega_M_0'] = omega_m 
        cosmo['omega_lambda_0'] = 1.0 - omega_m 
        cosmo['h'] = 0.676
        cosmo = cosmos.distance.set_omega_k_0(cosmo) 
        self.cosmo = cosmo 

        return cosmo 

    def survey_zlimits(self): 
        ''' Redshift limits of survey 
        '''
        
        # survey redshift limits 
        if self.catalog_name in ('lasdamasgeo', 'ldgdownnz'):     
            survey_zmin, survey_zmax = 0.16, 0.44

        elif self.catalog_name in ('tilingmock', 'qpm', 'patchy', 'nseries', 'bigmd'): 
            survey_zmin, survey_zmax = 0.43, 0.7 

        elif self.catalog_name in ('qso_bigmd'): 
            survey_zmin, survey_zmax = 0.9, 2.2 

        elif 'bigmd' in self.catalog_name:             
            survey_zmin, survey_zmax = 0.43, 0.7    

        elif 'cmass' in self.catalog_name:             

            if self.catalog_name == 'cmass': 
                survey_zmin, survey_zmax = 0.43, 0.7  

            elif 'cmasslowz' in self.catalog_name: 

                if '_high' in self.catalog_name: 
                    survey_zmin, survey_zmax = 0.5, 0.75    

                elif '_low' in self.catalog_name: 
                    survey_zmin, survey_zmax = 0.2, 0.5

            else: 
                raise NotImplementedError('CMASS or CMASSLOWZ combined sample')

        else: 
            raise NotImplementedError('Mock Catalog not included')

        return [survey_zmin, survey_zmax]

    def datacolumns(self): 
        ''' Columns of catalog data
        '''
        return (self.cat_col_dict[self.cat_col_dict_key])['cols']

    def datacols_fmt(self): 
        ''' Data format of columns of catalog data
        '''
        return (self.cat_col_dict[self.cat_col_dict_key])['fmts']

    def datacols_header(self): 
        ''' Header string that describes data columsn
        '''
        return (self.cat_col_dict[self.cat_col_dict_key])['hdrs']
