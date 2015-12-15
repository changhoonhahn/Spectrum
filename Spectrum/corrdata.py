''' 

Code to handle corrected galaxy data 

'''
import warnings
# --- Data class ---
from data import Data

class CorrData(Data): 
    def __init__(self, cat_corr, **kwargs): 
        '''
        Child class of Data describing CORRECTED galaxy catalog of 
        simulations or BOSS data. Corresponds to an ASCII file with 
        galaxy catalog.

        Parameters
        ----------
        cat_corr :  
            Catalog correction Dictionary 
        
        '''
        # correction class dictionary 
        try: 
            self.corrclass_dict.keys()
        except AttributeError:
            self.corrclass_dict = None 
       
        if self.corrclass_dict is not None:
            corr_name = cat_corr['correction']['name']

            if corr_name not in self.corrclass_dict.keys():
                cat_corr['correction'] = {'name': 'default'}
            else: 
                self.corrclass = self.corrclass_dict[corr_name](cat_corr, **kwargs)
        else: 
            try: 
                cat_corr['correction']['name']
            except KeyError: 
                cat_corr['correction'] = {'name': 'default'}
        
        super(CorrData, self).__init__(cat_corr, **kwargs)
        
    def build(self): 
        '''
        Calculate galaxy catalog
        '''
        corr_name = self.cat_corr['correction']['name'].lower()
        if corr_name not in self.corrclass_dict.keys(): 
            super(CorrData, self).build()
        else: 
            (self.corrclass).build()

        return None 

    def file(self): 
        ''' 
        Name of ASCII file of Data/Random catalogy
        '''
        try: 
            file_name = (self.corrclass).file()
            return file_name

        except AttributeError: 
            return super(CorrData, self).file()
        
if __name__ == '__main__':

    for i_mock in xrange(1,85): 
        for corr in ['true', 'upweight', 'photoz']:
            cat_corr = {
                    'catalog': {'name': 'nseries', 'n_mock': i_mock}, 
                    'correction': {'name': 'photoz'}
                    }
            corrclass = Data('data', cat_corr, clobber=True)
            corrclass.build()