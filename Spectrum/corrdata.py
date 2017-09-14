''' 

Code to handle corrected galaxy data 

'''
import warnings
# --- Data class ---
from data import Data

class CorrData(Data): 
    def __init__(self, catinfo, **kwargs): 
        ''' Child class of Data describing CORRECTED galaxy catalog of 
        simulations or BOSS data. Corresponds to an ASCII file with 
        galaxy catalog.
        '''
        # correction class dictionary 
        try: 
            self.corrclass_dict.keys()
        except AttributeError:
            self.corrclass_dict = None 
       
        if self.corrclass_dict is not None:
            corr_name = cat_corr['correction']['name']

            if corr_name in self.corrclass_dict.keys():
                self.corrclass = self.corrclass_dict[corr_name](catinfo, **kwargs)
            else: 
                raise ValueError
        else: 
            pass
            #raise NotImplementedError 

        super(CorrData, self).__init__(catinfo, **kwargs)
        self.file_name = self.file()
        
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
       
'''
if __name__ == '__main__':
    for i_mock in xrange(1,85): 
        cat_corr = {
                'catalog': {'name': 'qpm', 'n_mock': i_mock}, 
                'correction': {'name': 'true'}
                }
        corrclass = Data(cat_corr, clobber=True)
        corrclass.build()
'''
