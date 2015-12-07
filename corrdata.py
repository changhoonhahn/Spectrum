''' 

Code to handle corrected galaxy data 

'''
import warnings
# --- Data class ---
from data import Data
# --- Corrections ---
#from corrections.true import TrueCorr
#from corrections.dlospeak import DlospeakCorr
#from corrections.fibcollided import UpweightCorr
#from corrections.dlospeak_env import DlospeakEnvCorr
#from corrections.dlospeak_flex import DlospeakFlexCorr
#from corrections.dlospeak_known import DlospeakKnownCorr
#from corrections.dlospeak_photoz import DlospeakPhotozCorr
#from corrections.dlospeak_shuffle import DlospeakShuffleCorr
#from corrections.dlospeak_tailonly import DlospeakTailonlyCorr
#from corrections.dlospeak_peakonly import DlospeakPeakonlyCorr
#from corrections.photoz_corr import PhotozCorr
#from corrections.fourier_tophat import FourierTophatCorr

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
        corrclass_dict = None
        #'true': TrueCorr,
        #'upweight': UpweightCorr, 
        #'photoz': PhotozCorr,
        #'dlospeak': DlospeakCorr, 
        #'dlospeakenv': DlospeakEnvCorr, 
        #'dlospeakphotoz': DlospeakPhotozCorr,
        #'dlospeakknown': DlospeakKnownCorr,
        #'dlospeak.flex': DlospeakFlexCorr,
        #'dlospeak.shuffle': DlospeakShuffleCorr,
        #'dlospeak.tailonly': DlospeakTailonlyCorr, 
        #'dlospeak.peakonly': DlospeakPeakonlyCorr,
        #'fourier_tophat': FourierTophatCorr
        #}
        
        if corrclass_dict is not None: 
            corr_name = (self.cat_corr['correction'])['name']
            if corr_name not in corrclass_dict.keys():
                raise NameError()

            self.corrclass = corrclass_dict[corr_name](cat_corr, **kwargs)
        else: 
            pass
        
        super(CorrData, self).__init__(cat_corr, **kwargs)
        
    def build(self): 
        '''
        Calculate galaxy catalog
        '''
        corr_name = self.cat_corr['correction']['name'].lower()
        if corr_name not in corrclass_dict.keys(): 
            super(CorrData, self).build()
        else: 
            (self.corrclass).build()

        return None 

    def file(self): 
        """ Name of ASCII file of Data/Random catalogy
        """
        if self.corr_str is None: 
            warnings.warn("Correction specifier is None so the default will be calculated")

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
