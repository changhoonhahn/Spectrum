'''
'''
import numpy as np 
# -- local -- 
import fft as Fft


def buildfft(catalog): 
    catinfo = {'type': 'random', 
            'catalog': {'name': catalog, 'n_mock': 1}, 
            'spec': {'P0': 20000, 'Lbox': 3600, 'Ngrid':360, 'ell': 0}} 

    efft = Fft.Fft(catinfo) 
    print efft.data_file
    print efft.file_name 
    efft.build()
    return None 

if __name__=='__main__': 
    buildfft('nseries') 
