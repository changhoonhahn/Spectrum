'''
'''
import numpy as np 
# -- local -- 
import fft as Fft

def fft(catalog): 
    catinfo = {'type': 'data', 'catalog': {'name': catalog, 'n_mock': 1}} 

    efft = Fft.Fft(catinfo) 
    print efft.data_file
    print efft.file_name 
   

if __name__=='__main__': 
    fft('nseries') 
