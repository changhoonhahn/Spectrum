'''
'''
import numpy as np 

import env
from fort import Fcode

def compile(catalog): 
    catinfo = {'catalog': {'name': catalog, 'n_mock': 1}, 'spec': {'P0': 20000, 'Lbox': 3600, 'Ngrid':360, 'ell': 0}} 
    fort = Fcode('fft', catinfo) 
    print fort.code
    print fort.exe
    fort.compile()
    return None 


if __name__=='__main__': 
    compile('nseries') 
