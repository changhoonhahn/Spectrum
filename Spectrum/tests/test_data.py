'''

test data.py

'''
import numpy as np 
# -- local -- 
import env 
import data as Data 


def read_data(catalog): 
    ''' Test data.py 
    '''
    catinfo = {'catalog': {'name': catalog, 'n_mock': 1}} 
    D = Data.Data(catinfo) 
    print D.file_name  
    D.read() 
    print D.ra[:10]
    return None 


def read_random(catalog): 
    catinfo = {'catalog': {'type': 'random', 'name': catalog, 'n_mock': 1}} 
    R = Data.Data(catinfo) 
    print R.file_name  
    R.read() 
    return None


if __name__=='__main__': 
    read_random('nseries')
