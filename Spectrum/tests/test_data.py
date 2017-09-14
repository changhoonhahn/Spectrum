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
    D = Data.Data(catinfo, cosmo='fid') 
    print D.file_name  
    D.read() 
    print D.ra[:10]



if __name__=='__main__': 
    read_data('nseries')
