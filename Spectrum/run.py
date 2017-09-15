'''
'''
import os
import sys 
from spec import Spec

from interruptible_pool import InterruptiblePool as Pewl


def buildPk(params): 
    ''' calculate p(k) using pythons multiprocessing  
    '''
    catalog, n_mock = params
    catinfo = {'catalog': {'name': catalog, 'n_mock': n_mock}, 
            'spec': {'type': 'pk', 'P0': 20000, 'Lbox': 3600, 'Ngrid':480, 'ell': 2}} 

    pkay = Spec(catinfo) 
    if not os.path.isfile(pkay.file_name): 
        print pkay.file_name 
        pkay.build()
    return None 


if __name__=="__main__": 
    Nthreads = int(sys.argv[1])
    print 'running on ', Nthreads, ' threads'
    pool = Pewl(processes=Nthreads)
    mapfn = pool.map

    nmock0 = int(sys.argv[2])
    nmock1 = int(sys.argv[3])
    arglist = [['qpm', i_mock] for i_mock in range(nmock0, nmock1+1)]

    mapfn(buildPk, [arg for arg in arglist])
    pool.close()
    pool.terminate()
    pool.join() 
