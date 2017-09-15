'''
'''
import numpy as np

import env 
import util as UT
from spec import Spec


def buildPk(catalog): 
    catinfo = {'catalog': {'name': catalog, 'n_mock': 1}, 
            'spec': {'type': 'pk', 'P0': 20000, 'Lbox': 3600, 'Ngrid':480, 'ell': 2}} 
    pkay = Spec(catinfo) 
    print pkay.file_name 
    pkay.build()
    return None 


if __name__=="__main__": 
    buildPk('qpm') 
