"""

Manage directories for Fibercollisions Project

"""

def direc(type, cat_corr, **kwargs): 
    """ Get data directories
    """

    if type not in ('data', 'fft', 'spec'): 
        raise ValueError()

    catdict = (cat_corr)['catalog']
    catname = catdict['name']

    localdir = '/mount/riachuelo1/hahn/'    # riachuelo local directory

    if type == 'fft': 
        localdir = '/mount/chichipio2/hahn/'

    if type == 'data':
        typedir = 'data/'
    elif type == 'fft': 
        typedir = 'FFT/'
    elif type == 'spec': 
        typedir = 'power/'
    
    if catname == 'lasdamasgeo': 
        catdir = 'LasDamas/Geo/'
    elif catname == 'tilingmock': 
        catdir = 'tiling_mocks/'
    elif catname == 'qpm': 
        catdir = 'QPM/dr12d/'
    elif catname == 'nseries': 
        catdir = 'Nseries/'
    elif catname == 'patchy': 
        catdir = 'PATCHY/dr12/v6c/'
    elif 'bigmd' in catname:
        catdir = 'BigMD/'
    elif 'cmass' in catname: 
        if catname == 'cmass': 
            catdir = 'CMASS/'
        elif 'cmasslowz' in catname: 
            catdir = 'CMASS/dr12v5/'
    else: 
        raise NotImplementedError()

    dir = ''.join([
        localdir, 
        typedir, 
        catdir])

    return dir 
