"""

Plot power/bispectrum

Author(s): ChangHoon Hahn

"""
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# --- Local --- 
from defutility.plotting import prettyplot 
from defutility.plotting import prettycolors 
from spec import Spec
from average import AvgSpec


def plot_pk_comp(cat_corrs, n_mock, ell=0, type='ratio', **kwargs): 
    ''' Plot comparison of average power spectrum monopole or quadrupole (avg(P(k))) 
    for multiple a list of catalog and correction specifications. Main use is to 
    compare the effects of fiber collisions correction method. However, it can be
    used to compare any power spectra as long as cat_corr dictionary is specified.

    --------------------------------------------------------------------------
    Paramters
    --------------------------------------------------------------------------
    cat_corrs : list of catalog correction dictionary 
    n_mock : number of mocks 
    ell : ell-th component of multipole decomposition 
    type : 'regular' compares the actual P2(k) values, 'ratio' compares the ratio 
    with the true P2(k), 'residual' compares the difference with the true P2(k) 
    
    --------------------------------------------------------------------------
    Notes
    --------------------------------------------------------------------------
    * Long ass code with a lot of idiosyncracies.
    * Make sure k values agree with each other. 
    
    --------------------------------------------------------------------------
    Example
    --------------------------------------------------------------------------
    cat_corrs = [ 
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'true'}},
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'upweight'}},
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'dlospeak', 'fit': 'gauss', 'sigma': 3.9, 'fpeak': 0.68}} 
            ] 

    plot_pk_comp(cat_corrs, 84, quad=False, type='Pk')
    plot_pk_comp(cat_corrs, 84, quad=False, type='ratio')

    '''

    if 'Ngrid' in kwargs.keys():
        Ngrid = kwargs['Ngrid']
    else: 
        Ngrid = 360 

    if isinstance(n_mock, int): 
        n_mock_list = [ n_mock for i in xrange(len(cat_corrs)) ] 
    else: 
        if len(n_mock) != len(cat_corrs): 
            raise ValueError()
        else: 
            n_mock_list = n_mock

    corr_str = ''

    prettyplot()                         # set up plot 
    pretty_colors = prettycolors()
    
    if 'figsize' in kwargs.keys(): 
        fig = plt.figure(1, kwargs['figsize'])
    else: 
        fig = plt.figure(1, figsize=(7, 8)) # set up figure 
    sub = fig.add_subplot(111)

    for i_corr, cat_corr in enumerate(cat_corrs):

        catdict = cat_corr['catalog']
        corrdict = cat_corr['correction']
        specdict = {
                'P0': 20000,
                'Lbox': 3600, 
                'Ngrid': Ngrid, 
                'ell': ell 
                }
        cat_corr_i = {
                'catalog': {'name': catdict['name'], 'n_mock': 1}, 
                'correction': corrdict, 
                'spec': specdict
                }

        avg_spec = AvgSpec(n_mock_list[i_corr], 'pk', cat_corr_i)
        avg_spec.read()
    
        spec_type = 'pk'
        spec_key = ''.join(['p', str(ell), 'k'])

        k_arr = avg_spec.k
        avg_pk = getattr(avg_spec, spec_key)
        
        if type == 'Pk':         # Compare P(k) to each other 

            sub.plot( 
                    k_arr, avg_pk, 
                    color = pretty_colors[i_corr + 1], 
                    label = plot_label(cat_corr),
                    lw = 4
                    ) 

        elif type == 'Pk_err':  # Compare P(k) with sample variance error bar

            pk_err = avg_spec.stddev()

            sub.errorbar( 
                    k_arr, avg_pk, 
                    yerr = [pk_err, pk_err], 
                    color = pretty_colors[i_corr + 1], 
                    label = plot_label(cat_corr),
                    fmt='--o'
                    ) 

        elif type == 'Pk_all': 
            if isinstance(n_mock_list[i_corr], int): 
                n_mock_list_i = range(1, n_mock_list[i_corr]+1)
            else: 
                n_mock_list_i = n_mock_list[i_corr]

            for i_mock in n_mock_list_i:
                k_i, spec_i_spec = avg_spec.spec_i(i_mock)

                sub.plot( 
                        k_i, spec_i_spec, 
                        color = '0.25', 
                        lw = 1 
                        ) 
            sub.plot( 
                    k_arr, avg_pk, 
                    color = pretty_colors[i_corr + 1], 
                    label = plot_label(cat_corr),
                    lw = 4
                    ) 
        
        elif type == 'kPk':         # Compare k^1.5 * P(k) with each other. Enhances the 
            # BAO signature? (Upon Paco's request). 

            kPk = k_arr**1.5 * avg_pk
                
            sub.scatter(
                    k_arr, kPk, 
                    color = pretty_colors[i_corr+1], 
                    label = plot_label(cat_corr)
                    )

        elif type == 'ratio':       # Compare the ratio of the power spectra (P/P_denom)

            if i_corr == 0 :        
                avg_pk_denom = avg_pk
                #denom_cat = catdict['name']
                denom_cat = corrdict['name']

            else: 
                sub.scatter(
                        k_arr, 
                        avg_pk/avg_pk_denom, 
                        color = pretty_colors[i_corr+1], 
                        label = plot_label(cat_corr)
                        )
                
                print plot_label(cat_corr)

                largescale = np.where(k_arr < 0.2)
                smallscale = np.where(k_arr > 0.2)
                print np.sum( np.abs((avg_pk/avg_pk_denom) - 1.0 ) )
                print 'Large scale k < 0.2'
                print np.sum( np.abs((avg_pk[largescale]/avg_pk_denom[largescale]) - 1.0 ) )
                print 'Small scale k > 0.2'
                print np.sum( np.abs((avg_pk[smallscale]/avg_pk_denom[smallscale]) - 1.0 ) )

            if corrdict['name'] == 'true': 

                pk_err = avg_spec.stddev()
        
                sub.plot( 
                        k_arr, 1.0 + pk_err/np.abs(avg_pk),  
                        color = 'k', 
                        lw = 2, 
                        ls = '-.', 
                        label = r"$\mathtt{1 + \Delta P^{true} (k) / P^{true}}$"
                        ) 
        
                sub.plot( 
                        k_arr, 1.0 + -1.0 * pk_err/np.abs(avg_pk),  
                        color = 'k',  
                        lw = 2, 
                        ls = '-.'
                        ) 

        elif type == 'l1_norm':

            if i_corr == 0 :        
                avg_pk_denom = avg_pk
                denom_cat = catdict['name']

            else: 
                sub.scatter(
                        k_arr, 
                        avg_pk - avg_pk_denom, 
                        color = pretty_colors[i_corr+1], 
                        label = plot_label(cat_corr)
                        )
                
                print plot_label(cat_corr)
                print (avg_pk-avg_pk_denom)[-10:]
        
        del avg_pk
        
        # Specify corrections for figure file name  
        if 'dlospeak' in corrdict['name']: 
            try:
                corr_str += ''.join([
                    catdict['name'], '_', 
                    corrdict['name'], '_', 
                    corrdict['fit'], '_',
                    '_sigma', str(corrdict['sigma']), 
                    'fpeak', str(corrdict['fpeak'])
                    ]) 
            except KeyError: 
                corr_str += ''.join([
                    catdict['name'], '_', 
                    corrdict['name'], '_', 
                    '_sigma', str(corrdict['sigma'])
                    ]) 
        elif corrdict['name'] == 'fourier_tophat': 
            corr_str += ''.join([
                catdict['name'], '_', 
                corrdict['name'],  
                '.fs', str(round(corrdict['fs'], 1)), 
                '.rc', str(round(corrdict['rc'], 2)), 
                '.kfit', str(round(corrdict['k_fit'], 2)), 
                '.kfixed', str(round(corrdict['k_fixed'], 2))
                ])
        else: 
            corr_str += ''.join([ 
                catdict['name'], '_', 
                corrdict['name']
                ]) 
    
    # Dictate the x-range and y-range of the plotting
    # based on type of comparison 
    if (type == 'Pk') or (type == 'Pk_err') or (type == 'Pk_all'):
        if 'yrange' in kwargs.keys(): 
            ylimit = kwargs['yrange'] 
            yytext = 10**.5*min(ylimit) 
        else: 
            ylimit = [10**2,10**5.5]
            yytext = 10**2.5
        
        if 'ylabel' in kwargs.keys(): 
            ylabel = kwargs['ylabel']
        else: 
            ylabel = r'$\mathtt{P_'+str(ell)+'(k)}$'

        if 'xscale' in kwargs.keys(): 
            sub.set_xscale(kwargs['xscale']) 
        else: 
            sub.set_xscale('log')

        if 'yscale' in kwargs.keys(): 
            sub.set_yscale(kwargs['yscale'])
        else: 
            sub.set_yscale('log')
        
        if type == 'Pk': 
            resid_str = ''
        elif type == 'Pk_err': 
            resid_str = '_err'
        elif type == 'Pk_all':
            resid_str = '_all'
    
    elif type == 'kPk':
        if 'yrange' in kwargs.keys(): 
            ylimit = kwargs['yrange'] 
            yytext = 10**.5*min(ylimit) 
        else: 
            ylimit = [10**0,10**2.0]
            yytext = 10**0.1

        if 'ylabel' in kwargs.keys(): 
            ylabel = kwargs['ylabel']
        else: 
            ylabel = r'$\mathtt{k^{1.5} P_0(k)}$'

        if 'xscale' in kwargs.keys(): 
            sub.set_xscale(kwargs['xscale']) 
        else: 
            sub.set_xscale('log')

        if 'yscale' in kwargs.keys(): 
            sub.set_yscale(kwargs['yscale'])
        else: 
            sub.set_yscale('log')

        resid_str = '_kPk'

    elif type == 'ratio': 
        if 'yrange' in kwargs.keys(): 
            ylimit = kwargs['yrange'] 
            yytext = 0.05 + min(ylimit) 
        else: 
            ylimit = [0.5, 1.5] 
            yytext = 0.55

        ylabel = ''.join([
            r"$\mathtt{\overline{P_", str(ell), "(k)}/\overline{P_", str(ell), r"(k)_{\rm{", denom_cat, "}}}}$"
            ])
        
        if 'xscale' in kwargs.keys(): 
            sub.set_xscale(kwargs['xscale']) 
        else: 
            sub.set_xscale('log')

        if 'yscale' in kwargs.keys(): 
            sub.set_yscale(kwargs['yscale'])
        else: 
            pass

        resid_str = '_ratio'

    elif type == 'residual':
        if 'yrange' in kwargs.keys(): 
            ylimit = kwargs['yrange'] 
            yytext = np.mean(ylimit) 
        else: 
            ylimit = [0.0, 5.0]
            yytext = 2.5

        ylabel = ''.join([
            r"$\mathtt{|\overline{P_", str(ell), "(k)} - \overline{P_", str(ell), r"(k)_{\rm{True}}}|/\Delta P_", str(ell), "}$"
            ])

        if 'xscale' in kwargs.keys(): 
            sub.set_xscale(kwargs['xscale']) 
        else: 
            sub.set_xscale('log')

        if 'yscale' in kwargs.keys(): 
            sub.set_yscale(kwargs['yscale'])
            if (kwargs['yscale'] == 'log') and \
                    ('yrange' not in kwargs.keys()):
                ylimit = [10**-3, 10**1] 
        else: 
            pass

        resid_str = '_residual'

    elif type == 'l1_norm': 

        if 'yrange' in kwargs.keys(): 
            ylimit = kwargs['yrange'] 
            yytext = 0.05 + min(ylimit) 
        else: 
            ylimit = None 
            yytext = 0.55

        ylabel = ''.join([
            r"$\mathtt{\overline{P_", str(ell), "(k)} - \overline{P_", str(ell), r"(k)_{\rm{", denom_cat, "}}}}$"
            ])
        
        if 'xscale' in kwargs.keys(): 
            sub.set_xscale(kwargs['xscale']) 
        else: 
            sub.set_xscale('log')

        if 'yscale' in kwargs.keys(): 
            sub.set_yscale(kwargs['yscale'])
        else: 
            pass

        resid_str = '_l1norm'
    else: 
        raise NotImplementedError('asdfasdfasdf') 

    if 'xrange' in kwargs.keys():   # specify x-range 
        sub.set_xlim(kwargs['xrange']) 
        yxtext = 1.5*min(kwargs['xrange'])
    else: 
        sub.set_xlim([10**-3,10**0])
        yxtext = 1.5*10**-3

    if type == 'ratio': 
        sub.axhline(y = 1.0, lw=2, ls='--', c='k')

    sub.set_ylim(ylimit)
    sub.set_xlabel('k (h/Mpc)', fontsize=20)
    sub.set_ylabel(ylabel, fontsize=20)
    
    # Display the number of mocks for given catalog so that
    # I know how many mocks the P(k) is averaged over.
    n_mock_text = '\n'.join([
                ' '.join([
                    str(n_mock_list[ii]), 
                    ((cat_corrs[ii])['catalog'])['name'].upper()
                    ]) 
                for ii in xrange(len(n_mock_list))
                ])
    sub.text(yxtext, yytext, n_mock_text)

    sub.legend(scatterpoints=1, loc='upper left', prop={'size':14})
    
    try: 
        n_mock_str = '_'.join([str(nm) for nm in n_mock]) 
    except TypeError:
        n_mock_str = str(n_mock) 

    fig_name = ''.join([
        spec_key, '_', 
        n_mock_str, 
        'mock_fibcoll_', 
        corr_str, 
        resid_str, 
        '_comparison_Ngrid', 
        str(Ngrid), 
        '.png'
        ])     

    fig_dir = '/home/users/hahn/powercode/FiberCollisions/figure/'
    print fig_name
    fig.savefig(
            ''.join([fig_dir, fig_name]), 
            bbox_inches="tight"
            )
    #plt.show()
    plt.close()

    return None 

def plot_delpoverp_comp(cat_corrs, n_mock, ell=0, **kwargs): 
    ''' 
    Plot comparison of delta P/ avg(P) for multiple lists of catalog and correction 
    specifications. 

    --------------------------------------------------------------------------
    Paramters
    --------------------------------------------------------------------------
    cat_corrs : list of catalog correction dictionary 
    n_mock : number of mocks 
    ell : ell-th component of multipole decomposition 
    
    --------------------------------------------------------------------------
    Notes
    --------------------------------------------------------------------------
    * Long ass code with a lot of idiosyncracies.
    * Make sure k values agree with each other. 
    
    --------------------------------------------------------------------------
    Example
    --------------------------------------------------------------------------
    cat_corrs = [ 
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'true'}},
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'upweight'}},
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'dlospeak', 'fit': 'gauss', 'sigma': 3.9, 'fpeak': 0.68}} 
            ] 

    '''

    if 'Ngrid' in kwargs.keys():
        Ngrid = kwargs['Ngrid']
    else: 
        Ngrid = 360 

    if isinstance(n_mock, int): 
        n_mock_list = [ n_mock for i in xrange(len(cat_corrs)) ] 
    else: 
        if len(n_mock) != len(cat_corrs): 
            raise ValueError()
        else: 
            n_mock_list = n_mock

    corr_str = ''

    prettyplot()                         # set up plot 
    pretty_colors = prettycolors()
    
    fig = plt.figure(1, figsize=(14, 8)) # set up figure 
    sub = fig.add_subplot(111)

    for i_corr, cat_corr in enumerate(cat_corrs):

        catdict = cat_corr['catalog']
        corrdict = cat_corr['correction']
        specdict = {
                'P0': 20000,
                'Lbox': 3600, 
                'Ngrid': Ngrid, 
                'ell': ell 
                }
        cat_corr_i = {
                'catalog': {'name': catdict['name'], 'n_mock': 1}, 
                'correction': corrdict, 
                'spec': specdict
                }

        avg_spec = AvgSpec(n_mock_list[i_corr], 'pk', cat_corr_i)
        avg_spec.read()
    
        spec_type = 'pk'
        spec_key = ''.join(['p', str(ell), 'k'])

        k_arr = avg_spec.k
        avg_pk = getattr(avg_spec, spec_key)
        pk_err = avg_spec.stddev()
        
        sub.plot( 
                k_arr, pk_err/np.abs(avg_pk),  
                color = pretty_colors[i_corr + 1], 
                lw = 4
                ) 
        
        # Specify corrections for figure file name  
        if 'dlospeak' in corrdict['name']: 
            try:
                corr_str += ''.join([
                    catdict['name'], '_', 
                    corrdict['name'], '_', 
                    corrdict['fit'], '_',
                    '_sigma', str(corrdict['sigma']), 
                    'fpeak', str(corrdict['fpeak'])
                    ]) 
            except KeyError: 
                corr_str += ''.join([
                    catdict['name'], '_', 
                    corrdict['name'], '_', 
                    '_sigma', str(corrdict['sigma'])
                    ]) 
        else: 
            corr_str += ''.join([ 
                catdict['name'], '_', 
                corrdict['name']
                ]) 
    
    # x-axis
    if 'xscale' in kwargs.keys(): 
        sub.set_xscale(kwargs['xscale']) 
    else: 
        sub.set_xscale('log')

    if 'xrange' in kwargs.keys():   # specify x-range 
        sub.set_xlim(kwargs['xrange']) 
        yxtext = 1.5*min(kwargs['xrange'])
    else: 
        sub.set_xlim([10**-3,10**0])
        yxtext = 1.5*10**-3
    sub.set_xlabel('k (h/Mpc)', fontsize=20)

    # y-axis
    if 'ylabel' in kwargs.keys(): 
        ylabel = kwargs['ylabel']
    else: 
        ylabel = r'$\mathtt{\Delta P_'+str(ell)+'(k)/|\overline{P_'+str(ell)+'}|}$'
    if 'yrange' in kwargs.keys():   # specify x-range 
        sub.set_ylim(kwargs['yrange']) 
    else: 
        sub.set_ylim([-1.,1.0])
    sub.set_xlabel('k (h/Mpc)', fontsize=20)
    sub.set_ylabel(ylabel, fontsize=20)
    
    # Display the number of mocks for given catalog so that
    # I know how many mocks the P(k) is averaged over.
    n_mock_text = '\n'.join([
                ' '.join([
                    str(n_mock_list[ii]), 
                    ((cat_corrs[ii])['catalog'])['name'].upper()
                    ]) 
                for ii in xrange(len(n_mock_list))
                ])
    sub.text(yxtext, 0.05, n_mock_text)

    sub.legend(scatterpoints=1, loc='upper left', prop={'size':14})
    
    try: 
        n_mock_str = '_'.join([str(nm) for nm in n_mock]) 
    except TypeError:
        n_mock_str = str(n_mock) 

    plt.figtext(0, 0, "testing \n testing")

    fig_name = ''.join([
        'del', spec_key, 'over', spec_key, '_', 
        n_mock_str, 
        'mock_fibcoll_', 
        corr_str, 
        '_comparison_Ngrid', 
        str(Ngrid), 
        '.png'
        ])     

    fig_dir = '/home/users/hahn/powercode/FiberCollisions/figure/'
    fig.savefig(
            ''.join([fig_dir, fig_name]), 
            bbox_inches="tight"
            )

    #plt.show()
    plt.close()

    return None 

def plot_bk_comp(cat_corrs, n_mock, k='max', type='bk_ratio', **kwargs): 
    ''' Plot comparison of average bispectrum (B(k) or Q(k)) for multiple a list 
    of catalog and correction specifications. Main use is to compare the effects 
    of fiber collisions correction method. However, it can be used to compare any 
    power spectra as long as cat_corr dictionary is specified.

    --------------------------------------------------------------------------
    Paramters
    --------------------------------------------------------------------------
    cat_corrs : list of catalog correction dictionary 
    n_mock : number of mocks 
    type : 'regular' compares the actual P2(k) values, 'ratio' compares the ratio 
    with the true P2(k), 'residual' compares the difference with the true P2(k) 
    
    --------------------------------------------------------------------------
    Notes
    --------------------------------------------------------------------------
    * Long ass code with a lot of idiosyncracies.
    * Make sure k values agree with each other. 
    
    --------------------------------------------------------------------------
    Example
    --------------------------------------------------------------------------
    cat_corrs = [ 
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'true'}},
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'upweight'}},
            {'catalog': {'name': 'nseries'}, 'correction': {'name': 'dlospeak', 'fit': 'gauss', 'sigma': 3.9, 'fpeak': 0.68}} 
            ] 

    plot_bk_comp(cat_corrs, 84, quad=False, type='Pk')
    plot_bk_comp(cat_corrs, 84, quad=False, type='ratio')

    '''

    if 'Ngrid' in kwargs.keys():
        Ngrid = kwargs['Ngrid']
    else: 
        Ngrid = 360 

    if isinstance(n_mock, int): 
        n_mock_list = [ n_mock for i in xrange(len(cat_corrs)) ] 
    else: 
        if len(n_mock) != len(cat_corrs): 
            raise ValueError()
        else: 
            n_mock_list = n_mock

    corr_str = ''

    prettyplot()                         # set up plot 
    pretty_colors = prettycolors()
    
    if 'figsize' in kwargs.keys(): 
        fig = plt.figure(1, kwargs['figsize'])
    else: 
        fig = plt.figure(1, figsize=(7, 8)) # set up figure 
    sub = fig.add_subplot(111)

    for i_corr, cat_corr in enumerate(cat_corrs):

        catdict = cat_corr['catalog']
        try: 
            corrdict = cat_corr['correction']
        except KeyError: 
            corrdict = {} 
        specdict = {'P0': 20000, 'Lbox': 3600, 'Ngrid': Ngrid}
        cat_corr_i = {'catalog': {'name': catdict['name'], 'n_mock': 1}, 'correction': corrdict, 'spec': specdict}

        avg_spec = AvgSpec(n_mock_list[i_corr], 'bk', cat_corr_i)
        avg_spec.read()
    
        spec_type = 'bk'
        if 'bk' in type.lower(): 
            spec_key = 'avg_bk'
        elif 'qk' in type.lower():
            spec_key = 'avg_qk'
        else: 
            raise ValueError('Bk or Qk has to be specified in type')

        k1_arr = avg_spec.k1
        k2_arr = avg_spec.k2
        k3_arr = avg_spec.k3
        avg_bisp_k = getattr(avg_spec, spec_key)

        if k == 'max': 
            k_arr = np.max([k1_arr, k2_arr, k3_arr], axis=0)
        elif k == 'avg': 
            k_arr = np.mean([k1_arr, k2_arr, k3_arr], axis=0)
        elif k == 'i_triangle': 
            k_arr = np.arange(len(k1_arr))
        
        if 'err' in type:       # Compare B/Q(k) with sample variance error bar

            bk_err, qk_err = avg_spec.stddev()
            if 'bk' in type.lower(): 
                spec_err = bk_err
            elif 'qk' in type.lower(): 
                spec_err = qk_err

            sub.errorbar( 
                    k_arr, 
                    avg_bisp_k, 
                    yerr = [spec_err, spec_err], 
                    color = pretty_colors[i_corr + 1], 
                    label = plot_label(cat_corr),
                    fmt='--o'
                    ) 
        elif type == 'ratio':       # Compare the ratio of the power spectra (P/P_denom)

            if i_corr == 0 :        
                avg_bisp_denom = avg_bisp_k

                denom_cat = catdict['name']
                denom_corr = corrdict['name']

            else: 
                sub.scatter(
                        k_arr, 
                        avg_bisp_k/avg_bisp_denom, 
                        color = pretty_colors[i_corr+1], 
                        label = plot_label(cat_corr)
                        )
                
                print plot_label(cat_corr)

            if corrdict['name'] == 'true': 
            
                bk_err, qk_err = avg_spec.stddev()
                if 'bk' in type.lower(): 
                    spec_err = bk_err
                elif 'qk' in type.lower(): 
                    spec_err = qk_err

                sub.plot( 
                        k_arr, 1.0 + spec_err/np.abs(avg_bisp_k),  
                        color = 'k', 
                        lw = 2, 
                        ls = '-.', 
                        label = r"$\mathtt{1 + \Delta B^{true}/ B^{true}}$"
                        ) 
        
                sub.plot( 
                        k_arr, 1.0 - spec_err/np.abs(avg_bisp_k),  
                        color = 'k',  
                        lw = 2, 
                        ls = '-.'
                        ) 
        
        else:
            sub.scatter( 
                    k_arr, avg_bisp_k, 
                    color = pretty_colors[i_corr + 1], 
                    lw = 4
                    ) 

    # Dictate the x-range and y-range of the plotting
    # based on type of comparison 
    if 'bk' in type.lower(): 
        if 'ratio' not in type.lower(): 
            ylimit = [10**7,5*10**11]
        else: 
            ylimit = [0.5, 2.0]
        
        if 'ylabel' in kwargs.keys(): 
            ylabel = kwargs['ylabel']
        else: 
            if 'ratio' not in type.lower(): 
                ylabel = r'$\mathtt{B_(k_{'+k+'})}$'
            else: 
                ylabel = ''.join([
                    r"$\mathtt{\overline{B(k_{", k, "})}", 
                    "/", 
                    "\overline{B(k_{",k,r"})_{\rm{", denom_cat, 
                    "}}}}$"
                ])

    elif 'qk' in type.lower(): 
        if 'ratio' not in type.lower(): 
            ylimit = [10**7,5*10**11]
        else: 
            ylimit = [0.5, 2.0]

        if 'ylabel' in kwargs.keys(): 
            ylabel = kwargs['ylabel']
        else: 
            if 'ratio' not in type.lower(): 
                ylabel = r'$\mathtt{Q_(k_{'+k+'})}$'
            else: 
                ylabel = ''.join([
                    r"$\mathtt{\overline{Q(k_{", k, "})}", 
                    "/", 
                    "\overline{Q(k_{",k,r"})_{\rm{", denom_cat, 
                    "}}}}$"
                    ])
    
    if 'yrange' in kwargs.keys(): 
        ylimit = kwargs['yrange'] 
    #yytext = 10**.5*min(ylimit) 
        
    if 'xscale' in kwargs.keys(): 
        sub.set_xscale(kwargs['xscale']) 
    else: 
        if k != 'i_triangle': 
            sub.set_xscale('log')
        
    if 'yscale' in kwargs.keys(): 
        sub.set_yscale(kwargs['yscale'])
    else: 
        sub.set_yscale('log')
        
    resid_str = ''
    if 'err' in type: 
        resid_str = '_err'
    elif 'ratio' in type: 
        resid_str = '_ratio'

    if 'xrange' in kwargs.keys():   # specify x-range 
        sub.set_xlim(kwargs['xrange']) 
        #yxtext = 1.5*min(kwargs['xrange'])
    else: 
        if k == 'i_triangle':
            sub.set_xlim([0, np.max(k_arr)])
        else:
            sub.set_xlim([np.min(k_arr), 2.0*np.max(k_arr)])
        #yxtext = 1.5*10**-3

    if 'ratio' in type: 
        sub.axhline(y = 1.0, lw=2, ls='--', c='k')

    sub.set_ylim(ylimit)
    if k == 'i_triangle': 
        sub.set_xlabel('Triangles', fontsize=20)
    else:
        sub.set_xlabel(r'$k_{'+k+'}$ (h/Mpc)', fontsize=20)
    sub.set_ylabel(ylabel, fontsize=20)
    
    try: 
        n_mock_str = '_'.join([str(nm) for nm in n_mock]) 
    except TypeError:
        n_mock_str = str(n_mock) 

    fig_name = ''.join([
        spec_key, '_', 
        n_mock_str, 
        'mock_fibcoll_', 
        corr_str, 
        resid_str, 
        '_comparison_Ngrid', 
        str(Ngrid), 
        '.png'
        ])     

    fig_dir = '/home/users/hahn/powercode/FiberCollisions/figure/'
    print fig_name
    #fig.savefig(
    #        ''.join([fig_dir, fig_name]), 
    #        bbox_inches="tight"
    #        )
    plt.show()
    plt.close()

    return None 

def plot_label(cat_corr): 
    """ Labels for plotting that specify the catalog and 
    correction names. Written for plot_pk_comp but can be used
    elsewhere.
    """
    catdict = cat_corr['catalog'] 
    corrdict = cat_corr['correction']

    if 'dlospeak' in corrdict['name']: 
        
        try:
            label = ' '.join([
                catdict['name'].upper(), 
                ''.join([
                    corrdict['name'].upper(), ':', 
                    str(corrdict['sigma']), ',', str(corrdict['fpeak'])
                    ])
                ])
        except KeyError: 
            label = ' '.join([
                catdict['name'].upper(), 
                ''.join([
                    corrdict['name'].upper(), ':', 
                    str(corrdict['sigma'])
                    ])
                ])

    elif corrdict['name'] in ('true', 'upweight'): 

        label = ' '.join([ 
            catdict['name'].upper(), 
            corrdict['name'].upper()
            ])

    elif corrdict['name'] == 'fourier_tophat':

        label = ''.join([
            'FOURIER TOPHAT:', 
            '$f_s = '+str(round(corrdict['fs'], 1))+'$,', 
            '$r_c = '+str(round(corrdict['rc'], 2))+'$', '\n',
            '$k_{fit} = '+str(round(corrdict['k_fit'], 2))+'$,', 
            '$k_{fixed} = '+str(round(corrdict['k_fixed'], 2))+'$'
            ])
        print label

    else: 
        raise NotImplementedError()

    return label 

if __name__=='__main__': 

    cat_corrs = [ 
            {
                'catalog': {'name': 'nseries'}, 
                'correction': {'name': 'upweight'}
                },
            {
                'catalog': {'name': 'nseries'}, 
                'correction': {'name': 'true'}
                },
            {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'fourier_tophat', 'fs': 1.0, 'rc': 0.43, 'k_fit': 0.7, 'k_fixed': 0.84}
            }
            ]
    #plot_pk_comp(cat_corrs, 20, Ngrid=960, ell=2, type='Pk')#, yrange=[0.0, 2.0], xrange=[10**-1, 10**0.])
    #plot_pk_comp(cat_corrs, 20, Ngrid=960, ell=2, type='ratio', yrange=[0.0, 2.0], figsize=[14,8])#, xrange=[10**-1, 10**0.])
    cat_corrs = [{'catalog': {'name': 'nseries'}}]
    plot_bk_comp(cat_corrs, 1, Ngrid=360, k='max', type='bk', figsize=[14,8])
    plot_bk_comp(cat_corrs, 1, Ngrid=360, k='avg', type='bk', figsize=[14,8])
    plot_bk_comp(cat_corrs, 1, Ngrid=360, k='i_triangle', type='bk', figsize=[14,8])
