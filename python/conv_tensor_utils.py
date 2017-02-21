import scipy.sparse as sps
import numpy as np
# from dolfin import dx, grad, inner

import dolfin_navier_scipy.data_output_utils as dou


def linearzd_quadterm(H, linv, retparts=False, hlstr=None):
    """ compute the matrices that represent the linearized convection


    Parameters:
    ---
    hlstr : str, optional
        name of location from where to load or where to store the \
        the wanted data, if `None` nothing is loaded or stored, \
        defaults to `None

    """
    try:
        if hlstr is None:
            raise IOError()
        if retparts:
            H1L = dou.load_spa(hlstr + '_H1L.mtx')
            H2L = dou.load_spa(hlstr + '_H2L.mtx')
        else:
            HL = dou.load_spa(hlstr + '.mtx')
        print 'loaded `hlmat`'
    except IOError:
        print 'assembling hlmat ...'
    nv = linv.size
    if retparts:
        H1L = H * (sps.kron(sps.eye(nv), linv))
        H2L = H * (sps.kron(linv, sps.eye(nv)))
        return H1L, H2L
    else:
        HL = H * (sps.kron(sps.eye(nv), linv) + sps.kron(linv, sps.eye(nv)))
        return HL


def eva_quadterm(H, v):
    ''' function to evaluate `H*kron(v, v)` without forming `kron(v, v)`
    '''
    NV = v.size
    hvv = np.zeros((NV, 1))
    # with dou.Timer('eva: h*vv'):
    for k, vi in enumerate(v):
        hviv = H[:, k*NV:(k+1)*NV]*(vi[0]*v)
        hvv = hvv + hviv
    return np.array(hvv)
