from __future__ import print_function 
import scipy.sparse as sps
import numpy as np
# from dolfin import dx, grad, inner



def linearzd_quadterm(H, linv, retparts=False, hlstr=None):
    """ compute the matrices `L1`, `L2` that represent the linearized convection

    `H(v, v) ~ L1*v + L2*v - H(linv, linv)`


    Parameters:
    ---
    H : (nv, nv*nv) sparse array
        the tensor (as a matrix) that evaluates the convection term
    linv : (nv, 1) numpy array
        the stat at which the linearization is about
    retparts: Boolean, optional
        whether to return the `L1` or `L2` separately, \
        defaults to `False`, i.e. `L1+L2` is returned
    hlstr : str, optional
        name of location from where to load or where to store the \
        the wanted data, if `None` nothing is loaded or stored, \
        defaults to `None`

    """
    try:
        import dolfin_navier_scipy.data_output_utils as dou
        if hlstr is None:
            raise IOError()
        if retparts:
            H1L = dou.load_spa(hlstr + '_H1L.mtx')
            H2L = dou.load_spa(hlstr + '_H2L.mtx')
        else:
            HL = dou.load_spa(hlstr + '.mtx')
        print('loaded `hlmat`')
    except (IOError, ImportError) as e:
        print('assembling hlmat ...')
    nv = linv.size
    if retparts:
        try:
            H1L = H * (sps.kron(sps.eye(nv), linv))
            H2L = H * (sps.kron(linv, sps.eye(nv)))
        except TypeError:  # for earlier scipys
            H1L = H * (sps.kron(sps.eye(nv, nv), linv))
            H2L = H * (sps.kron(linv, sps.eye(nv, nv)))
        return H1L, H2L
    else:
        try:
            HL = H * (sps.kron(sps.eye(nv), linv) + sps.kron(linv, sps.eye(nv)))
        except TypeError:  # for earlier scipys
            HL = H * (sps.kron(sps.eye(nv, nv), linv) + \
                sps.kron(linv, sps.eye(nv, nv)))
        return HL


def eva_quadterm(H, v):
    ''' function to evaluate `H*kron(v, v)` without forming `kron(v, v)`

    Parameters:
    ---
    H : (nv, nv*nv) sparse array
        the tensor (as a matrix) that evaluates the convection term

    '''

    NV = v.size
    hvv = np.zeros((NV, 1))
    for k, vi in enumerate(v):
        hviv = H[:, k*NV:(k+1)*NV]*(vi[0]*v)
        hvv = hvv + hviv
    return np.array(hvv)
