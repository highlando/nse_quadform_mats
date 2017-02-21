import dolfin
import os
import numpy as np
import datetime
import scipy.io
import scipy.sparse.linalg as spsla

import dolfin_navier_scipy.dolfin_to_sparrays as dts
import dolfin_navier_scipy.data_output_utils as dou
import dolfin_navier_scipy.stokes_navier_utils as snu
import dolfin_navier_scipy.problem_setups as dnsps

import distr_control_fenics.cont_obs_utils as cou

dolfin.parameters.linear_algebra_backend = 'Eigen'


def comp_exp_nsbili(problemname='drivencavity', N=10, bccontrol=False,
                    mddir='pathtodatastorage', Re=1, palpha=1, visu=False):

    # the Reynoldsnumber will be multiplied to the matrices later

    femp, stokesmatsc, rhsd_vfrc, rhsd_stbc = \
        dnsps.get_sysmats(problem=problemname, N=N, Re=Re,
                          bccontrol=bccontrol, ppin=None)

    invinds = femp['invinds']
    A, J, M = stokesmatsc['A'], stokesmatsc['J'], stokesmatsc['M']
    fv, fp = rhsd_vfrc['fvc'], rhsd_vfrc['fpr']
    fv_diff, fp_div = rhsd_stbc['fv'], rhsd_stbc['fp']
    invinds = femp['invinds']
    NV = invinds.shape[0]

    # #######################
    soldict = stokesmatsc  # containing A, J, JT
    soldict.update(femp)  # adding V, Q, invinds, diribcs
    soldict.update(rhsd_vfrc)  # adding fvc, fpr

    soldict.update(fv=fv+fv_diff, fp=fp+fp_div,
                   nu=femp['nu'])
    if bccontrol:
        Arob = stokesmatsc['Arob']
        Brob = stokesmatsc['Brob']
        soldict.update(A=A+1./palpha*Arob)
        # TODO: general task -TODO- what about the `nu` TODO
        # no addition to rhs for the uncontrolled, i.e. `u=0`, case

    # compute the uncontrolled steady state STOKES solution
    vp_ss_stokes, _ = snu.\
        solve_steadystate_nse(return_vp=True,
                              clearprvdata=True,
                              vel_pcrd_stps=0, vel_nwtn_stps=0,
                              data_prfx='results/', paraviewoutput=True,
                              vfileprfx='results/v_stokes',
                              pfileprfx='results/p_stokes',
                              **soldict)
    v_ss_stokes, p_ss_stokes = vp_ss_stokes[:NV], vp_ss_stokes[NV:]

    # compute the uncontrolled steady state Navier-Stokes solution
    vp_ss_nse, list_norm_nwtnupd = snu.\
        solve_steadystate_nse(return_vp=True,
                              clearprvdata=True,
                              data_prfx='results/', paraviewoutput=True,
                              vfileprfx='results/v', pfileprfx='results/p',
                              **soldict)
    v_ss_nse, p_ss_nse = vp_ss_nse[:NV], vp_ss_nse[NV:]

    manual_paraview_data(vvec=v_ss_nse, V=femp['V'], invinds=invinds,
                         diribcs=femp['diribcs'], fstr='nsesstest')

    # #######################

    data_prfx = problemname  # + '__noRe__'
    NU, NY = 3, 4

    # specify in what spatial direction Bu changes. The remaining is constant
    if problemname == 'drivencavity':
        uspacedep = 0
    elif problemname == 'cylinderwake':
        uspacedep = 1

    # output
    cwd = os.getcwd()
    try:
        os.chdir(mddir)
    except OSError:
        raise Warning('need "' + mddir + '" subdir for storing the data')
    os.chdir(cwd)

    #
    # Control mats
    #
    contsetupstr = data_prfx + '__NV{0}NU{1}NY{2}'.format(NV, NU, NY)
    if bccontrol:
        contsetupstr = contsetupstr + '__bccontrol'

    # get the control and observation operators
    try:
        b_mat = dou.load_spa(mddir + contsetupstr + '__b_mat')
        u_masmat = dou.load_spa(mddir + contsetupstr + '__u_masmat')
        print 'loaded `b_mat`'
    except IOError:
        print 'computing `b_mat`...'
        b_mat, u_masmat = cou.get_inp_opa(cdcoo=femp['cdcoo'], V=femp['V'],
                                          NU=NU, xcomp=uspacedep)
        dou.save_spa(b_mat, mddir + contsetupstr + '__b_mat')
        dou.save_spa(u_masmat, mddir + contsetupstr + '__u_masmat')
    try:
        mc_mat = dou.load_spa(mddir + contsetupstr + '__mc_mat')
        y_masmat = dou.load_spa(mddir + contsetupstr + '__y_masmat')
        print 'loaded `c_mat`'
    except IOError:
        print 'computing `c_mat`...'
        mc_mat, y_masmat = cou.get_mout_opa(odcoo=femp['odcoo'],
                                            V=femp['V'], NY=NY)
        dou.save_spa(mc_mat, mddir + contsetupstr + '__mc_mat')
        dou.save_spa(y_masmat, mddir + contsetupstr + '__y_masmat')

    # restrict the operators to the inner nodes
    mc_mat = mc_mat[:, invinds][:, :]
    b_mat = b_mat[invinds, :][:, :]

    c_mat = apply_massinv(y_masmat, mc_mat, output='sparse')
    # TODO: right choice of norms for y
    #       and necessity of regularization here
    #       by now, we go on number save

    # the pressure observation mean over a small domain
    if problemname == 'cylinderwake':
        podcoo = dict(xmin=0.6,
                      xmax=0.64,
                      ymin=0.18,
                      ymax=0.22)
    elif problemname == 'drivencavity':
        podcoo = dict(xmin=0.45,
                      xmax=0.55,
                      ymin=0.7,
                      ymax=0.8)
    else:
        podcoo = femp['odcoo']

    # description of the control and observation domains
    dmd = femp['cdcoo']
    xmin, xmax, ymin, ymax = dmd['xmin'], dmd['xmax'], dmd['ymin'], dmd['ymax']
    velcondomstr = 'vel control domain: [{0}, {1}]x[{2}, {3}]\n'.\
        format(xmin, xmax, ymin, ymax)
    dmd = femp['odcoo']
    xmin, xmax, ymin, ymax = dmd['xmin'], dmd['xmax'], dmd['ymin'], dmd['ymax']
    velobsdomstr = 'vel observation domain: [{0}, {1}]x[{2}, {3}]\n'.\
        format(xmin, xmax, ymin, ymax)
    dmd = podcoo
    xmin, xmax, ymin, ymax = dmd['xmin'], dmd['xmax'], dmd['ymin'], dmd['ymax']
    pobsdomstr = 'pressure observation domain: [{0}, {1}]x[{2}, {3}]\n'.\
        format(xmin, xmax, ymin, ymax)

    pcmat = cou.get_pavrg_onsubd(odcoo=podcoo, Q=femp['Q'])

    cdatstr = snu.get_datastr_snu(time=None, meshp=N, nu=None, Nts=None)

    (coors, xinds,
     yinds, corfunvec) = dts.get_dof_coors(femp['V'], invinds=invinds)
    if bccontrol:
        bcctrstr = ' Boundary control is realized Robin penalization.' +\
            ' The setup is described in the accompanying preprint.' +\
            ' The implementation is in `https://github.com/highlando/' +\
            ' dolfin_navier_scipy/dolfin_navier_scipy/problem_setups.py` \n\n'
    else:
        bcctrstr = ''

    ctrl_visu_str = \
        ' the (distributed) control setup is as follows \n' +\
        ' B maps into the domain of control --' +\
        velcondomstr +\
        ' -- the first half of the columns' +\
        'actuate in x-direction, the second in y direction \n' +\
        ' Cv measures averaged velocities in the domain of observation' +\
        velobsdomstr +\
        ' Cp measures the averaged pressure' +\
        ' in the domain of pressure observation: ' +\
        pobsdomstr +\
        ' the first components are in x, the last in y-direction \n\n' +\
        bcctrstr +\
        ' Visualization: \n\n' +\
        ' `coors`   -- array of (x,y) coordinates in ' +\
        ' the same order as v[xinds] or v[yinds] \n' +\
        ' `xinds`, `yinds` -- indices of x and y components' +\
        ' of v = [vx, vy] -- note that indexing starts with 0\n' +\
        ' for testing use corfunvec wich is the interpolant of\n' +\
        ' f(x,y) = [x, y] on the grid \n\n' +\
        'Created in `get_exp_nsmats.py` ' +\
        '(see https://github.com/highlando/dolfin_navier_scipy) at\n' +\
        datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
    # TODO: this github address is outdated
    hstr = mddir + problemname + '_N{0}_hmat'.format(N)
    if bccontrol:
        hstr = hstr + '__bccontrol'
    try:
        hmat = dou.load_spa(hstr)
        print 'loaded `hmat`'
    except IOError:
        print 'assembling hmat ...'
        hmat = dts.ass_convmat_asmatquad(W=femp['V'], invindsw=invinds)
        dou.save_spa(hmat, hstr)

    zerv = np.zeros((NV, 1))
    L, fv_conv, rhsbc_convbc = \
        snu.get_v_conv_conts(prev_v=zerv, V=femp['V'], invinds=invinds,
                             diribcs=femp['diribcs'], Picard=False)

    L1, _, _ = \
        snu.get_v_conv_conts(prev_v=zerv, V=femp['V'], invinds=invinds,
                             diribcs=femp['diribcs'], Picard=True)

    L2 = L - L1

    if bccontrol:
        equationstring = ' `M v.dt + Av + 1/a*Arob v+ [L1+L2]v + H*kron(v,v) - J^Tp = 1/a*Brob*ubc + Bu + fv + fv_diff + fv_conv` \n' +\
            ' and  `Jv = -fp_div`$ \n' +\
            ' where `a` is a penalization parameter' +\
            ' (here set to `a={0}`)'.format(palpha) +\
            ' and `ubc` is the boundary control \n\n'
    else:
        equationstring = ' `M v.dt + Av + [L1+L2]v + H*kron(v,v) - J^Tp = Bu + fv + fv_diff + fv_conv` \n' +\
            ' and  `Jv = -fp_div`$ \n\n'

    infostr = 'These are the coefficient matrices of the quadratic ' +\
        'formulation of the Navier-Stokes Equations \n for the ' +\
        problemname + ' to be used as \n\n' +\
        equationstring +\
        ' note that we have assumed that the viscosity `nu=1`. \n' +\
        ' to change to a different `nu`, one needs to premultiply' +\
        ' the diffusion operator `A` and `fv_diff` accordingly ' +\
        ' then, the Reynoldsnumber has to determined from the setup ' +\
        ' in this case the relation `Re = {0}/nu` may be used\n\n'.\
        format(femp['Re']) + ctrl_visu_str

    savematdict = dict(A=A, M=M, H=hmat, J=J, L1=L1, L2=L2, fv=fv, fp=fp,
                       fv_conv=-fv_conv, fv_diff=fv_diff, fp_div=fp_div,
                       p_ss_nse=p_ss_nse, v_ss_nse=v_ss_nse,
                       p_ss_stokes=p_ss_stokes, v_ss_stokes=v_ss_stokes,
                       B=b_mat, Cv=c_mat, Cp=pcmat, info=infostr,
                       contsetupstr=contsetupstr, datastr=cdatstr, coors=coors,
                       xinds=xinds, yinds=yinds, corfunvec=corfunvec)

    matssstr = '__mats_NV{0}_Re{1}'.format(NV, femp['Re'])
    if bccontrol:
        savematdict.update(Brob=Brob, Arob=Arob)
        matssstr = matssstr + '_bccontrol_palpha{0}'.format(palpha)
    scipy.io.savemat(mddir + data_prfx + matssstr, savematdict,
                     do_compression=True)
    print 'saved mats to: ' + mddir + data_prfx + matssstr
    # NOTE: no case distinction between boundary controlled and uncontrolled
    # case as all boundaries will be initialized with zero
    # and nonzero values are applied last
    if visu:
        jsn_vslztn_dct(diribcs=femp['diribcs'], V=femp['V'], Q=femp['Q'],
                       invinds=invinds, N=N,
                       fstring=mddir+'visualization_cylN{0}.jsn'.format(N))

    return


def get_vtx_dofs(V, whichdim=0):
    mesh = V.mesh()
    dofmap = V.dofmap()
    nvertices = mesh.ufl_cell().num_vertices()
    # For each cell these are indices of the cell dofs which
    # correspond to vertices=entities of dim 0
    indices = [dofmap.tabulate_entity_dofs(0, i)[whichdim]
               for i in range(nvertices)]

    vertex_2_dof = dict()
    [vertex_2_dof.update(dict(vd for vd in
                         zip(cell.entities(0),
                             dofmap.cell_dofs(cell.index())[indices])))
     for cell in dolfin.cells(mesh)]

    # For check: the dof values at vertices match those
    # extracted from the array
    vertex_indices, vertex_dofs = map(list, zip(*vertex_2_dof.iteritems()))

    return vertex_dofs


def manual_paraview_data(vvec=None, V=None, invinds=None, diribcs=None,
                         fstr=None):
    import itertools

    # vfun = dolfin.Function(V)
    if diribcs is None:
        vaux = vvec
    else:
        vaux = np.zeros((V.dim(), 1))
        # fill in the boundary values
        for bc in diribcs:
            bcdict = bc.get_boundary_values()
            vaux[bcdict.keys(), 0] = bcdict.values()
        vaux[invinds] = vvec

    vxvtxdofs = get_vtx_dofs(V, whichdim=0)  # the x-dims
    vyvtxdofs = get_vtx_dofs(V, whichdim=1)  # the y-dims

    if fstr is None:
        for xvtx, yvtx in itertools.izip(vxvtxdofs, vyvtxdofs):
            print vaux[xvtx], vaux[yvtx], 0.
    else:
        wfile = file(fstr, 'w')
        for xvtx, yvtx in itertools.izip(vxvtxdofs, vyvtxdofs):
            wfile.write('{0} {1} {2} '
                        .format(vaux[xvtx][0], vaux[yvtx][0], 0.))

    return


def jsn_vslztn_dct(diribcs=None, V=None, N=None, invinds=None,
                   fstring=None, Q=None):
    ''' data for writing paraview output

    Notes:
    ---
    Todos:
     * pressure pinning
    '''

    import json

    vdim = V.dim()
    bclist = []
    for bc in diribcs:
        bcdict = bc.get_boundary_values()
        bclist.append(bcdict)

    vxvtxdofs = get_vtx_dofs(V, whichdim=0)  # the x-dims
    vyvtxdofs = get_vtx_dofs(V, whichdim=1)  # the y-dims
    pvtxdofs = get_vtx_dofs(Q, whichdim=0)  # the y-dims

    # convert to integers (fixes issues with json serialization)
    vxvtxdofs = [int(vtxdf) for vtxdf in vxvtxdofs]
    vyvtxdofs = [int(vtxdf) for vtxdf in vyvtxdofs]
    pvtxdofs = [int(ptxdf) for ptxdf in pvtxdofs]

    vtuheader_v = open('vel_vtuheader_cyl_N{0}.txt'.format(N), 'r').read()
    vtuheader_p = open('p_vtuheader_cyl_N{0}.txt'.format(N), 'r').read()

    vtufooter_v = '</DataArray> </PointData> </Piece>' +\
        ' </UnstructuredGrid> </VTKFile>'
    vtufooter_p = '</DataArray> </PointData> </Piece>' +\
        ' </UnstructuredGrid> </VTKFile>'

    vsldct = dict(bclist=bclist, invinds=invinds.tolist(), vdim=vdim,
                  vtufooter_v=vtufooter_v, vtuheader_v=vtuheader_v,
                  vtufooter_p=vtufooter_p, vtuheader_p=vtuheader_p,
                  vyvtxdofs=vyvtxdofs, vxvtxdofs=vxvtxdofs, pvtxdofs=pvtxdofs)
    jsfile = open(fstring, mode='w')
    jsfile.write(json.dumps(vsldct))
    print 'Visualization data dumped to json file: ', fstring


def apply_massinv(M, rhsa, output=None):
    """ Apply the inverse of mass or any other spd matrix

    to a rhs array
    TODO: by now just a wrapper for spsla.spsolve
    change e.g. to CG

    Parameters
    ----------
    M : (N,N) sparse matrix
        symmetric strictly positive definite
    rhsa : (N,K) ndarray array or sparse matrix
        array the inverse of M is to be applied to
    output : string, optional
        set to 'sparse' if rhsa has many zero columns
        to get the output as a sparse matrix

    Returns
    -------
    , : (N,K) ndarray or sparse matrix
        the inverse of `M` applied to `rhsa`

    """

    if output == 'sparse':
        colinds = rhsa.tocsr().indices
        colinds = np.unique(colinds)
        rhsa_cpy = rhsa.tolil()
        for col in colinds:
            rhsa_cpy[:, col] = np.atleast_2d(spsla.spsolve(M,
                                             rhsa_cpy[:, col])).T
        return rhsa_cpy

    else:
        mlusolve = spsla.factorized(M.tocsc())
        try:
            mirhs = np.copy(rhsa.todense())
        except AttributeError:
            mirhs = np.copy(rhsa)

        for ccol in range(mirhs.shape[1]):
            mirhs[:, ccol] = mlusolve(mirhs[:, ccol])

        return mirhs

if __name__ == '__main__':
    mddir = '../data/'
    # comp_exp_nsbili(problemname='drivencavity', N=10, mddir=mddir)
    # comp_exp_nsbili(problemname='drivencavity', N=20, mddir=mddir)
    # comp_exp_nsbili(problemname='drivencavity', N=30, mddir=mddir)
    # comp_exp_nsbili(problemname='cylinderwake', N=1, mddir=mddir, Re=80)
    # comp_exp_nsbili(problemname='cylinderwake', N=1, mddir=mddir, Re=40)
    # comp_exp_nsbili(problemname='cylinderwake', N=1, mddir=mddir)
    # comp_exp_nsbili(problemname='cylinderwake', N=3, mddir=mddir)
    comp_exp_nsbili(problemname='cylinderwake', N=2,
                    mddir=mddir, bccontrol=True, palpha=1)
    # comp_exp_nsbili(problemname='cylinderwake', N=2, Re=40,
    #                 mddir=mddir, bccontrol=True, palpha=1e-3)
