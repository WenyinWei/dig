import numpy as np

import MDSplus
from dig import MDSplus_server_IP

def rcentr_bcentr_to_Bt(rcentr, bcentr, R, Z):
    """from efit rcentr, bcentr and R,Z mesh to Bt

    Args:
        rcentr (_type_): _description_
        bcentr (_type_): _description_
        R (_type_): _description_
        Z (_type_): _description_

    Returns:
        Bt (np.ndarray): 2D B toroidal field in [iZ, iR]
    """
    Bt = np.ones( (R.size, Z.size) )
    Bt = rcentr * bcentr * Bt / R[None,:]
    return Bt

def R_Z_psirz_to_BR_BZ(R, Z, psirz):
    """input R, Z are 1D mesh, psirz[iZ,iR] is psi 2D distribution

    Args:
        R (_type_): _description_
        Z (_type_): _description_
        psirz (_type_): _description_

    Returns:
        BR, BZ (np.ndarray): 2D mesh [iZ, iR] in matlab index style 
    """
    psirz_ = psirz.T
    from numpy import gradient
    dpsidR, dpsidZ = gradient(psirz_.T, R, Z)
    dpsidR, dpsidZ = dpsidR.T, dpsidZ.T
    RR = np.matmul( np.ones_like(Z)[:,None], R[None,:] )
    BZ = -dpsidR/RR 
    BR = dpsidZ/RR 
    return BR, BZ
        

def get_EAST_EFIT_BR_BZ_BPhi(shotnum:int, tpoints:list=None):
    """ Return the list of BR, BZ, Bt on the specified tpoints. If tpoints unspecified, all BR, BZ, Bt on time points given by efit would be calculated. 
    Note the first index is for R while the second one is for R, 
    BR = BRs[0]
    BR[ R_index, Z_index ]

    Args:
        shotnum (int): [description]
        tpoints (list of float): [description]

    Returns:
        R: on which R grid is the BR, BZ, Bt
        Z: on which Z grid is the BR, BZ, Bt
        tpoints (list of np.ndarray): on which time points the BR,BZ,Bt in BRs, BZs, Bts are
        BRs (list of np.ndarray): B R component in cylindrical coordinate, 
        BZs (list of np.ndarray): B Z component in cylindrical coordinate
        Bts (list of np.ndarray): B Phi component (toroidal field) in cylindrical coordinate
    """
    R0 = 1.75 # Make sure this is right for your fusion machine
    if tpoints is None:
        whether_use_time_of_efit = True
    else:
        whether_use_time_of_efit = False
    mds_conn = MDSplus.connection.Connection(MDSplus_server_IP)

    # mds_conn.openTree("pcs_east", shotnum)
    # tdRsep = mds_conn.get("dim_of(\\idsdrsep)").T
    # dRsep_0 = mds_conn.get("\\idsdrsep").T
    # dRsep = np.interp(tpoint, tdRsep, dRsep_0) 
    # mds_conn.closeTree("pcs_east", shotnum)

    # -------------------------------------------------------------------------
    mds_conn.openTree("efit_east", shotnum)



    tefit = mds_conn.get("\\time").data().T  # s
    psi = mds_conn.get("\\psirz").data().T
    R = mds_conn.get("\\R").data().T
    Z = mds_conn.get("\\Z").data().T
    BRs, BZs, Bts = [], [], []
    if whether_use_time_of_efit:
        tpoints = tefit
        for i in range( len(tpoints) ):
            psi_ = psi[:,:,i].T

            from numpy import gradient
            dpsidR, dpsidZ = gradient(psi_.T, R, Z)
            dpsidR, dpsidZ = dpsidR.T, dpsidZ.T 
            RR = np.matmul( np.ones_like(Z)[:,None], R[None,:] )
            BZ = -dpsidR/RR 
            BR = dpsidZ/RR 

            BRs.append(BR.T)
            BZs.append(BZ.T)
    else: # does nto use time given by efit directly
        for tpoint in tpoints:
            I = np.searchsorted(tefit, tpoint, side='right')
            tslice_1 = tefit[I] 
            tslice_2 = tefit[I+1] 
            trange = tslice_2 - tslice_1

            # yq95 = mds_conn.get("\\Q95").data().T
            # q95 = yq95[I]
            psi_1 = psi[:,:,I].T
            psi_2 = psi[:,:,I+1].T

            psi_ = psi_1 * ((tslice_2 - tpoint) / trange) + psi_2 *  ((tpoint - tslice_1) / trange) # simple linear weighting

            # aminor = mds_conn.get("\\aminor").data().T/100  # m
            # aminor = aminor[I] 
            # Rmaxis = mds_conn.get("\\Rmaxis").data().T  # m
            # rmaxis = Rmaxis[I] 
            # Zmaxis = mds_conn.get("\\Zmaxis").data().T  # m
            # zmaxis = Zmaxis[I] 

            from numpy import gradient
            dpsidR, dpsidZ = gradient(psi_.T, R, Z)
            dpsidR, dpsidZ = dpsidR.T, dpsidZ.T
            RR = np.matmul( np.ones_like(Z)[:,None], R[None,:] )
            BZ = -dpsidR/RR 
            BR = dpsidZ/RR 

            BRs.append(BR.T)
            BZs.append(BZ.T)

    mds_conn.closeTree("efit_east", shotnum)

    #-------------------It
    mds_conn.openTree("east_1", shotnum)
    try:
        temp = mds_conn.get("\\focs4").data().T
    except:
        mds_conn.openTree("pcs_east", shotnum) 
        temp = mds_conn.get("\\it").data().T
        mds_conn.closeTree("pcs_east", shotnum) 
    it= np.mean(temp) # Positive current means anti-clockwise direction viewed from above
    Bt0 = it/4086.0 # 4086 is EAST-specific
    for _ in tpoints: # considered as constant now
        Bt = np.ones_like(BR)
        Bt = R0 * Bt0 * Bt / R[None,:]
        Bts.append(Bt.T)

    mds_conn.closeTree("east_1", shotnum)
    mds_conn.closeAllTrees()
    mds_conn.disconnect()
    return R, Z, tpoints, BRs, BZs, Bts

def get_EAST_EFIT_psi(shotnum:int, tpoints:float=None):
    """psi distribution interpolated on the specified time points.

    psi = psis[0]
    psi[ R_index, Z_index ]
    

    Args:
        shotnum (int): [description]
        tpoints (list of float): [description]

    Returns:
        R: on which R grid is the BR, BZ, Bt
        Z: on which Z grid is the BR, BZ, Bt
        psis (list of np.ndarray): a list of psi[Z_ind, R_ind] distribution, the time of which are corresponding to the specified tpoints (if not specified, all is given)
    """
    if tpoints is None:
        whether_use_time_of_efit = True
    else:
        whether_use_time_of_efit = False

    mds_conn = MDSplus.connection.Connection(MDSplus_server_IP)
    mds_conn.openTree("efit_east", shotnum)

    tefit = mds_conn.get("\\time").data().T  # s
    psi = mds_conn.get("\\psirz").data().T
    psis = []
    if whether_use_time_of_efit:
        tpoints = tefit
        for i in range( len(tpoints) ):
            psi_ = psi[:,:,i].T
            psis.append(psi_)
    else:
        for tpoint in tpoints:
            I = np.searchsorted(tefit, tpoint, side='right')
            tslice_1 = tefit[I] 
            tslice_2 = tefit[I+1] 
            trange = tslice_2 - tslice_1

            psi_1 = psi[:,:,I].T
            psi_2 = psi[:,:,I+1].T
            psi_ = psi_1 * ((tslice_2 - tpoint) / trange) + psi_2 *  ((tpoint - tslice_1) / trange) # simple linear weighting

            psis.append(psi_.T)

    R = mds_conn.get("\\R").data().T
    Z = mds_conn.get("\\Z").data().T

    mds_conn.closeAllTrees()
    mds_conn.disconnect()

    return R, Z, tpoints, psis

