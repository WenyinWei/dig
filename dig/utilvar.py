import numpy as np

import MDSplus
from dig import MDSplus_server_IP

def time_linear_weighting(tdataseries, dataseries, tpoint):
    I = np.searchsorted(tdataseries, tpoint, ) - 1
    tslice_1 = tdataseries[I] 
    tslice_2 = tdataseries[I+1] 
    trange = tslice_2 - tslice_1
    return dataseries[I] * ((tslice_2 - tpoint) / trange) + dataseries[I+1] *  ((tpoint - tslice_1) / trange) # simple linear weighting

def rcentr_bcentr_to_Bt(rcentr, bcentr, R, Z):
    """from efit rcentr, bcentr and R,Z mesh to Bt

    Args:
        rcentr (_type_): _description_
        bcentr (_type_): _description_
        R (_type_): _description_
        Z (_type_): _description_

    Returns:
        Bt (np.ndarray): 2D B toroidal field in [iR, iZ]
    """
    Bt = np.ones( (R.size, Z.size) )
    Bt = rcentr * bcentr * Bt / R[None,:]
    return Bt.T

        

def get_EFIT_BR_BZ_BPhi(machine:str, shotnum:int, tpoints:list=None):
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
    if machine != "EAST":
        raise NotImplementedError("Only EAST is supported now")
    
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
            psi_ = psi[:,:,i].T # psi_ in [iZ, iR], psi in [iR, iZ, iT]

            from numpy import gradient
            # Eq. (3) W Zwingmann et al Equilibrium analysis of tokamak discharges with anisotropic pressure Plasma Phys. Control. Fusion (2001) vol. 43 1441
            dpsidR, dpsidZ = gradient(psi_.T, R, Z) # dpsidR, dpsidZ in [iR, iZ]
            dpsidR, dpsidZ = dpsidR.T, dpsidZ.T 
            RR = np.matmul( np.ones_like(Z)[:,None], R[None,:] )
            BZ = -dpsidR/RR # [iZ, iR]
            BR = dpsidZ/RR 

            BRs.append(BR.T) # BR.T in [iR, iZ]
            BZs.append(BZ.T)
    else: # does nto use time given by efit directly
        for tpoint in tpoints:
            I = np.searchsorted(tefit, tpoint, ) - 1 
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

def get_EFIT_psi(machine:str, shotnum:int, tpoints:float=None):
    """psi distribution interpolated on the specified time points.

    psi = psis[0]
    psi[ R_index, Z_index ]
    

    Args:
        shotnum (int): [description]
        tpoints (list of float): [description]

    Returns:
        R: on which R grid is the BR, BZ, Bt
        Z: on which Z grid is the BR, BZ, Bt
        psis (list of np.ndarray): a list of psi[R_ind, Z_ind] distribution, the time of which are corresponding to the specified tpoints (if not specified, all is given)
    """
    if machine != "EAST":
        raise NotImplementedError("Only EAST is supported now")
    
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
            I = np.searchsorted(tefit, tpoint, ) - 1
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


def get_tRMP_IRMP(machine:str, shotnum:int, mds_treedict=None):
    """To get EAST RMP (resonant magnetic perturbation) coil time and current sequences.

    Args:
        shotnum (int): Shot number on EAST

    Returns:
        mds_treedict (dict of dict): mds_treedict["east"] has the keys "tRMP", "IRMPU1", ... , "IRMPU8", "IRMPL1", ... , "IRMPL8", all the values of which are resepectively a 1D numpy array.
        
    Plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1,1)
        for i in range(1, 8+1):
            ax.plot( mds_treedict["east"][f"tRMP"], mds_treedict["east"][f"IRMPU{i}"], label=f"IRMPU{i}" )
        fig.legend()
    """
    if machine != "EAST":
        raise NotImplementedError("Only EAST is supported now")
    
    mds_conn = MDSplus.connection.Connection(MDSplus_server_IP)

    tree = "east"
    mds_conn.openTree(tree, shotnum)
    if mds_treedict is None:
        mds_treedict = dict()
    mds_treedict[tree] = dict()
    mds_treedict[tree][f"tRMP"] = mds_conn.get(f"dim_of(\\IRMPL1)").data()
    for i in range(1, 8 +1):
        mds_treedict[tree][f"IRMPL{i}"] = mds_conn.get(f"\\IRMPL{i}").data()
        mds_treedict[tree][f"IRMPU{i}"] = mds_conn.get(f"\\IRMPU{i}").data()
    mds_conn.closeTree(tree, shotnum)

    mds_conn.closeAllTrees()
    mds_conn.disconnect()
    
    return mds_treedict
    
    
def get_limiter(machine:str, shotnum:int, mds_treedict=None):
    if machine != "EAST":
        raise NotImplementedError("Only EAST is supported now")
    
    mds_conn = MDSplus.connection.Connection(MDSplus_server_IP)

    tree = "efit_east"
    if mds_treedict is None:
        mds_treedict = dict()
    mds_treedict[tree] = dict()
    mds_conn.openTree(tree, shotnum)
    mds_treedict[tree]["lim"] = mds_conn.get("\\lim").data()
    mds_conn.closeTree(tree, shotnum)

    mds_conn.closeAllTrees()
    mds_conn.disconnect()
    
    return mds_treedict
    
def get_tPF_IPF(machine:str, shotnum:int, mds_treedict=None):
    if machine != "EAST":
        raise NotImplementedError("Only EAST is supported now")

    mds_conn = MDSplus.connection.Connection(MDSplus_server_IP)

    tree = "east"
    if mds_treedict is None:
        mds_treedict = dict()
    mds_treedict[tree] = dict()
    mds_conn.openTree(tree, shotnum)
    mds_treedict[tree]["tPF"] = mds_conn.get(f"dim_of(\\PF1P)").data()
    for i in range(1, 8+1):
        mds_treedict[tree][f"IPF{i}"] = mds_conn.get(f"\\PF{i}P").data()
    mds_treedict[tree]["IPF9"] = mds_treedict[tree]["IPF7"]
    mds_treedict[tree]["IPF10"] = mds_treedict[tree]["IPF8"]
    for i in range(11, 14+1):
        mds_treedict[tree][f"IPF{i}"] = mds_conn.get(f"\\PF{i-2}P").data()
    mds_conn.closeTree(tree, shotnum)

    mds_conn.closeAllTrees()
    mds_conn.disconnect()
    
    return mds_treedict
    