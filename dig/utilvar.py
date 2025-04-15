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
    psi = mds_conn.get("\\psirz").data().T # poloidal flux over 2π, indexes [iR, iZ, iT]
    R = mds_conn.get("\\R").data().T
    Z = mds_conn.get("\\Z").data().T
    FPOL = mds_conn.get('\\FPOL').data().T # Poloidal current, Fpol = R * B_phi, indexes  [ipsi, it]
    PSIMAG = mds_conn.get('\\SSIMAG').data().T # psi at the magnetic axis, index [it]
    PSIBRY = mds_conn.get('\\SSIBRY').data().T # psi at the boundary, index [it]
    ZXPT1 = mds_conn.get('\\ZXPT1').data().T # Lower X point Z coordinate, index [it]
    ZXPT2 = mds_conn.get('\\ZXPT2').data().T # Upper X point Z coordinate, index [it]
    from numpy import gradient
    # Compute the grid spacing
    dR = np.gradient(R)
    dZ = np.gradient(Z)

    BRs, BZs, Bts = [], [], []
    if whether_use_time_of_efit:
        tpoints = tefit
        for i in range( len(tpoints) ):
            psi_ = psi[:,:,i] 

            # Compute the first derivative of ψ with respect to R
            ppsi_pR = np.gradient(psi_, axis=0) / dR[:, None]
            # Compute the first derivative of ψ with respect to Z
            ppsi_pZ = np.gradient(psi_, axis=1) / dZ[None, :]

            BZ = ppsi_pR/R[:,None]
            BR = - ppsi_pZ/R[:,None] 

            BRs.append(BR)
            BZs.append(BZ)
    else: # does nto use time given by efit directly
        for tpoint in tpoints:
            I = np.searchsorted(tefit, tpoint, ) - 1 
            tslice_1 = tefit[I] 
            tslice_2 = tefit[I+1] 
            assert tslice_1 <= tpoint <= tslice_2, f"Time point {tpoint} is not between {tslice_1} and {tslice_2}"
            trange = tslice_2 - tslice_1

            # yq95 = mds_conn.get("\\Q95").data().T
            # q95 = yq95[I]
            psi_1 = psi[:,:,I]
            psi_2 = psi[:,:,I+1]

            psi_ = psi_1 * ((tslice_2 - tpoint) / trange) + psi_2 *  ((tpoint - tslice_1) / trange) # simple linear weighting

            # aminor = mds_conn.get("\\aminor").data().T/100  # m
            # aminor = aminor[I] 
            # Rmaxis = mds_conn.get("\\Rmaxis").data().T  # m
            # rmaxis = Rmaxis[I] 
            # Zmaxis = mds_conn.get("\\Zmaxis").data().T  # m
            # zmaxis = Zmaxis[I] 

            # Compute the first derivative of ψ with respect to R
            ppsi_pR = np.gradient(psi_, axis=0) / dR[:, None]
            # Compute the first derivative of ψ with respect to Z
            ppsi_pZ = np.gradient(psi_, axis=1) / dZ[None, :]

            BZ = ppsi_pR/R[:,None]
            BR = - ppsi_pZ/R[:,None] 

            BRs.append(BR)
            BZs.append(BZ)

    from scipy.interpolate import interp1d
    if whether_use_time_of_efit:
        for it in range( len(tpoints) ):
            psi_arr_ = np.linspace(PSIMAG[it], PSIBRY[it], num=129, endpoint=True)
            FPOL_RZ = interp1d(psi_arr_, FPOL[:, it], bounds_error=False, fill_value=(FPOL[0,it], FPOL[-1,it]))(psi[:,:,it])
            for iZ in range(Z.size):
                if (Z[iZ] < ZXPT1[it]) or (Z[iZ] > ZXPT2[it]):
                    FPOL_RZ[:,iZ] = FPOL[-1,it]
            Bt = FPOL_RZ[:,:] / R[:,None]
            Bts.append(Bt)
    else: # does not use time given by efit directly
        for tpoint in tpoints:
            I = np.searchsorted(tefit, tpoint, ) - 1 
            tslice_1 = tefit[I] 
            tslice_2 = tefit[I+1] 
            assert tslice_1 <= tpoint <= tslice_2, f"Time point {tpoint} is not between {tslice_1} and {tslice_2}"
            trange = tslice_2 - tslice_1
            PSIMAG_ = PSIMAG[I] * ((tslice_2 - tpoint) / trange) + PSIMAG[I+1] * ((tpoint - tslice_1) / trange) 
            PSIBRY_ = PSIBRY[I] * ((tslice_2 - tpoint) / trange) + PSIBRY[I+1] * ((tpoint - tslice_1) / trange) 
            psi_arr_ = np.linspace(PSIMAG_, PSIBRY_, num=129, endpoint=True)
            FPOL_ = FPOL[:,I] * ((tslice_2 - tpoint) / trange) + FPOL[:,I+1] * ((tpoint - tslice_1) / trange) 
            psi_ = psi[:,:,I] * ((tslice_2 - tpoint) / trange) + psi[:,:,I+1] * ((tpoint - tslice_1) / trange)
            FPOL_RZ = interp1d(psi_arr_, FPOL_, bounds_error=False, fill_value=(FPOL_[0], FPOL_[-1]))(psi_)
            ZXPT1_ = ZXPT1[I] * ((tslice_2 - tpoint) / trange) + ZXPT1[I+1] * ((tpoint - tslice_1) / trange)
            ZXPT2_ = ZXPT2[I] * ((tslice_2 - tpoint) / trange) + ZXPT2[I+1] * ((tpoint - tslice_1) / trange)
            for iZ in range(Z.size):
                if (Z[iZ] < ZXPT1_) or (Z[iZ] > ZXPT2_):
                    FPOL_RZ[:,iZ] = FPOL_[-1]
            Bt = FPOL_RZ[:,:]/ R[:,None] 
            Bts.append(Bt)
    mds_conn.closeTree("efit_east", shotnum)

    mds_conn.closeAllTrees()
    mds_conn.disconnect()
    return R, Z, tpoints, BRs, BZs, Bts

def get_EFIT_JR_JZ_JPhi(machine:str, shotnum:int, tpoints:list=None):
    """ Return the list of JR, JZ, JPhi on the specified tpoints. If tpoints unspecified, all JR, JZ, JPhi on time points given by efit would be calculated. 
    Note the first index is for R while the second one is for R,
    JR = JRs[0]
    JR[ R_index, Z_index ]
    Args:
        shotnum (int): [description]
        tpoints (list of float): [description]
    Returns:
        R: on which R grid is the JR, JZ, JPhi
        Z: on which Z grid is the JR, JZ, JPhi
        tpoints (list of np.ndarray): on which time points the JR,JZ,JPhi in JRs, JZs, JPhis are 
        JRs (list of np.ndarray): J R component in cylindrical coordinate,
        JZs (list of np.ndarray): J Z component in cylindrical coordinate
        JPhis (list of np.ndarray): J Phi component (toroidal field) in cylindrical coordinate
    """
    from scipy.interpolate import interp1d
    from dig.math import axisymmetric_laplace

    R, Z, _, BRs, BZs, BPhis = get_EFIT_BR_BZ_BPhi(machine, shotnum, tpoints)
    R, Z, _, psiRZs = get_EFIT_psi(machine, shotnum, tpoints)
    if machine != "EAST":
        raise NotImplementedError("Only EAST is supported now")
    
    mu_0 = 4 * np.pi * 1e-7 # vacuum permeability in [H/m]
    if tpoints is None:
        whether_use_time_of_efit = True
    else:
        whether_use_time_of_efit = False
    mds_conn = MDSplus.connection.Connection(MDSplus_server_IP)
    mds_conn.openTree("efit_east", shotnum)

    tefit = mds_conn.get("\\time").data().T  # s
    JRs, JZs, JPhis = [], [], []
    PPRIME = mds_conn.get('\\PPRIME').data().T # P' = dP/dpsi, indexes [ipsi, it]
    FFPRIM = mds_conn.get('\\FFPRIM').data().T # FF' = F * dF/dpsi, indexes [ipsi, it]
    FPOL = mds_conn.get('\\FPOL').data().T # Poloidal current, Fpol = R * B_phi, indexes [ipsi, it]
    FPRIM = FFPRIM / FPOL # F' = FF' / F, indexes [ipsi, it]

    PSIBRY = mds_conn.get('\\SSIBRY').data().T # psi at the boundary, index [it]
    PSIMAG = mds_conn.get('\\SSIMAG').data().T # psi at the magnetic axis, index [it]

    ZXPT1 = mds_conn.get('\\ZXPT1').data().T # Lower X point Z coordinate
    ZXPT2 = mds_conn.get('\\ZXPT2').data().T # Upper X point Z coordinate

    if whether_use_time_of_efit:
        tpoints = tefit
        for it in range( len(tpoints) ):
            psi_arr_ = np.linspace(PSIMAG[it], PSIBRY[it], num=129, endpoint=True)
            FPRIM_interpolator = interp1d(psi_arr_, FPRIM[:,it], bounds_error=False, fill_value=(FPRIM[0,it], 0.0))
            FPRIM_on_grid = FPRIM_interpolator(psiRZs[it][:,:])
            for iZ in range(Z.size):
                if (Z[iZ] < ZXPT1[it]) or (Z[iZ] > ZXPT2[it]):
                    FPRIM_on_grid[:,iZ] = 0.0

            JRs.append( (1.0 / mu_0) * BRs[it] * FPRIM_on_grid)
            JZs.append( (1.0 / mu_0) * BZs[it] * FPRIM_on_grid)

            JPhi = - (1/mu_0) * axisymmetric_laplace(psiRZs[it], R, Z) / R[:,None]
            JPhis.append( JPhi )

    else: # does not use time given by efit directly
        for it, tpoint in enumerate(tpoints):
            I = np.searchsorted(tefit, tpoint, ) - 1 
            tslice_1 = tefit[I] 
            tslice_2 = tefit[I+1] 
            assert tslice_1 <= tpoint <= tslice_2, f"Time point {tpoint} is not between {tslice_1} and {tslice_2}"
            Delta_t = tslice_2 - tslice_1
            PSIMAG_ = PSIMAG[I] * ((tslice_2 - tpoint) / Delta_t) + PSIMAG[I+1] * ((tpoint - tslice_1) / Delta_t) 
            PSIBRY_ = PSIBRY[I] * ((tslice_2 - tpoint) / Delta_t) + PSIBRY[I+1] * ((tpoint - tslice_1) / Delta_t) 
            psi_arr_ = np.linspace(PSIMAG_, PSIBRY_, num=129, endpoint=True)
            FPRIM_1 = FPRIM[:,I]
            FPRIM_2 = FPRIM[:,I+1]
            FPRIM_ = FPRIM_1 * ((tslice_2 - tpoint) / Delta_t) + FPRIM_2 *  ((tpoint - tslice_1) / Delta_t) 
            FPRIM_interpolator = interp1d(psi_arr_, FPRIM_, bounds_error=False, fill_value=(FPRIM_[0], 0.0))
            FPRIM_RZ_ = FPRIM_interpolator(psiRZs[it][:,:])
            for iZ in range(Z.size):
                if (Z[iZ] < ZXPT1[I]) or (Z[iZ] > ZXPT2[I]):
                    FPRIM_RZ_[:,iZ] = 0.0
            JRs.append( (1.0 / mu_0) * BRs[it] * FPRIM_RZ_)
            JZs.append( (1.0 / mu_0) * BZs[it] * FPRIM_RZ_)
            JPhi = - (1/mu_0) * axisymmetric_laplace(psiRZs[it], R, Z) / R[:,None]
            JPhis.append( JPhi )
    return R, Z, tpoints, JRs, JZs, JPhis

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
    psi = mds_conn.get("\\psirz").data().T # poloidal flux over 2π, indexes [iR, iZ, iT]
    psis = []
    if whether_use_time_of_efit:
        tpoints = tefit
        for i in range( len(tpoints) ):
            psi_ = psi[:,:,i]
            psis.append(psi_)
    else:
        for tpoint in tpoints:
            I = np.searchsorted(tefit, tpoint, ) - 1
            tslice_1 = tefit[I] 
            tslice_2 = tefit[I+1] 
            trange = tslice_2 - tslice_1

            psi_1 = psi[:,:,I]
            psi_2 = psi[:,:,I+1]
            psi_ = psi_1 * ((tslice_2 - tpoint) / trange) + psi_2 *  ((tpoint - tslice_1) / trange) # simple linear weighting

            psis.append(psi_)

    R = mds_conn.get("\\R").data().T
    Z = mds_conn.get("\\Z").data().T

    mds_conn.closeAllTrees()
    mds_conn.disconnect()

    return R, Z, tpoints, psis

def get_EFIT_RZMAGIS(machine:str, shotnum:int, tpoints:list=None):
    """Fetch the R, Z coordinates of axis data from EAST EFIT tree.

    Args:
        machine (str): Machine name
        shotnum (int): Shot number
        tpoints (list, optional): Time points to evaluate the data. Defaults to None, if so, the tpoints takes the EFIT time sequence.
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
    Rmagis = mds_conn.get("\\RMAGIS").data().T
    Zmagis = mds_conn.get("\\ZMAGIS").data().T

    Rmagis_list = []
    Zmagis_list = []
    if whether_use_time_of_efit:
        tpoints = tefit
        Rmagis_list = list(Rmagis)
        Zmagis_list = list(Zmagis)
    else:
        for tpoint in tpoints:
            I = np.searchsorted(tefit, tpoint, ) - 1
            tslice_1 = tefit[I] 
            tslice_2 = tefit[I+1] 
            trange = tslice_2 - tslice_1

            Rmagis_1 = Rmagis[I]
            Rmagis_2 = Rmagis[I+1]
            Zmagis_1 = Zmagis[I]
            Zmagis_2 = Zmagis[I+1]

            Rmagis_val = Rmagis_1 * ((tslice_2 - tpoint) / trange) + Rmagis_2 *  ((tpoint - tslice_1) / trange) # simple linear weighting
            Zmagis_val = Zmagis_1 * ((tslice_2 - tpoint) / trange) + Zmagis_2 *  ((tpoint - tslice_1) / trange) # simple linear weighting

            Rmagis_list.append(Rmagis_val)
            Zmagis_list.append(Zmagis_val)
    
    mds_conn.closeAllTrees()
    mds_conn.disconnect()

    return None, None, tpoints, Rmagis_list, Zmagis_list


def get_EFIT_qpsi(machine:str, shotnum:int, tpoints:list=None):
    """Fetch the qpsi data from EAST EFIT tree.

    Args:
        machine (str): Machine name
        shotnum (int): Shot number
        tpoints (list, optional): Time points to evaluate the data. Defaults to None, if so, the tpoints takes the EFIT time sequence.
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
    qpsi = mds_conn.get("\\qpsi").data().T # qpsi, indexes [ipsi, iT]

    qpsis = []
    if whether_use_time_of_efit:
        tpoints = tefit
        for i in range( len(tpoints) ):
            qpsi_ = qpsi[:,i]
            qpsis.append(qpsi_)
    else:
        for tpoint in tpoints:
            I = np.searchsorted(tefit, tpoint, ) - 1
            tslice_1 = tefit[I] 
            tslice_2 = tefit[I+1] 
            trange = tslice_2 - tslice_1

            qpsi_1 = qpsi[:,I]
            qpsi_2 = qpsi[:,I+1]
            qpsi_ = qpsi_1 * ((tslice_2 - tpoint) / trange) + qpsi_2 *  ((tpoint - tslice_1) / trange)
            qpsis.append(qpsi_)
    mds_conn.closeAllTrees()
    mds_conn.disconnect()
    return None, None, tpoints, qpsis


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
    