import warnings
import numpy as np

import MDSplus
from dig import MDSplus_server_IP

from deprecated import deprecated

def twin2id(t, twin):
    if twin:
        ta = twin[0]
        tb = twin[1]
        Ia = np.argmin( np.abs(t-ta) )
        Ib = np.argmin( np.abs(t-tb) )+1
    else:
        # Ia = np.min( np.abs(t-0) )
        Ia = 0
        Ib = len(t)
    return [Ia, Ib]

def get_MDSplus_t_signal(shotnum:int, signal:str, tree:str, twin:list = [], ):
    """read signal from MDSplus according to the specified shotnum, tree name and signal name.

    Args:
        shotnum (int): which shot you want to get the data from.
        signal (str): which signal of the shot you want.
        tree (str): which MDSplus tree is the signal on. 
        twin (list, optional): the wanted time range of the data (time between, twin). Note that we crop the signal array according to the *downloaded* time array so a small time range won't fasten the download. Defaults to [].

    Returns:
        [np.ndarray, np.ndarray]: [the time array in twin, the signal array in twin]
    """    
    if 'east' in tree.lower() or 'east_1' in tree.lower():
        removebackground = True

    mds_conn = MDSplus.connection.Connection(MDSplus_server_IP)
    mds_conn.openTree(tree, shotnum)

    try:
        if 'efit_east' in tree.lower() or 'efitrt_east' in tree.lower():
            t = mds_conn.get("\\time").data().T
        else:
            t = mds_conn.get(f"dim_of(\\{signal})").data().T
    except Exception as e:
        print()
        print('%----------------------------------------------%')
        print(e)
        warnings.warn(f"#{shotnum} {signal} signal not found!!")
        print('%----------------------------------------------%')
        print()
        tx, x = None, None
        return tx, x
        
    x = mds_conn.get(f"\\{signal}").data().T # add '\' before signal string
    
    if twin:
        ta = twin[0]
        tb = twin[1]
        Ia = np.argmin( np.abs(t-ta) )
        Ib = np.argmin( np.abs(t-tb) )+1
    else:
        # Ia = np.min( np.abs(t-0) )
        Ia = 0
        Ib = len(t)
    
    tx = t[Ia:Ib]
    if 'efit_east' in tree.lower() or 'pcs_east' in tree.lower():
        x = x[Ia:Ib]
    else:
        x = x[Ia:Ib]

    mds_conn.closeTree(tree, shotnum)
    mds_conn.disconnect()

    return tx, x

def get_MDSplus_signal(shotnum:int, signal:str, tree:str, twinid:list = [], ):
    """read signal from EAST MDSplus according to the specified shotnum, tree name and signal name.

    Args:
        shotnum (int): which shot you want to get the data from.
        signal (str): which signal of the shot you want.
        tree (str): which MDSplus tree is the signal on. 
        twin (list, optional): the wanted time range of the data (time between, twin). Note that we crop the signal array according to the *downloaded* time array so a small time range won't fasten the download. Defaults to [].
        removebackground (bool, optional): whether or not to remove the background. Defaults to False.

    Returns:
        [np.ndarray, np.ndarray]: [the time array in twin, the signal array in twin]
    """
    if 'east' in tree.lower() or 'east_1' in tree.lower():
        removebackground = True

    mds_conn = MDSplus.connection.Connection(MDSplus_server_IP)
    mds_conn.openTree(tree, shotnum)
    x = mds_conn.get('\\' + signal).data().T # add '\' before signal string
    mds_conn.closeAllTrees()
    mds_conn.disconnect()

    if twinid:
        Ia, Ib = twinid[0], twinid[1]
        x = x[Ia:Ib]
    return None, x
