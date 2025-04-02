import numpy as np

def axisymmetric_laplace(psiRZ, R, Z):
    """
    Compute the axisymmetric Laplace operator:
    R * d/dR (1/R * dψ/dR) + d²ψ/dZ²

    Parameters:
    - psiRZ: 2D numpy array of shape [iR, iZ], the ψ values.
    - R: 1D numpy array of shape [iR], the R grid points.
    - Z: 1D numpy array of shape [iZ], the Z grid points.

    Returns:
    - result: 2D numpy array of shape [iR, iZ], the result of the operator.
    """
    # Compute the grid spacing
    dR = np.gradient(R)
    dZ = np.gradient(Z)

    # Compute the first derivative of ψ with respect to R
    ppsi_pR = np.gradient(psiRZ, axis=0) / dR[:, None]

    # Compute the second derivative of ψ with respect to R
    p2psi_pR2 = np.gradient(ppsi_pR, axis=0) / dR[:, None]

    # Compute the first derivative of ψ with respect to Z
    ppsi_pZ = np.gradient(psiRZ, axis=1) / dZ[None, :]

    # Compute the second derivative of ψ with respect to Z
    p2psi_pZ2 = np.gradient(ppsi_pZ, axis=1) / dZ[None, :]

    return p2psi_pR2 - ppsi_pR / R[:,None] + p2psi_pZ2