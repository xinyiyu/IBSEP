import os, sys, copy
import pandas as pd
import numpy as np
import logging
from tqdm import tqdm
from scipy import linalg
import scipy.stats as st
from scipy.stats import norm
import gc
import warnings
warnings.filterwarnings('ignore') 

## IBSEP estimator
def TensorBlue(betas, ses, props, Omega, Sigma, ld):
    """
    Runs the core method to combine results and generate final, combined summary statistics
    :param betas: Mx(P+1) matrix of beta values (M = # of SNPs, P = # of cell types)
    :param ses: Mx(P+1) matrix of se values 
    :param props: P array of cell type proportion mean values
    :param Omega: PxP genetic covariance matrix
    :param Sigma: (P+1)x(P+1) matrix of ldsc intercepts
    :param ld: M or Mx1 ldscore
    :return: Tuple containing:
                 1) Result ndarray of betas (MxP)
                 2) Result ndarray of beta standard errors (MxP)
    """

    M, P_ = betas.shape
    P = P_ - 1
    for i in range(P):
        if Omega[i, i] == 0:
            Omega[i, i] = 1e-16
    ld = ld.squeeze() # (M,)
    omega = ld[:, np.newaxis, np.newaxis] * Omega[np.newaxis, :, :] # Omega_j for all snps (M, P, P)
    ses_diag = np.eye(P+1) * ses[:, np.newaxis, :] # (M, P+1, P+1)
    sigma = ses_diag @ Sigma @ ses_diag # SCS for all snps, (M, P+1, P+1) 
    A = np.vstack([np.eye(P), props]) # (P+1, P)

    # Create a 3D matrix, M rows of Px1 column vectors with shape (M, P, 1)
    d_indices = np.arange(P)
    # print(omega.shape, d_indices)
    omega_diag = omega[:, d_indices, d_indices][:, :, np.newaxis] # (M, P, 1)
    omega_pp_scaled = np.divide(omega, omega_diag) # Normalize Omega_j, (M, P, P)
    lamb = A[np.newaxis, :, :] @ omega # lambda for all snps (M, P+1, P)
    lamb_diag = omega[:, d_indices, d_indices][:, np.newaxis, :] # (M, 1, P)
    lamb_pp_scaled = np.divide(lamb, lamb_diag)  # by colomn, (M, P+1, P)
    lamb_pp_scaled = np.transpose(lamb_pp_scaled, (0, 2, 1)) # (M, P, P+1)

    # Produce center matrix in steps (product of omega terms, add omega and sigma, then invert)
    center_matrix_inv = omega[:, np.newaxis, :, :] - omega_pp_scaled[:, :, :, np.newaxis] * omega[:, :, np.newaxis, :] # Lambda^-1 (M, P, P, P)
    center_matrix_inv = A[np.newaxis, np.newaxis, :, :] @ center_matrix_inv @ np.transpose(A[np.newaxis, np.newaxis, :, :], (0, 1, 3, 2)) # (M, P, P+1, P+1)
    center_matrix_inv += sigma[:, np.newaxis, :, :] # (M, P, P+1, P+1)
    center_matrix = np.linalg.inv(center_matrix_inv) # Inverts each slice separately
    del center_matrix_inv  # Clean up the inverse matrix to free space
    gc.collect()

    # Calculate lambda * center_matrix
    left_product = np.matmul(lamb_pp_scaled[:, :, np.newaxis, :], center_matrix) # (M, P, 1, P+1)
    del center_matrix  # Clean up the center matrix to free space
    gc.collect()

    # Calculate denominator (M x P x 1 x 1)
    denom = np.matmul(left_product, lamb_pp_scaled[:, :, :, np.newaxis])
    denom_recip = np.reciprocal(denom)
    denom_recip_view = denom_recip.view()
    denom_recip_view.shape = (M, P)

    # Calculate numerator (M x P x 1 x 1))
    left_product_view = left_product.view().reshape(M, P, P + 1)
    numer = np.matmul(left_product_view, betas[:, :, np.newaxis]) # (M, P, P+1) x (M, P+1, 1) = (M, P, 1)
    numer_view = numer.view().reshape(M, P)

    # Calculate result betas and standard errors
    new_betas = denom_recip_view * numer_view
    new_beta_ses = np.sqrt(denom_recip_view)

    return new_betas, new_beta_ses, left_product_view


## run IBSEP
def run_ibsep(betas, ses, props, ld, Omega=None, Sigma=None):
    """
    :param betas: Mx(P+1)
    :param ses: Mx(P+1)
    :param props: P
    :param ld: M
    :return: blue estimator, ses_blue
    """
    P = len(props)
    M = len(ld)
    
    # zscores
    zscores = betas / ses # (M, P+1)
    
    # blue estimator
    betas_blue, ses_blue = TensorIBSEP(betas, ses, props, Omega, Sigma, ld) # (M, P)
    return betas_blue, ses_blue

