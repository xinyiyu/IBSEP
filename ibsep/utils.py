import os, sys, copy
import pandas as pd
import numpy as np
from pandas_plink import read_plink1_bin
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

def ReadPlink(plink_file, dtype=np.float32):
    Genotype = read_plink1_bin(plink_file+".bed", plink_file+".bim", plink_file+".fam", verbose=False)
    Genotype = Genotype.astype(np.int8)
    G_geno = Genotype.values
    G_geno[np.isnan(G_geno)] = 2
    G_geno = 2-G_geno
    return G_geno.astype(dtype)

def ScaleGenotype(G_geno,dtype=np.float32):
    G_geno_mean = G_geno.mean(axis=0)
    G_geno_sd = G_geno.var(axis=0)**.5
    G_geno = (G_geno-G_geno_mean)/G_geno_sd
    return G_geno.astype(dtype)

def GenotypeMoment(G_geno):
    G_geno_mean = G_geno.mean(axis=0)
    G_geno_sd = G_geno.var(axis=0)**.5
    return G_geno_mean, G_geno_sd

def restric_corr(A,maxc=0.95):
    d = A.shape[0]
    for i in range(d):
        A[i,i] = A[i,i] if A[i,i]>0 else 1e-12 
    for i in range(d):
        for j in range(i+1,d):
            if A[i,i]<=1e-12 or A[j,j]<=1e-12:
                A[i,j] = A[j,i] = 0
                continue
            if A[i,j] > maxc*(A[i,i]*A[j,j])**.5:
                A[i,j] = A[j,i] = maxc*(A[i,i]*A[j,j])**.5
            if A[i,j] < -maxc*(A[i,i]*A[j,j])**.5:
                A[i,j] = A[j,i] = -maxc*(A[i,i]*A[j,j])**.5
    return A

def make_omega_pd(Omega, Omega_se, shrink=0.9):
    Omega = Omega.squeeze()
    K = Omega.shape[-1]
    eigs = np.linalg.eigvals(Omega)
    diag_Omega = np.diag(Omega) * np.eye(K)
    off_Omega = Omega - diag_Omega
    diag_Omega_se = np.diag(Omega_se) * np.eye(K)
    off_Omega_se = Omega_se - diag_Omega_se
    j = 0
    while not np.all(eigs > 0):
        off_Omega = off_Omega * shrink
        off_Omega_se = off_Omega_se * shrink
        Omega = diag_Omega + off_Omega
        Omega_se = diag_Omega_se + off_Omega_se
        eigs = np.linalg.eigvals(Omega)
        j += 1
    print(j)
    return Omega, Omega_se, j
