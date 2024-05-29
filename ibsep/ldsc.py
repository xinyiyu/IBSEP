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

def get_coef_raw(ldscore_mat, zsc1, zsc2, max_int):
    """
    Obtain coefficient for heritability / genetic covariance
    """
   
    # get dimension
    nsnp,ncoef = ldscore_mat.shape[0],ldscore_mat.shape[1]

    # set up variables
    nprod = np.sqrt((zsc1['N'].values)*(zsc2['N'].values))
    zprod = (zsc1['Z'].values)*(zsc2['Z'].values)
    score_sum = np.sum(ldscore_mat[:,:-1], axis=1)
    
    # get averages
    mean_ldscore = np.mean(score_sum)
    mean_zprod = np.mean(zprod)
    mean_nprod = np.mean(nprod)

    # estimate intercept
    zprod_sorted = np.sort(zprod)
    idx = int(float(nsnp)*0.95)
    intercept = np.mean(zprod_sorted[0:idx])
    if intercept > max_int:
        intercept = max_int

    # get raw coef
    coef = (mean_zprod-intercept) / (mean_nprod*mean_ldscore)

    return coef, intercept


def get_pred(coef, x, n1, n2, intercept):
    
    nprod = np.sqrt(n1*n2)
    score_sum = np.sum(x[:,:-1], axis=1)
    pred = coef*score_sum*nprod + intercept

    return pred.astype(np.float32)


def update_weights(w1, w2, wx, pred1, pred2, predx):
   
    var1 = 2.0*np.square(pred1)
    var2 = 2.0*np.square(pred2)
    var12 = np.sqrt(var1*var2)/2.0 + np.square(predx)

    new_w1 = w1 / var1
    new_w2 = w2 / var2
    new_wx = wx / var12

    return new_w1, new_w2, new_wx


def regression(x, y, constrain_intercept,subtract):
    """
    Perform least square regression with block jack knife
    """

    # perform regression
    xtx = np.dot(x.T, x)
    xty = np.dot(x.T, y)
    ncoef = xtx.shape[0]

    # obtain coefficient
    n = len(y)
    coef = np.zeros(ncoef, dtype=np.float32)
    if constrain_intercept:
        coef[:-1] = np.linalg.solve(xtx[:-1,:-1], xty[:-1])
        coef[-1] = subtract        
        k = x.shape[1] - 1
    else:
        coef = np.linalg.solve(xtx, xty)
        k = x.shape[1]
        
    ## calculate standard errors
    # estimated residual variance 
    err = y - x @ coef
    sigsq = np.dot(err, err) / (n - k)
    # covariance matrix of betas
    cov_beta = sigsq * np.linalg.inv(xtx)
    coef_se = np.sqrt(np.diag(cov_beta))

    return coef, coef_se


def get_coef(score, zsc1, zsc2, w, constrain_intercept, subtract):
    """
    Obtain coefficient for heritability / genetic covariance
    """
    
    # set up the regression
    nsnp = score.shape[0]
    nprod = np.sqrt((zsc1['N'].values)*(zsc2['N'].values))
    zprod = (zsc1['Z'].values)*(zsc2['Z'].values)
    if constrain_intercept:
        zprod -= subtract
    score[:,:-1] = score[:,:-1] * nprod[:,np.newaxis]
    
    # scale the matrix to improve matrix condition
    nbar = np.mean(nprod)
    score[:,:-1] /= nbar
    
    # apply the weight
    score_w = score * np.sqrt(w[:,np.newaxis])
    zprod_w = zprod * np.sqrt(w)
    
    # run regression
    coef, coef_se = regression(score_w, zprod_w, constrain_intercept, subtract)
    coef[:-1] /= nbar
    coef_se[:-1] /= nbar
    
    return coef, coef_se
  

def estimate_h2(sumstat1,ldscore1,reg_w1=1,constrain_intercept=False,int1=1):
    reg_w1[reg_w1<0] = 0.01
    # tau: [coeficient, intercept]
    ldscore1 = np.hstack((ldscore1,np.ones((ldscore1.shape[0],1))))
    tau1,int1_ = get_coef_raw(ldscore1,sumstat1,sumstat1,1)
    n1 = sumstat1['N'].values
    pred1 = get_pred(tau1, ldscore1, n1, n1, int1_)
    var1 = 2.0*np.square(pred1)
    weight1_ = reg_w1 / var1
    tau1, tau1_se = get_coef(ldscore1, sumstat1, sumstat1,
       weight1_, constrain_intercept, int1)
    return tau1, tau1_se


def estimate_gc(sumstat1,sumstat2,ldscore1,ldscore2,ldscorex,reg_w1=1,reg_w2=1,constrain_intercept=False,int1=1,int2=1,intx=0):
    reg_w1[reg_w1<0] = 0.01
    reg_w2[reg_w2<0] = 0.01
    reg_wx = np.sqrt(reg_w1*reg_w2)

    ldscore1 = np.hstack((ldscore1,np.ones((ldscore1.shape[0],1))))
    ldscore2 = np.hstack((ldscore2,np.ones((ldscore2.shape[0],1))))
    ldscorex = np.hstack((ldscorex,np.ones((ldscorex.shape[0],1))))

    # get initial estimate of coefs for constructing weights
    tau1,int1_ = get_coef_raw(ldscore1,sumstat1,sumstat1,1)
    tau2,int2_ = get_coef_raw(ldscore2,sumstat2,sumstat2,1)
    theta,intx_ = get_coef_raw(ldscorex,sumstat1,sumstat2,0)

    # predict chi-square and prodcut 
    n1,n2 = sumstat1['N'].values,sumstat2['N'].values
    pred1 = get_pred(tau1, ldscore1, n1, n1, int1_)
    pred2 = get_pred(tau2, ldscore2, n2, n2, int2_)
    predx = get_pred(theta, ldscorex, n1, n2, intx_)
    
    # update weight
    weight1_, weight2_, weightx_ = update_weights(reg_w1, reg_w2, reg_wx, pred1, pred2, predx)

    # get regression coefficients
    tau1, tau1_se = get_coef(ldscore1, sumstat1, sumstat1,
       weight1_, constrain_intercept, int1)
    tau2, tau2_se = get_coef(ldscore2, sumstat2, sumstat2,
        weight2_, constrain_intercept, int2)
    theta, theta_se = get_coef(ldscorex, sumstat1, sumstat2,
        weightx_, constrain_intercept, intx)

    return tau1,tau2,theta,tau1_se,tau2_se,theta_se


def ldscore_regression_h2(idx1,df_eas,ldscore,intercept=None):
    if intercept is None:
        #idx1 = (df_eas['Z']**2<30).values
        fit_step1, fit_step1_se = estimate_h2(df_eas.loc[idx1],ldscore[idx1],
                reg_w1 = 1/ldscore[idx1,0],constrain_intercept=False)
        fit_step2, fit_step2_se = estimate_h2(df_eas,ldscore,
                reg_w1 = 1/ldscore[:,0],constrain_intercept=True,int1=fit_step1[-1])
    else:
        fit_step2, fit_step2_se = estimate_h2(df_eas,ldscore,
                reg_w1 = 1/ldscore[:,0],constrain_intercept=True,int1=intercept)
    fit_step2[0] = fit_step2[0] if fit_step2[0]>0 else 1e-12 

    return fit_step2, fit_step2_se


def ldscore_regression_gc(idx1,df_eas,df_eur,ldscore1,ldscore2,ldscorete,intercept1=None,intercept2=None,interceptx=None):
    #idx1 = (df_eas['Z']**2<30).values & (df_eur['Z']**2<30).values
    fit_step1 = estimate_gc(df_eas.loc[idx1],df_eur.loc[idx1],ldscore1[idx1],ldscore2[idx1],ldscorete[idx1],
            reg_w1 = 1/ldscore1[idx1,0],reg_w2 = 1/ldscore2[idx1,0],constrain_intercept=False)
    if intercept1 is None:
        intercept1 = fit_step1[0][-1]
    if intercept2 is None:
        intercept2 = fit_step1[1][-1]
    if interceptx is None:
        interceptx = fit_step1[2][-1]
    fit_step2 = estimate_gc(df_eas,df_eur,ldscore1,ldscore2,ldscorete,
            reg_w1 = 1/ldscore1[:,0],reg_w2 = 1/ldscore2[:,0],constrain_intercept=True,
            int1=intercept1,int2=intercept2,intx=interceptx)
    
    h1_snp = fit_step2[0][0]
    h2_snp = fit_step2[1][0]
    
    h1_snp = h1_snp if h1_snp>0 else 1e-12 
    h2_snp = h2_snp if h2_snp>0 else 1e-12 
    gc_avg = fit_step2[2][0]/np.sqrt(h1_snp*h2_snp)
    gc_avg = gc_avg if gc_avg<1 else 0.99
    gc_avg = gc_avg if gc_avg>-1 else -0.99   
    h12_snp = np.sqrt(h1_snp*h2_snp)*gc_avg
    fit_step2[2][0] = h12_snp
    
    return fit_step2[2], fit_step2[-1]    

## estimate omega and sigma by ldsc
def run_ldscore_regressions_multi(sumstats,ldscore,tissue,intercept=None):

    K = len(sumstats)
    names = list(sumstats.keys())
    p_thres = st.norm.sf(30**.5/2) # corresponding to Z**2=30
    tmp = []
    for k, sumstat in sumstats.items():
        tmp.append(list(sumstat.PVAL > p_thres))
    if intercept is None:
        tmp.append(list(tissue.PVAL > p_thres))
    idx1 = np.array(tmp).all(axis=0)
    # print(f'cell idx1: {sum(idx1)}/{len(idx1)}')
    result_coefs = np.zeros((2, K+1, K+1))
    result_coefs_se = np.zeros((2, K+1, K+1))
    if intercept is not None:
        if intercept.shape[0] != intercept.shape[1] or intercept.shape[0] != K+1:
            raise ValueError('user definded intercept must be a square matrix with dimention = number of sumstats')

    ## Cell
    # Calculate diagonal coefficients first
    for i, (k, sumstat) in enumerate(sumstats.items()):
        if intercept is None:
            result_coefs[:,i,i], result_coefs_se[:,i,i] = ldscore_regression_h2(idx1,sumstat,ldscore)
        else:
            result_coefs[:,i,i], result_coefs_se[:,i,i] = ldscore_regression_h2(idx1,sumstat,ldscore,intercept[i,i])

    # Calculate each off-diagonal element
    for i in range(K):
        ki = names[i]
        for j in range(i+1, K):
            kj = names[j]
            if intercept is None:
                result_coefs[:,i,j], result_coefs_se[:,i,j] = ldscore_regression_gc(idx1,sumstats[ki],sumstats[kj],ldscore,ldscore,ldscore)
            else:
                result_coefs[:,i,j], result_coefs_se[:,i,j] = ldscore_regression_gc(idx1,sumstats[ki],sumstats[kj],ldscore,ldscore,ldscore,\
                        intercept[i,i],intercept[j,j],intercept[i,j])
            result_coefs[:,j,i] = result_coefs[:,i,j] 
            result_coefs_se[:,j,i] = result_coefs_se[:,i,j]

    ## Tissue
    # Diagonal
    if intercept is None:
        result_coefs[:,K,K], result_coefs_se[:,K,K] = ldscore_regression_h2(idx1,sumstat,ldscore)

    # Off-diagonal
    if intercept is None:
        for i in range(K):
            ki = names[i]
            result_coefs[:,i,K], result_coefs_se[:,i,K] = ldscore_regression_gc(idx1,sumstats[ki],tissue,ldscore,ldscore,ldscore)
            result_coefs[:,K,i] = result_coefs[:,i,K]
            result_coefs_se[:,K,i] = result_coefs_se[:,i,K]

    return result_coefs[:-1], result_coefs[-1], result_coefs_se[:-1]  

