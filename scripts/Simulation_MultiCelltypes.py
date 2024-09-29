import os, sys, time, glob, argparse

threads = '1'
os.environ["OMP_NUM_THREADS"] = threads
os.environ["OPENBLAS_NUM_THREADS"] = threads
os.environ["MKL_NUM_THREADS"] = threads
os.environ["VECLIB_MAXIMUM_THREADS"] = threads
os.environ["NUMEXPR_NUM_THREADS"] = threads

from datetime import datetime
from tqdm import tqdm
import numpy as np
import pandas as pd
import scipy.stats as st
from multiprocessing import Pool
from ibsep.core import *
from ibsep import ldsc
from ibsep.utils import ReadPlink, ScaleGenotype
import statsmodels.api as sm

def make_pd_1(Omega, Omega_se, shrink=0.9):
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

def make_pd_2(Omega):
    eigval, eigvec = np.linalg.eig(Omega)
    print(eigval.min())
    K = Omega.shape[-1]
    eigval = np.maximum(eigval, 1e-12)
    Omega = eigvec @ (eigval * np.eye(K)) @ eigvec.T
    # force omega to be real
    Omega = np.real(Omega)
    return Omega

# first 1000 SNPs
def get_genotype_ldscore(snp_num):
    
    ## genotype and ldscore
    plink_file = 'data/UKB_20k_chr20/UKB_20k_hm3_rmAmbi_chr20'
    G = ReadPlink(plink_file)
    bim = pd.read_csv(f'{plink_file}.bim', sep='\t', header=None)
    bim.insert(0, 'bim_index', bim.index.values)
    # add window id
    ld_df = pd.read_csv('data/LDscoresEUR-EAS/ldsc_annot_EUR_EAS_1mb_TGP_hm3_chr20_std_pop1.gz', sep='\t')
    # merge genotype and ldscore
    merge = pd.merge(bim, ld_df, how='inner', left_on=1, right_on='SNP')
    G_merge = G[:, merge.bim_index.values]
    # select one window
    G_win = G_merge[:1100, :snp_num]
    ld = merge['base'].values[:snp_num]
    ld = ld[:, None] 
    print(snp_num, G_win.shape)
    
    return G_win, ld
    
def generate_data(h1sq, h2sq, rg, pi1, Nc, Nt, pcau, delta, r, save_dir, seed, G=None, ieqtl=False):
    
    np.random.seed(int(r * 1000 + seed))
    assert G is not None, 'Please provide G!'
    M = G.shape[1]

        
    ## genotype
    Nc = 100
    Nt = 1000
    Gc = G[:Nc, :]
    Gt = G[Nc:, :]
    Xc = (Gc - np.mean(Gc, axis=0)) / np.std(Gc, axis=0)
    Xt = (Gt - np.mean(Gt, axis=0)) / np.std(Gt, axis=0)

    if metric == 'power':
        ## hyperparameters
        Omega = np.eye(6)*h2sq
        Omega[0,0] = h1sq
        Omega_sub = np.ones(5*5).reshape(5,5)*np.sqrt(h2sq*h2sq)*rg
        Omega_sub[np.diag_indices_from(Omega_sub)] = h2sq
        Omega[1:,1:] = Omega_sub
        Omega[0,1:] = Omega[1:,0] = np.sqrt(h1sq*h2sq)*rg
        Omega = Omega / (pcau * M)
        #Omega = np.array([[h1sq, np.sqrt(h1sq*h2sq)*rg], [np.sqrt(h1sq*h2sq)*rg, h2sq]]) / (pcau * M)
        ## cell grex
        beta1, beta2, beta3, beta4, beta5, beta6 = np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M)
        cau_ids = np.random.choice(np.arange(M), int(pcau*M))
        cau_ids = sorted(cau_ids)
        beta_cau = np.random.multivariate_normal(mean=np.zeros(6), cov=Omega, size=int(pcau*M))
    elif metric == 'type_1_error':
        cutoff = 0.85
        Omega1 = np.eye(6)*h2sq
        Omega1[0,0] = h1sq
        Omega1_sub = np.ones(5*5).reshape(5,5)*np.sqrt(h2sq*h2sq)*rg
        Omega1_sub[np.diag_indices_from(Omega1_sub)] = h2sq
        Omega1[1:,1:] = Omega1_sub
        Omega1 = Omega1 / (pcau * M)

        Omega2 = np.ones(6*6).reshape(6,6)*np.sqrt(h2sq*h2sq)*rg
        Omega2[np.diag_indices_from(Omega2)] = h2sq
        Omega2 = Omega2 / (pcau * M)
        ## cell grex
        causal_snps_num = int(pcau*M)
        causal_snps_num1 = int(causal_snps_num*cutoff)
        beta1, beta2, beta3, beta4, beta5, beta6 = np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M)
        cau_ids = np.random.choice(np.arange(M), causal_snps_num)
        cau_ids = sorted(cau_ids)
        beta_cau1 = np.random.multivariate_normal(mean=np.zeros(6), cov=Omega1, size=causal_snps_num1)
        beta_cau2 = np.random.multivariate_normal(mean=np.zeros(6), cov=Omega2, size=causal_snps_num-causal_snps_num1)
        beta_cau = np.vstack((beta_cau1,beta_cau2))
    
    beta1[cau_ids] = beta_cau[:, 0]
    beta2[cau_ids] = beta_cau[:, 1]
    beta3[cau_ids] = beta_cau[:, 2]
    beta4[cau_ids] = beta_cau[:, 3]
    beta5[cau_ids] = beta_cau[:, 4]
    beta6[cau_ids] = beta_cau[:, 5]
    y1 = Xc @ beta1 + np.sqrt(1 - h1sq) * np.random.randn(Nc)
    y2 = Xc @ beta2 + np.sqrt(1 - h2sq) * np.random.randn(Nc)
    y3 = Xc @ beta3 + np.sqrt(1 - h2sq) * np.random.randn(Nc)
    y4 = Xc @ beta4 + np.sqrt(1 - h2sq) * np.random.randn(Nc)
    y5 = Xc @ beta5 + np.sqrt(1 - h2sq) * np.random.randn(Nc)
    y6 = Xc @ beta6 + np.sqrt(1 - h2sq) * np.random.randn(Nc)
    
    ## tissue grex
    pi1_ind = np.random.beta(pi1*delta, (1-pi1)*delta, Nt)
    pi1s_remain = [_ for _ in pi1s if _ != pi1]
    pi2_ind = np.random.beta(pi1s_remain[0]*delta, (1-pi1s_remain[0])*delta, Nt)
    pi3_ind = np.random.beta(pi1s_remain[1]*delta, (1-pi1s_remain[1])*delta, Nt)
    pi4_ind = np.random.beta(pi1s_remain[2]*delta, (1-pi1s_remain[2])*delta, Nt)
    pi5_ind = np.random.beta(pi1s_remain[3]*delta, (1-pi1s_remain[3])*delta, Nt)
    pi6_ind = np.random.beta(pi1s_remain[4]*delta, (1-pi1s_remain[4])*delta, Nt)

    pis_ind = np.stack((pi1_ind,pi2_ind,pi3_ind,pi4_ind,pi5_ind,pi6_ind),axis=1)
    pis_ind = pis_ind/pis_ind.sum(1)[:,None]

    pi1_ind = pis_ind[:,0]
    pi2_ind = pis_ind[:,1]
    pi3_ind = pis_ind[:,2]
    pi4_ind = pis_ind[:,3]
    pi5_ind = pis_ind[:,4]
    pi6_ind = pis_ind[:,5]
    
    pi1_mean = np.mean(pi1_ind)
    pi2_mean = np.mean(pi2_ind)
    pi3_mean = np.mean(pi3_ind)
    pi4_mean = np.mean(pi4_ind)
    pi5_mean = np.mean(pi5_ind)
    pi6_mean = np.mean(pi6_ind)
    
    yt = pi1_ind * (Xt @ beta1) + pi2_ind * (Xt @ beta2) + \
        pi3_ind * (Xt @ beta3) + pi4_ind * (Xt @ beta4) + \
        pi5_ind * (Xt @ beta5) + pi6_ind * (Xt @ beta6) + \
        np.sqrt(1 - pi1_mean**2*h1sq - pi2_mean**2*h2sq  \
               - pi3_mean**2*h2sq - pi4_mean**2*h2sq - pi5_mean**2*h2sq \
               - pi6_mean**2*h2sq) * np.random.randn(Nt)
    
    ## sumstats
    b1_hat, b2_hat, b3_hat, b4_hat, b5_hat, b6_hat, bt_hat = np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M)
    se1_hat, se2_hat, se3_hat, se4_hat, se5_hat, se6_hat, set_hat = np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M), np.zeros(M)
    if ieqtl:
        bt_int, set_int = np.zeros(M), np.zeros(M)
        zt_int, pvalt_int = np.zeros(M), np.zeros(M)
    for j in range(M):
        b1_hat[j] = np.dot(Xc[:,j], y1) / np.dot(Xc[:,j], Xc[:,j]) 
        b2_hat[j] = np.dot(Xc[:,j], y2) / np.dot(Xc[:,j], Xc[:,j]) 
        b3_hat[j] = np.dot(Xc[:,j], y3) / np.dot(Xc[:,j], Xc[:,j]) 
        b4_hat[j] = np.dot(Xc[:,j], y4) / np.dot(Xc[:,j], Xc[:,j]) 
        b5_hat[j] = np.dot(Xc[:,j], y5) / np.dot(Xc[:,j], Xc[:,j]) 
        b6_hat[j] = np.dot(Xc[:,j], y6) / np.dot(Xc[:,j], Xc[:,j]) 
        bt_hat[j] = np.dot(Xt[:,j], yt) / np.dot(Xt[:,j], Xt[:,j]) 
        e1 = y1 - Xc[:,j] * b1_hat[j]
        e2 = y2 - Xc[:,j] * b2_hat[j]
        e3 = y3 - Xc[:,j] * b3_hat[j]
        e4 = y4 - Xc[:,j] * b4_hat[j]
        e5 = y5 - Xc[:,j] * b5_hat[j]
        e6 = y6 - Xc[:,j] * b6_hat[j]
        et = yt - Xt[:,j] * bt_hat[j]
        se1_hat[j] = np.sqrt(np.dot(e1, e1) / (Nc * np.dot(Xc[:,j], Xc[:,j])))
        se2_hat[j] = np.sqrt(np.dot(e2, e2) / (Nc * np.dot(Xc[:,j], Xc[:,j])))
        se3_hat[j] = np.sqrt(np.dot(e3, e3) / (Nc * np.dot(Xc[:,j], Xc[:,j])))
        se4_hat[j] = np.sqrt(np.dot(e4, e4) / (Nc * np.dot(Xc[:,j], Xc[:,j])))
        se5_hat[j] = np.sqrt(np.dot(e5, e5) / (Nc * np.dot(Xc[:,j], Xc[:,j])))
        se6_hat[j] = np.sqrt(np.dot(e6, e6) / (Nc * np.dot(Xc[:,j], Xc[:,j])))
        set_hat[j] = np.sqrt(np.dot(et, et) / (Nt * np.dot(Xt[:,j], Xt[:,j])))
        if ieqtl:
            X = np.stack([pi1_ind, Xt[:,j], pi1_ind*Xt[:,j]], axis=1)
            model = sm.OLS(yt, X)
            res = model.fit()
            bt_int[j] = res.params[-1]
            set_int[j] = res.bse[-1]
            zt_int[j] = res.tvalues[-1]
            pvalt_int[j] = res.pvalues[-1]
    z1 = b1_hat / se1_hat
    z2 = b2_hat / se2_hat
    z3 = b3_hat / se3_hat
    z4 = b4_hat / se4_hat
    z5 = b5_hat / se5_hat
    z6 = b6_hat / se6_hat
    zt = bt_hat / set_hat
    pval1 = st.norm.sf(abs(z1)) * 2
    pval2 = st.norm.sf(abs(z2)) * 2
    pval3 = st.norm.sf(abs(z3)) * 2
    pval4 = st.norm.sf(abs(z4)) * 2
    pval5 = st.norm.sf(abs(z5)) * 2
    pval6 = st.norm.sf(abs(z6)) * 2
    pvalt = st.norm.sf(abs(zt)) * 2

    ## save data
    # snp info
    snp_data = pd.DataFrame({'SNP': [f'S{str(x).zfill(len(str(M)))}' for x in np.arange(M)],
                             'Beta1': beta1,
                             'Beta2': beta2,
                             'Beta3': beta3,
                             'Beta4': beta4,
                             'Beta5': beta5,
                             'Beta6': beta6})
    # cell gene expression
    cell_grex = pd.DataFrame({'ID': [f'C{str(x).zfill(len(str(Nc)))}' for x in np.arange(Nc)],
                              'y1': y1,
                              'y2': y2,
                              'y3': y3,
                              'y4': y4,
                              'y5': y5,
                              'y6': y6})
    # tissue gene expression and individual pi's
    tissue_grex = pd.DataFrame({'ID': [f'T{str(x).zfill(len(str(Nt)))}' for x in np.arange(Nt)],
                                'yt': yt,
                                'pi1': pi1_ind,
                                'pi2': pi2_ind,
                                'pi3': pi3_ind,
                                'pi4': pi4_ind,
                                'pi5': pi5_ind,
                                'pi6': pi6_ind})
    # cell1 sumstats
    cell1 = pd.DataFrame({'SNP': np.arange(M),
                          'Beta': b1_hat,
                          'SE': se1_hat,
                          'Zscore': z1,
                          'Pval': pval1,
                          'Sig': (pval1 < 0.05).astype('int32'),
                          'Causal': (beta1 != 0).astype('int32')})
    # cell2 sumstats
    cell2 = pd.DataFrame({'SNP': np.arange(M),
                          'Beta': b2_hat,
                          'SE': se2_hat,
                          'Zscore': z2,
                          'Pval': pval2,
                          'Sig': (pval2 < 0.05).astype('int32'),
                          'Causal': (beta2 != 0).astype('int32')})
    # cell3 sumstats
    cell3 = pd.DataFrame({'SNP': np.arange(M),
                          'Beta': b3_hat,
                          'SE': se3_hat,
                          'Zscore': z3,
                          'Pval': pval3,
                          'Sig': (pval3 < 0.05).astype('int32'),
                          'Causal': (beta3 != 0).astype('int32')})
    # cell4 sumstats
    cell4 = pd.DataFrame({'SNP': np.arange(M),
                          'Beta': b4_hat,
                          'SE': se4_hat,
                          'Zscore': z4,
                          'Pval': pval4,
                          'Sig': (pval4 < 0.05).astype('int32'),
                          'Causal': (beta4 != 0).astype('int32')})
    # cell5 sumstats
    cell5 = pd.DataFrame({'SNP': np.arange(M),
                          'Beta': b5_hat,
                          'SE': se5_hat,
                          'Zscore': z5,
                          'Pval': pval5,
                          'Sig': (pval5 < 0.05).astype('int32'),
                          'Causal': (beta5 != 0).astype('int32')})
    # cell6 sumstats
    cell6 = pd.DataFrame({'SNP': np.arange(M),
                          'Beta': b6_hat,
                          'SE': se6_hat,
                          'Zscore': z6,
                          'Pval': pval6,
                          'Sig': (pval6 < 0.05).astype('int32'),
                          'Causal': (beta6 != 0).astype('int32')})
    # tissue sumstats
    tissue = pd.DataFrame({'SNP': np.arange(M),
                           'Beta': bt_hat,
                           'SE': set_hat,
                           'Zscore': zt,
                           'Pval': pvalt,
                           'Sig': (pvalt < 0.05).astype('int32')})
    if ieqtl:
        # ieqtl sumstats
        tissue_ieqtl = pd.DataFrame({'SNP': np.arange(M),
                                     'Beta': bt_int,
                                     'SE': set_int,
                                     'Zscore': zt_int,
                                     'Pval': pvalt_int,
                                     'Sig': (pvalt_int < 0.05).astype('int32')})

    snp_data.to_csv(os.path.join(save_dir, f'snp-data-{r}.csv'), sep='\t', index=False)
    cell_grex.to_csv(os.path.join(save_dir, f'cell-grex-{r}.csv'), sep='\t', index=False)
    tissue_grex.to_csv(os.path.join(save_dir, f'tissue-grex-{r}.csv'), sep='\t', index=False)
    cell1.to_csv(os.path.join(save_dir, f'cell1-sumstats-{r}.csv'), sep='\t', index=False)
    cell2.to_csv(os.path.join(save_dir, f'cell2-sumstats-{r}.csv'), sep='\t', index=False)
    cell3.to_csv(os.path.join(save_dir, f'cell3-sumstats-{r}.csv'), sep='\t', index=False)
    cell4.to_csv(os.path.join(save_dir, f'cell4-sumstats-{r}.csv'), sep='\t', index=False)
    cell5.to_csv(os.path.join(save_dir, f'cell5-sumstats-{r}.csv'), sep='\t', index=False)
    cell6.to_csv(os.path.join(save_dir, f'cell6-sumstats-{r}.csv'), sep='\t', index=False)
    tissue.to_csv(os.path.join(save_dir, f'tissue-sumstats-{r}.csv'), sep='\t', index=False)
    if ieqtl:
        tissue_ieqtl.to_csv(os.path.join(save_dir, f'tissue-ieqtl-sumstats-{r}.csv'), sep='\t', index=False)

def simu(r, save_dir, true_Omega, omega_shrink, max_hsq, method=['blue', 'tram'], zero_corr=False, trun_corr=False):
   
    assert not(zero_corr & trun_corr), 'Not allowed setting zero_corr=True and trun_corr=True at the same time!'
        
    ## true Omega or ldsc estimated Omega
    strOmega = 'trueOmega' if true_Omega else 'ldscOmega'

    cell_types = ['cell1','cell2','cell3','cell4','cell5','cell6']
    K = len(cell_types)

    sumstats = {}
    for cell_type in cell_types:
        sumstat = pd.read_csv(os.path.join(save_dir, f'{cell_type}-sumstats-{r}.csv'), sep='\t')
        sumstat['Z'] = sumstat.Beta / sumstat.SE
        sumstat['N'] = Nc
        sumstat['PVAL'] = sumstat['Pval'] 
        sumstats[cell_type] = sumstat
    tissue = pd.read_csv(os.path.join(save_dir, f'tissue-sumstats-{r}.csv'), sep='\t')
    tissue_grex = pd.read_csv(os.path.join(save_dir, f'tissue-grex-{r}.csv'), sep='\t')
    tissue['Z'] = tissue.Beta / sumstat.SE
    tissue['N'] = Nt
    tissue['PVAL'] = tissue['Pval'] 
    M = len(ld)
    
    ## inputs
    zs_ = {k: v.Z.values for k, v in sumstats.items()} | {'tissue': tissue.Z.values}
    betas_ = {k: v.Beta.values for k, v in sumstats.items()} | {'tissue': tissue.Beta.values}
    ses_ = {k: v.SE.values for k, v in sumstats.items()} | {'tissue': tissue.SE.values}
    zs = np.stack(list(zs_.values()), axis=1)
    betas = np.stack(list(betas_.values()), axis=1)
    ses = np.stack(list(ses_.values()), axis=1)
    # print(f'Data shapes: zscores: {zs.shape}, betas: {betas.shape}, ses: {ses.shape}')

    pi1 = np.mean(tissue_grex['pi1'].values)
    pi2 = np.mean(tissue_grex['pi2'].values)
    pi3 = np.mean(tissue_grex['pi3'].values)
    pi4 = np.mean(tissue_grex['pi4'].values)
    pi5 = np.mean(tissue_grex['pi5'].values)
    pi6 = np.mean(tissue_grex['pi6'].values)
    
    props = [pi1, pi2, pi3, pi4, pi5, pi6]

    ## Omega and Sigma
    assert method in ['blue', 'tram'], "method should be 'blue' or 'tram'!"
    # fix Sigma to be identity matrix
    Sigma = np.eye(K+1) if method == 'blue' else np.eye(K)
    if true_Omega:
        h1sq = float(save_dir.split('/')[-1].split('-')[0][5:])
        h2sq = float(save_dir.split('/')[-1].split('-')[1][5:])  
        rg = float(save_dir.split('/')[-1].split('-')[2][3:])
        Omega = np.eye(K)*h2sq
        Omega[0,0] = h1sq
        Omega = Omega / M    
        np.save(os.path.join(save_dir, f'Omega-{method}-{strOmega}-{r}.npy'), Omega)
    else:
        Omega, _, Omega_se = ldsc.run_ldscore_regressions_tram_multi(sumstats,ld,tissue,intercept=np.eye(K+1))
        
        ## numerical operations
        Omega = Omega.squeeze()[:K,:K]
        Omega_se = Omega_se.squeeze()[:K,:K]
        oriOmega = Omega.copy() * M
        oriOmega_se = Omega_se.copy() * M
        for i in range(len(Omega)):
            Omega[i, i] = min(Omega[i, i], max_hsq/M) # max heritability: 0.8
        Omega = make_positive_definite(Omega, 0.95) # max correlation: 0.95
    
        ## hij test
        if trun_corr:
            n0 = int(K * (K - 1) / 2)
            n = 0
            for i in range(0, K):
                for j in range(i+1, K):
                    hij_z = Omega[i, j] / Omega_se[i, j]
                    hij_p = st.norm.sf(abs(hij_z)) * 2
                    if method == 'blue':
                        if hij_p >= 0.01:
                            Omega[i, j] = 1e-32
                            Omega[j, i] = 1e-32
                            n += 1
                    else:
                        if hij_p >= 0.05:
                            Omega[i, j] = 1e-32
                            Omega[j, i] = 1e-32
                            n += 1
        if zero_corr:
            for i in range(0, K):
                for j in range(i+1, K):
                    Omega[i, j] = 1e-32
                    Omega[j, i] = 1e-32

        ## make pd
        filename = f'Omega-{method}-{strOmega}-{r}.npz'
        if not np.all(np.linalg.eigvals(Omega) > 0):
            if omega_shrink < 1:
                Omega, Omega_se, n_shrink = make_pd_1(Omega, Omega_se, omega_shrink)
                shrinkOmega = Omega.copy() * M
                shrinkOmega_se = Omega_se.copy() * M
                np.savez(os.path.join(out_dir, filename), oriOmega=oriOmega, oriOmega_se=oriOmega_se, 
                         shrinkOmega=shrinkOmega, shrinkOmega_se=shrinkOmega_se,  
                         omega_shrink=omega_shrink, n_shrink=n_shrink)
            else:
                Omega = make_pd_2(Omega)
                pdOmega = Omega.copy() * M
                np.savez(os.path.join(out_dir, filename), oriOmega=oriOmega, oriOmega_se=oriOmega_se, 
                         pdOmega=pdOmega)
        else:
            np.savez(os.path.join(out_dir, filename), oriOmega=oriOmega, oriOmega_se=oriOmega_se)
        
    ## calculate ours or tram's blue estimators
    if method == 'blue':
        betas_tensor, ses_tensor = run_ours(betas, ses, props, ld, Omega=Omega, Sigma=Sigma)
    elif method == 'tram':
        betas_tensor, ses_tensor = run_tram(betas[:,:K], ses[:,:K], ld, Omega=Omega, Sigma=Sigma)
        
    zs_tensor = betas_tensor / ses_tensor
    pval_tensor = st.norm.sf(abs(zs_tensor)) * 2
    cau_tensor = (pval_tensor < 0.05).astype('int32')
    
    ## save results
    for i, cell_type in enumerate(cell_types):

        celli_meta = pd.DataFrame({'SNP': np.arange(M),
                               'Beta': betas_tensor[:, i],
                               'SE': ses_tensor[:, i],
                               'Zscore': zs_tensor[:, i],
                               'Pval': pval_tensor[:, i],
                               'Sig': cau_tensor[:, i]})
        
        filename = f'{cell_type}-{method}-{strOmega}-{r}.csv'
        if trun_corr:
            filename = '-'.join(filename.split('-')[:-1] + ['truncorr', filename.split('-')[-1]])
        if zero_corr:
            filename = '-'.join(filename.split('-')[:-1] + ['zerocorr', filename.split('-')[-1]])
        celli_meta.to_csv(os.path.join(save_dir, filename), sep='\t', index=False)

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--out-root', type=str, help='root directory for saving simulation data and results', required=True)
    parser.add_argument('--metric', type=str, choices=['type_1_error', 'power'], help='type_1_error or power experiment', required=True)
    parser.add_argument('--method', type=str, choices=['blue', 'tram'], help='use ours (blue) or tram', required=True)
    parser.add_argument('--snps', type=int, help='SNPs num', default=1000, required=True)
    parser.add_argument('--delta', type=int, help='delta', default=5)
    parser.add_argument('--true-Omega', action='store_true', help='use true Omega')
    parser.add_argument('--parallel', action='store_true', help='run replications parallelly')
    parser.add_argument('--omega-shrink', type=float, help='omega shrinkage factor', default=0.9)
    parser.add_argument('--max-hsq', type=float, help='omega shrinkage factor', default=0.8)
    parser.add_argument('--zero-corr', action='store_true', help='force genetic correlation to be zero (h12=0)')
    parser.add_argument('--trun-corr', action='store_true', help='truncate insignificant h12 to be zero (h12_insig=0)')
    parser.add_argument('--ieqtl', action='store_true', help='perform ieqtl')
    args = parser.parse_args()
    
    ## simulation setting
    seed = 0 
    Nc = 100
    Nt = 1000
    h2sqs = [0.01, 0.05, 0.1]
    pi1s = [0.5, 0.24, 0.12, 0.08, 0.04, 0.02] # mean proportion of interested cell type 
    pcaus = [0.005,0.02] 
    delta = args.delta
    nrep = 100
    start = 0
    reps = np.arange(start+1, start+nrep+1)
    metric = args.metric
    true_Omega = args.true_Omega
    method = args.method 
    parallel = args.parallel
    zero_corr = args.zero_corr
    trun_corr = args.trun_corr
    out_root = args.out_root
    snp_num = args.snps
    ieqtl = args.ieqtl
    omega_shrink = args.omega_shrink
    max_hsq = args.max_hsq
    
    # different hyperparameters
    strgeno = 'real_geno'
    if metric == 'type_1_error':
        h1sqs = [0.0]
        rgs = [0.0,0.3,0.6,0.9]
        out_root = os.path.join(out_root, strgeno, 'type_1_error')
    elif metric == 'power':
        h1sqs = [0.05]
        h2sqs = [0.05]
        # h1sqs = [0.5] # debug
        # h2sqs = [0.5] # debug
        rgs = [0.0, 0.3, 0.6, 0.9]
        out_root = os.path.join(out_root, strgeno, 'power')
    
    ## genotype and ldscore
    # use window 0
    snp_nums = [snp_num]
    for snp_num in snp_nums:
        G, ld = get_genotype_ldscore(snp_num)

        ## generate data and save
        t0 = time.time()
        print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] ***** Start generating data (window {snp_num}) *****")
        for h1sq in h1sqs:
            for h2sq in h2sqs:
                for rg in rgs:
                    for pi1 in pi1s:
                        for pcau in pcaus:
                            out_dir = os.path.join(out_root, f'h1sq={h1sq}-h2sq={h2sq}-rg={rg}-pi1={pi1}-win{snp_num}-Nc{Nc}-Nt{Nt}-pcau{pcau}-seed{seed}')
                            os.makedirs(out_dir, exist_ok=True)
                            stream = tqdm(reps, desc=f'win={snp_num}, h1sq={h1sq}, h2sq={h2sq}, rg={rg}, pi1={pi1}, pcau={pcau}', ascii=' >=')
                            for r in stream:
                                if os.path.exists(os.path.join(out_dir, f'tissue-sumstats-{r}.csv')):
                                    continue
                                else:
                                    generate_data(h1sq, h2sq, rg, pi1, Nc, Nt, pcau, delta, r, out_dir, seed, G, ieqtl)
        print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] ***** Finish generating data (window {snp_num}) *****")
        print(f"Generating data time (unparallel): {time.time() - t0:.2f}s")

        ## perform ours/tram
        if parallel:
            t0 = time.time()
            print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] ***** Start perform method: {method} (window {snp_num}) *****")
            for h1sq in h1sqs:
                for h2sq in h2sqs:
                    for rg in rgs:
                        for pi1 in pi1s:
                            for pcau in pcaus:
                                out_dir = os.path.join(out_root, f'h1sq={h1sq}-h2sq={h2sq}-rg={rg}-pi1={pi1}-win{snp_num}-Nc{Nc}-Nt{Nt}-pcau{pcau}-seed{seed}')
                                inp_args = []
                                for r in reps:
                                    inp_args.append((r, out_dir, true_Omega, omega_shrink, max_hsq, method, zero_corr, trun_corr))
                                try:
                                    with Pool(len(reps)) as pool:
                                        pool.starmap(simu, inp_args)
                                        pool.close()
                                        pool.join()
                                except Exception as e:
                                    print(e)
            print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] ***** Finish performing method: {method} (window {snp_num}) *****")
            print(f"Performing method time (parallel): {time.time() - t0:.2f}s")
        else:
            t0 = time.time()
            print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] ***** Start perform method: {method} (window {snp_num}) *****")
            for h1sq in h1sqs:
                for h2sq in h2sqs:
                    for rg in rgs:
                        for pi1 in pi1s:
                            for pcau in pcaus:
                                out_dir = os.path.join(out_root, f'h1sq={h1sq}-h2sq={h2sq}-rg={rg}-pi1={pi1}-win{snp_num}-Nc{Nc}-Nt{Nt}-pcau{pcau}-seed{seed}')
                                stream = tqdm(reps, desc=f'win={snp_num}, h1sq={h1sq}, h2sq={h2sq}, rg={rg}, pi1={pi1}, pcau={pcau}', ascii=' >=')
                                for r in stream:
                                    simu(r, out_dir, true_Omega, omega_shrink, max_hsq, method, zero_corr, trun_corr)
            print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] ***** Finish performing method: {method} (window {snp_num}) *****")
            print(f"Performing method time (unparallel): {time.time() - t0:.2f}s")
