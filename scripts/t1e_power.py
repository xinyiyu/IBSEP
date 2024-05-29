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
    
def get_genotype_ldscore(win_id):
    
    ## genotype and ldscore
    plink_file = 'data/UKB_20k_chr20/UKB_20k_hm3_rmAmbi_chr20'
    G = ReadPlink(plink_file)
    bim = pd.read_csv(f'{plink_file}.bim', sep='\t', header=None)
    bim.insert(0, 'bim_index', bim.index.values)
    # add window id
    ld_df = pd.read_csv('data/LDscoresEUR-EAS/ldsc_annot_EUR_EAS_1mb_TGP_hm3_chr20_std_pop1.gz', sep='\t')
    winnames = ld_df.columns[4:].values.tolist()
    starts = [eval(x.split('_')[1]) for x in winnames]
    ends = [eval(x.split('_')[2]) for x in winnames]
    windows = []
    for i, row in ld_df.iterrows():
        bp = row.BP
        for win_id_, (s, e) in enumerate(zip(starts, ends)):
            if bp >= s and bp < e:
                windows.append(win_id_)
                break
    ld_df.insert(3, 'window', windows)
    # merge genotype and ldscore
    merge = pd.merge(bim, ld_df, how='inner', left_on=1, right_on='SNP')
    G_merge = G[:, merge.bim_index.values]
    # select one window
    win_name = winnames[win_id]
    G_win = G_merge[:1100, np.where(merge.window == win_id)[0]]
    ld = merge.loc[merge.window == win_id, win_name].values
    ld = ld[:, None] 
    print(win_id, G_win.shape)
    
    return G_win, ld
    
def generate_data(h1sq, h2sq, rg, pi1, Nc, Nt, pcau, r, save_dir, seed, G=None, ieqtl=False):
    
    np.random.seed(int(r * 1000 + seed))
    assert G is not None, 'Please provide G!'
    M = G.shape[1]

    ## hyperparameters
    Omega = np.array([[h1sq, np.sqrt(h1sq*h2sq)*rg], [np.sqrt(h1sq*h2sq)*rg, h2sq]]) / (pcau * M)
    
    ## genotype
    Nc = 100
    Nt = 1000
    Gc = G[:Nc, :]
    Gt = G[Nc:, :]
    Xc = (Gc - np.mean(Gc, axis=0)) / np.std(Gc, axis=0)
    Xt = (Gt - np.mean(Gt, axis=0)) / np.std(Gt, axis=0)

    ## cell grex
    beta1, beta2 = np.zeros(M), np.zeros(M)
    cau_ids = np.random.choice(np.arange(M), int(pcau*M))
    beta_cau = np.random.multivariate_normal(mean=np.zeros(2), cov=Omega, size=int(pcau*M))
    beta1[cau_ids] = beta_cau[:, 0]
    beta2[cau_ids] = beta_cau[:, 1]
    y1 = Xc @ beta1 + np.sqrt(1 - h1sq) * np.random.randn(Nc)
    y2 = Xc @ beta2 + np.sqrt(1 - h2sq) * np.random.randn(Nc)
    
    ## tissue grex
    delta = 50
    pi1_ind = np.random.beta(pi1*delta, (1-pi1)*delta, Nt)
    pi2_ind = 1 - pi1_ind
    pi1_mean = np.mean(pi1_ind)
    pi2_mean = np.mean(pi2_ind)
    yt = pi1_ind * (Xt @ beta1) + pi2_ind * (Xt @ beta2) + np.sqrt(1 - pi1_mean**2*h1sq - pi2_mean**2*h2sq - 2*pi1_mean*pi2_mean*rg*np.sqrt(h1sq*h2sq)) * np.random.randn(Nt)
    
    ## sumstats
    b1_hat, b2_hat, bt_hat = np.zeros(M), np.zeros(M), np.zeros(M)
    se1_hat, se2_hat, set_hat = np.zeros(M), np.zeros(M), np.zeros(M)
    if ieqtl:
        bt_int, set_int = np.zeros(M), np.zeros(M)
        zt_int, pvalt_int = np.zeros(M), np.zeros(M)
    for j in range(M):
        b1_hat[j] = np.dot(Xc[:,j], y1) / np.dot(Xc[:,j], Xc[:,j]) 
        b2_hat[j] = np.dot(Xc[:,j], y2) / np.dot(Xc[:,j], Xc[:,j]) 
        bt_hat[j] = np.dot(Xt[:,j], yt) / np.dot(Xt[:,j], Xt[:,j]) 
        e1 = y1 - Xc[:,j] * b1_hat[j]
        e2 = y2 - Xc[:,j] * b2_hat[j]
        et = yt - Xt[:,j] * bt_hat[j]
        se1_hat[j] = np.sqrt(np.dot(e1, e1) / (Nc * np.dot(Xc[:,j], Xc[:,j])))
        se2_hat[j] = np.sqrt(np.dot(e2, e2) / (Nc * np.dot(Xc[:,j], Xc[:,j])))
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
    zt = bt_hat / set_hat
    zs = np.stack([z1, z2], axis=1)
    pval1 = st.norm.sf(abs(z1)) * 2
    pval2 = st.norm.sf(abs(z2)) * 2
    pvalt = st.norm.sf(abs(zt)) * 2

    ## save data
    # snp info
    snp_data = pd.DataFrame({'SNP': [f'S{str(x).zfill(len(str(M)))}' for x in np.arange(M)],
                             'Beta1': beta1,
                             'Beta2': beta2})
    # cell gene expression
    cell_grex = pd.DataFrame({'ID': [f'C{str(x).zfill(len(str(Nc)))}' for x in np.arange(Nc)],
                              'y1': y1,
                              'y2': y2})
    # tissue gene expression and individual pi's
    tissue_grex = pd.DataFrame({'ID': [f'T{str(x).zfill(len(str(Nt)))}' for x in np.arange(Nt)],
                                'yt': yt,
                                'pi1': pi1_ind,
                                'pi2': pi2_ind})
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
    tissue.to_csv(os.path.join(save_dir, f'tissue-sumstats-{r}.csv'), sep='\t', index=False)
    if ieqtl:
        tissue_ieqtl.to_csv(os.path.join(save_dir, f'tissue-ieqtl-sumstats-{r}.csv'), sep='\t', index=False)

def simu(r, save_dir, true_Omega, omega_shrink, max_hsq, method=['blue', 'tram'], zero_corr=False, trun_corr=False):

    K = 2
    
    assert not(zero_corr & trun_corr), 'Not allowed setting zero_corr=True and trun_corr=True at the same time!'
        
    ## true Omega or ldsc estimated Omega
    strOmega = 'trueOmega' if true_Omega else 'ldscOmega'
    
    ## load data: need sumstats, pi's
    cell1 = pd.read_csv(os.path.join(save_dir, f'cell1-sumstats-{r}.csv'), sep='\t')
    cell2 = pd.read_csv(os.path.join(save_dir, f'cell2-sumstats-{r}.csv'), sep='\t')
    tissue = pd.read_csv(os.path.join(save_dir, f'tissue-sumstats-{r}.csv'), sep='\t')
    tissue_grex = pd.read_csv(os.path.join(save_dir, f'tissue-grex-{r}.csv'), sep='\t')
    
    ## inputs to ldsc and our method
    zs = np.stack([cell1['Zscore'].values, cell2['Zscore'].values], axis=1)
    betas = np.stack([cell1['Beta'].values, cell2['Beta'].values, tissue['Beta'].values], axis=1)
    ses = np.stack([cell1['SE'].values, cell2['SE'].values, tissue['SE'].values], axis=1)
    pi1_ind = tissue_grex['pi1'].values
    pi1 = np.mean(pi1_ind)
    props = [pi1, 1 - pi1]
    M = len(zs)

    ## Omega and Sigma
    assert method in ['blue', 'tram'], "method should be 'blue' or 'tram'!"
    # fix Sigma to be identity matrix
    Sigma = np.eye(3) if method == 'blue' else np.eye(2)
    if true_Omega:
        h1sq = float(save_dir.split('/')[-1].split('-')[0][5:])
        h2sq = float(save_dir.split('/')[-1].split('-')[1][5:])  
        rg = float(save_dir.split('/')[-1].split('-')[2][3:])
        Omega = np.array([[h1sq, rg*np.sqrt(h1sq*h2sq)], [rg*np.sqrt(h1sq*h2sq), h2sq]]) / M 
        np.save(os.path.join(save_dir, f'Omega-{method}-{strOmega}-{r}.npy'), Omega)
    else:
        Nc, Nt = 100, 1000
        sumstats_pop1_names = ['cell1']
        sumstats_pop2_names = ['cell2']
        sumstats_pop1 = {'cell1': pd.DataFrame({'Z': zs[:,0], 'N': Nc})}
        sumstats_pop2 = {'cell2': pd.DataFrame({'Z': zs[:,1], 'N': Nc})}
        Omega, Sigma_LD, Omega_se = ldsc.run_ldscore_regressions_tram(sumstats_pop1_names,sumstats_pop1,sumstats_pop2_names,\
                sumstats_pop2,ld,ld,ld,intercept=np.eye(2))

        ## numerical operations
        Omega = Omega.squeeze()
        Omega_se = Omega_se.squeeze()
        for i in range(len(Omega)):
            Omega[i, i] = min(Omega[i, i], max_hsq/M) # max heritability: 0.8
        Omega = make_positive_definite(Omega, 0.95) # max correlation: 0.95
        oriOmega = Omega.copy() * M
        oriOmega_se = Omega_se.copy() * M
    
        ## hij test
        if trun_corr:
            n0 = int(K * (K - 1) / 2)
            n = 0
            for i in range(0, K):
                for j in range(i+1, K):
                    hij_z = Omega[i, j] / Omega_se[i, j]
                    hij_p = st.norm.sf(abs(hij_z)) * 2
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
        betas_tensor, ses_tensor = run_tram(betas[:,:2], ses[:,:2], ld, Omega=Omega, Sigma=Sigma)
        
    zs_tensor = betas_tensor / ses_tensor
    pval_tensor = st.norm.sf(abs(zs_tensor)) * 2
    cau_tensor = (pval_tensor < 0.05).astype('int32')
    
    ## save results
    # ours/tram cell1 estimators
    cell1_meta = pd.DataFrame({'SNP': np.arange(M),
                               'Beta': betas_tensor[:, 0],
                               'SE': ses_tensor[:, 0],
                               'Zscore': zs_tensor[:, 0],
                               'Pval': pval_tensor[:, 0],
                               'Sig': cau_tensor[:, 0]})
    # ours/tram cell2 estimators
    cell2_meta = pd.DataFrame({'SNP': np.arange(M),
                               'Beta': betas_tensor[:, 1],
                               'SE': ses_tensor[:, 1],
                               'Zscore': zs_tensor[:, 1],
                               'Pval': pval_tensor[:, 1],
                               'Sig': cau_tensor[:, 1]})
    filename1 = f'cell1-{method}-{strOmega}-{r}.csv'
    filename2 = f'cell2-{method}-{strOmega}-{r}.csv'
    if trun_corr:
        filename1 = '-'.join(filename1.split('-')[:-1] + ['truncorr', filename1.split('-')[-1]])
        filename2 = '-'.join(filename2.split('-')[:-1] + ['truncorr', filename2.split('-')[-1]])
    if zero_corr:
        filename1 = '-'.join(filename1.split('-')[:-1] + ['zerocorr', filename1.split('-')[-1]])
        filename2 = '-'.join(filename2.split('-')[:-1] + ['zerocorr', filename2.split('-')[-1]])
    cell1_meta.to_csv(os.path.join(save_dir, filename1), sep='\t', index=False)
    cell2_meta.to_csv(os.path.join(save_dir, filename2), sep='\t', index=False)

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--out-root', type=str, help='root directory for saving simulation data and results', required=True)
    parser.add_argument('--metric', type=str, choices=['type_1_error', 'power'], help='type_1_error or power experiment', required=True)
    parser.add_argument('--method', type=str, choices=['blue', 'tram'], help='use ours (blue) or tram', required=True)
    parser.add_argument('--win', type=int, help='window id', default=0, required=True)
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
    # h2sqs = [0.8] # debug
    pi1s = [0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99]
    pcaus = [0.01, 0.02, 0.05] # lower pcau
    nrep = 5
    start = 0
    reps = np.arange(start+1, start+nrep+1)
    metric = args.metric
    true_Omega = args.true_Omega
    method = args.method 
    parallel = args.parallel
    zero_corr = args.zero_corr
    trun_corr = args.trun_corr
    out_root = args.out_root
    win_id = args.win
    ieqtl = args.ieqtl
    omega_shrink = args.omega_shrink
    max_hsq = args.max_hsq
    
    # different hyperparameters
    strgeno = 'real_geno'
    if metric == 'type_1_error':
        h1sqs = [0.0]
        rgs = [0.0]
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
    win_ids = [win_id]
    for win_id in win_ids:
        G, ld = get_genotype_ldscore(win_id)

        ## generate data and save
        t0 = time.time()
        print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] ***** Start generating data (window {win_id}) *****")
        for h1sq in h1sqs:
            for h2sq in h2sqs:
                for rg in rgs:
                    for pi1 in pi1s:
                        for pcau in pcaus:
                            out_dir = os.path.join(out_root, f'h1sq={h1sq}-h2sq={h2sq}-rg={rg}-pi1={pi1}-win{win_id}-Nc{Nc}-Nt{Nt}-pcau{pcau}-seed{seed}')
                            os.makedirs(out_dir, exist_ok=True)
                            stream = tqdm(reps, desc=f'win={win_id}, h1sq={h1sq}, h2sq={h2sq}, rg={rg}, pi1={pi1}, pcau={pcau}', ascii=' >=')
                            for r in stream:
                                if os.path.exists(os.path.join(out_dir, f'tissue-sumstats-{r}.csv')):
                                    continue
                                else:
                                    generate_data(h1sq, h2sq, rg, pi1, Nc, Nt, pcau, r, out_dir, seed, G, ieqtl)
        print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] ***** Finish generating data (window {win_id}) *****")
        print(f"Generating data time (unparallel): {time.time() - t0:.2f}s")

        ## perform ours/tram
        if parallel:
            t0 = time.time()
            print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] ***** Start perform method: {method} (window {win_id}) *****")
            for h1sq in h1sqs:
                for h2sq in h2sqs:
                    for rg in rgs:
                        for pi1 in pi1s:
                            for pcau in pcaus:
                                out_dir = os.path.join(out_root, f'h1sq={h1sq}-h2sq={h2sq}-rg={rg}-pi1={pi1}-win{win_id}-Nc{Nc}-Nt{Nt}-pcau{pcau}-seed{seed}')
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
            print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] ***** Finish performing method: {method} (window {win_id}) *****")
            print(f"Performing method time (parallel): {time.time() - t0:.2f}s")
        else:
            t0 = time.time()
            print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] ***** Start perform method: {method} (window {win_id}) *****")
            for h1sq in h1sqs:
                for h2sq in h2sqs:
                    for rg in rgs:
                        for pi1 in pi1s:
                            for pcau in pcaus:
                                out_dir = os.path.join(out_root, f'h1sq={h1sq}-h2sq={h2sq}-rg={rg}-pi1={pi1}-win{win_id}-Nc{Nc}-Nt{Nt}-pcau{pcau}-seed{seed}')
                                stream = tqdm(reps, desc=f'win={win_id}, h1sq={h1sq}, h2sq={h2sq}, rg={rg}, pi1={pi1}, pcau={pcau}', ascii=' >=')
                                for r in stream:
                                    simu(r, out_dir, true_Omega, omega_shrink, max_hsq, method, zero_corr, trun_corr)
            print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] ***** Finish performing method: {method} (window {win_id}) *****")
            print(f"Performing method time (unparallel): {time.time() - t0:.2f}s")
