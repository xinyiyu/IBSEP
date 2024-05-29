import os, sys, time, glob, argparse, traceback

threads = '4'
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


def run_experiment(data_dir, out_dir, gene_id, cell_type_df, pval_thres, tau, omega_shrink, max_hsq, method=['blue', 'tram'], trun_corr=True, zero_corr=False, return_w1=False, reg_int_ident=False):

    os.makedirs(out_dir, exist_ok=True)
    all_cell_types = cell_type_df.cell_type.values.tolist()
    sumstat_files = [x.split('/')[-1] for x in glob.glob(f'{data_dir}/{gene_id}*')]
    sumstat_files = [x for x in sumstat_files if x.replace(f'{gene_id}_', '').replace('.csv', '') in all_cell_types]
    cell_types = sorted([x.replace(f'{gene_id}_', '').replace('.csv', '') for x in sumstat_files])
    indices = [all_cell_types.index(x) for x in cell_types]
    K = len(cell_types)
    props = cell_type_df.iloc[indices, :]['prop'].values
    # props = cell_type_df.loc[cell_type_df.cell_type.isin(cell_types), 'prop'].values
    props = props * tau # tau: tissue weight
    # print(f'{K} cell types for {gene_id}:\n', pd.DataFrame({'cell_type': cell_types, 'prop': props}))
    # print(props)
        
    ## load data (TODO: merge sumstats with tissue in case different SNP order)
    sumstats = {}
    for cell_type in cell_types:
        sumstat = pd.read_csv(os.path.join(data_dir, f'{gene_id}_{cell_type}.csv'), sep='\t')
        sumstat['Z'] = sumstat.BETA / sumstat.SE
        sumstat['N'] = 1 / sumstat.SE ** 2
        sumstats[cell_type] = sumstat.sort_values(by='BP').reset_index(drop=True)
    tissue = pd.read_csv(os.path.join(data_dir, f'{gene_id}_gtex_withld.csv'), sep='\t')
    tissue = tissue.sort_values(by='BP').reset_index(drop=True)
    tissue['Z'] = tissue.BETA / sumstat.SE
    tissue['N'] = 1 / tissue.SE ** 2
    ld = tissue.LD.values[:,None]
    M = len(ld)
    # print(f'Loaded tissue sumstats, {len(sumstats)} cell sumstats and ldscore: {M} SNPs.')
    
    ## inputs
    zs_ = {k: v.Z.values for k, v in sumstats.items()} | {'gtex': tissue.Z.values}
    betas_ = {k: v.BETA.values for k, v in sumstats.items()} | {'gtex': tissue.BETA.values}
    ses_ = {k: v.SE.values for k, v in sumstats.items()} | {'gtex': tissue.SE.values}
    zs = np.stack(list(zs_.values()), axis=1)
    betas = np.stack(list(betas_.values()), axis=1)
    ses = np.stack(list(ses_.values()), axis=1)
    # print(f'Data shapes: zscores: {zs.shape}, betas: {betas.shape}, ses: {ses.shape}')

    ## Omega and Sigma
    if reg_int_ident:
        # intercept = np.eye(K+1) if method == 'blue' else np.eye(K)
        intercept = np.eye(K+1)
    else:
        intercept = None
    Omega, Sigma, Omega_se = ldsc.run_ldscore_regressions_tram_multi(sumstats,ld,tissue,intercept=intercept)
    if reg_int_ident:
        Sigma = np.eye(K+1) if method == 'blue' else np.eye(K)

    ## numerical operations
    Omega = Omega.squeeze()[:K,:K]
    Omega_se = Omega_se.squeeze()[:K,:K]
    oriOmega = Omega.copy() * M
    oriOmega_se = Omega_se.copy() * M
    for i in range(len(Omega)):
        Omega[i, i] = min(Omega[i, i], max_hsq/M) # max heritability: 0.8
    Omega = make_positive_definite(Omega, 0.95) # max correlation: 0.95

    if not reg_int_ident:
        filename = f'{gene_id}_Sigma_{method}_scs.npz'
        np.savez(os.path.join(out_dir, filename), Sigma=Sigma, cell_types=cell_types)

    ## hij test
    if trun_corr:
        n0 = int(K * (K - 1) / 2)
        n = 0
        for i in range(0, K):
            for j in range(i+1, K):
                hij_z = Omega[i, j] / Omega_se[i, j]
                hij_p = st.norm.sf(abs(hij_z)) * 2
                if hij_p >= pval_thres:
                    Omega[i, j] = 1e-32
                    Omega[j, i] = 1e-32
                    n += 1
    if zero_corr:
        for i in range(0, K):
            for j in range(i+1, K):
                Omega[i, j] = 1e-32
                Omega[j, i] = 1e-32

    ## make pd
    filename = f'{gene_id}_Omega_{method}.npz'
    if not np.all(np.linalg.eigvals(Omega) > 0):
        if omega_shrink < 1:
            Omega, Omega_se, n_shrink = make_pd_1(Omega, Omega_se, omega_shrink)
            shrinkOmega = Omega.copy() * M
            shrinkOmega_se = Omega_se.copy() * M
            np.savez(os.path.join(out_dir, filename), oriOmega=oriOmega, oriOmega_se=oriOmega_se, 
                     shrinkOmega=shrinkOmega, shrinkOmega_se=shrinkOmega_se, cell_types=cell_types, 
                     omega_shrink=omega_shrink, n_shrink=n_shrink)
        else:
            Omega = make_pd_2(Omega)
            pdOmega = Omega.copy() * M
            np.savez(os.path.join(out_dir, filename), oriOmega=oriOmega, oriOmega_se=oriOmega_se, 
                     pdOmega=pdOmega, cell_types=cell_types)
    else:
        np.savez(os.path.join(out_dir, filename), oriOmega=oriOmega, oriOmega_se=oriOmega_se, 
                 cell_types=cell_types)
    
    ## ours or tram
    if method == 'blue':
        if return_w1:
            betas_tensor, ses_tensor, w1_tensor = run_ours(betas, ses, props, ld, Omega=Omega[:K,:K], Sigma=Sigma, return_w1=True)
            np.save(os.path.join(out_dir, f'w1-{method}-{r}.npy'), w1_tensor)
        else:
            # print(betas.shape, ses.shape, props.shape, ld.shape, Omega.shape, Sigma.shape)
            betas_tensor, ses_tensor = run_ours(betas, ses, props, ld, Omega=Omega[:K,:K], Sigma=Sigma, return_w1=False)
    elif method == 'tram':
        betas_tensor, ses_tensor = run_tram(betas[:,:K], ses[:,:K], ld, Omega=Omega[:K,:K], Sigma=Sigma)
    # print(f'Calculated new [betas, ses] by {method}.')
        
    zs_tensor = betas_tensor / ses_tensor
    pval_tensor = st.norm.sf(abs(zs_tensor)) * 2
    
    ## save results
    for i, cell_type in enumerate(cell_types):
        sumstat = sumstats[cell_type]
        sumstat['BETA_BLUE'] = betas_tensor[:, i]
        sumstat['SE_BLUE'] = ses_tensor[:, i]
        sumstat['PVAL_BLUE'] = pval_tensor[:, i]
        filename = f'{gene_id}_{cell_type}_{method}.csv'
        if trun_corr:
            filename = f'{gene_id}_{cell_type}_{method}_truncorr.csv'
            filename = '_'.join(filename.split('_')[:-1] + [f'pval{pval_thres}', filename.split('_')[-1]])
        if zero_corr:
            filename = f'{gene_id}_{cell_type}_{method}_zerocorr.csv'
        if not reg_int_ident:
            filename = '_'.join(filename.split('_')[:-1] + ['scs', filename.split('_')[-1]])
        filename = '_'.join(filename.split('_')[:-1] + [f'tau{tau}', filename.split('_')[-1]])
        sumstat.to_csv(os.path.join(out_dir, filename), sep='\t', index=False)
    # print(f'Saved result sumstats.')

    return K, cell_types, props
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', type=str, help='directory saving sumstats and ldscore', required=True)
    parser.add_argument('--out-dir', type=str, help='directory for saving results', required=True)
    parser.add_argument('--method', type=str, choices=['blue', 'tram'], help='use ours (blue) or tram', required=True)
    parser.add_argument('--cell-type-props', type=str, help='cell type proportion file, a dataframe', required=True)
    parser.add_argument('--gene-id', type=str, help='gene to be processed')
    parser.add_argument('--pval-thres', type=float, help='threshold of hij test', default=0.001)
    parser.add_argument('--omega-shrink', type=float, help='omega shrinkage factor', default=0.9)
    parser.add_argument('--max-hsq', type=float, help='omega shrinkage factor', default=0.8)
    parser.add_argument('--tau', type=float, help='tissue summary weight', default=1)
    parser.add_argument('--parallel', action='store_true', help='parallel')
    parser.add_argument('--zero-corr', action='store_true', help='force hij to be zero')
    parser.add_argument('--trun-corr', action='store_true', help='truncate insignificant hij to be zero (hij_insig=0)')
    parser.add_argument('--return-w1', action='store_true', help="return w1 of our method")
    parser.add_argument('--reg-int-ident', action='store_true', help="identity intercept matrix (recommended)")
    args = parser.parse_args()

    print(f'pval-thres: {args.pval_thres}, {type(args.pval_thres)}')

    ## load cell type proportion
    chrom = args.data_dir.split('/')[-1].replace('chr', '')
    cell_type_df = pd.read_csv(args.cell_type_props, sep='\t')
    cell_type_df.cell_type = cell_type_df.cell_type.apply(lambda x: x.replace(' ', '_'))
    cell_type_df = cell_type_df.sort_values(by='cell_type').reset_index(drop=True)
    cell_types = cell_type_df.cell_type.values.tolist()
    props = cell_type_df.prop.values
    print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] {len(cell_types)} cell types in total:")
    print(cell_type_df)
    if args.gene_id is not None:
        gene_id = args.gene_id
        print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] CHR {chrom} Gene id: {args.gene_id}")
        ## perform ours/tram
        t0 = time.time()
        print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] Start {args.method} ({gene_id})")
        run_experiment(args.data_dir, args.out_dir, args.gene_id, cell_type_df, args.pval_thres, args.tau, args.omega_shrink, args.max_hsq, args.method, args.trun_corr, args.zero_corr, args.return_w1, args.reg_int_ident)
        print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] Finish {args.method} ({gene_id})")
        print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] Total time: {time.time() - t0:.1f}s")
    else:
        files = glob.glob(f'{args.data_dir}/*gtex.csv')
        gene_ids = [x.split('/')[-1].split('_')[0] for x in files]
        print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] CHR {chrom} {len(gene_ids)} genes.")
        ## perform ours/tram
        t0 = time.time()
        stream = tqdm(gene_ids, ascii=' >=')
        for gene_id in stream:
            try:
                K, current_cell_types, current_props = run_experiment(args.data_dir, args.out_dir, gene_id, cell_type_df, args.pval_thres, args.tau, args.omega_shrink, args.max_hsq, args.method, args.trun_corr, args.zero_corr, args.return_w1, args.reg_int_ident)
            except:
                print(f'!!! Error {gene_id} !!!')
                traceback.print_exc()
            abr_cell_props = [(x[:3], round(y, 3)) for x, y in zip(current_cell_types, current_props)]
            stream.set_description(f'{gene_id}, {K} cell types: {abr_cell_props}')
        print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] CHR {chrom} Total time: {time.time() - t0:.1f}s")

