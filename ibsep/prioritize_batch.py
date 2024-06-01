import os, sys, time, glob, argparse, traceback
from datetime import datetime
from tqdm import tqdm
import numpy as np
import pandas as pd
import scipy.stats as st
from ibsep import core, ldsc, utils

def run_experiment(data_dir, out_dir, gene_id, cell_type_df, pval_thres, trun_corr=True, zero_corr=False, reg_int_ident=True):

    os.makedirs(out_dir, exist_ok=True)
    all_cell_types = cell_type_df.cell_type.values.tolist()
    sumstat_files = [x.split('/')[-1] for x in glob.glob(f'{data_dir}/{gene_id}*')]
    sumstat_files = [x for x in sumstat_files if x.replace(f'{gene_id}_', '').replace('.csv', '') in all_cell_types]
    cell_types = sorted([x.replace(f'{gene_id}_', '').replace('.csv', '') for x in sumstat_files])
    indices = [all_cell_types.index(x) for x in cell_types]
    K = len(cell_types)
    props = cell_type_df.iloc[indices, :]['prop'].values
        
    ## load data (TODO: merge sumstats with tissue in case different SNP order)
    common_snps = []
    sumstats = {}
    for cell_type in cell_types:
        sumstat = pd.read_csv(os.path.join(data_dir, f'{gene_id}_{cell_type}.csv'), sep='\t')
        sumstat.drop_duplicates(subset=['SNP'], inplace=True)
        sumstat['Z'] = sumstat.BETA / sumstat.SE
        sumstat['N'] = 1 / sumstat.SE ** 2
        sumstats[cell_type] = sumstat
        common_snps.append(sumstat.SNP.unique().tolist())
    tissue = pd.read_csv(os.path.join(data_dir, f'{gene_id}_gtex_withld.csv'), sep='\t')
    tissue.drop_duplicates(subset=['SNP'], inplace=True)
    tissue['Z'] = tissue.BETA / sumstat.SE
    tissue['N'] = 1 / tissue.SE ** 2
    common_snps.append(tissue.SNP.unique().tolist())
    # keep common snps
    common_snps = list(set.intersection(*map(set, common_snps)))
    for cell_type in cell_types:
        sumstats[cell_type] = sumstats[cell_type].loc[sumstats[cell_type].SNP.isin(common_snps)].sort_values(by='BP').reset_index(drop=True)
    tissue = tissue.loc[tissue.SNP.isin(common_snps)].sort_values(by='BP').reset_index(drop=True)
    ld = tissue.LD.values[:,None]
    M = len(ld)
    
    ## inputs
    zs_ = {k: v.Z.values for k, v in sumstats.items()} | {'gtex': tissue.Z.values}
    betas_ = {k: v.BETA.values for k, v in sumstats.items()} | {'gtex': tissue.BETA.values}
    ses_ = {k: v.SE.values for k, v in sumstats.items()} | {'gtex': tissue.SE.values}
    zs = np.stack(list(zs_.values()), axis=1)
    betas = np.stack(list(betas_.values()), axis=1)
    ses = np.stack(list(ses_.values()), axis=1)

    ## Omega and Sigma
    intercept = np.eye(K+1) if reg_int_ident else None
    Omega, Sigma, Omega_se = ldsc.run_ldscore_regressions_multi(sumstats,ld,tissue,intercept=intercept)
    if reg_int_ident:
        Sigma = np.eye(K+1)

    ## numerical operations
    Omega = Omega.squeeze()[:K,:K]
    Omega_se = Omega_se.squeeze()[:K,:K]
    oriOmega = Omega.copy() * M
    oriOmega_se = Omega_se.copy() * M
    for i in range(len(Omega)):
        Omega[i, i] = min(Omega[i, i], 0.8/M) # max heritability: 0.8
    Omega = utils.restrict_corr(Omega, 0.95) # max correlation: 0.95

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
    filename = f'{gene_id}_Omega_Sigma.npz'
    if not np.all(np.linalg.eigvals(Omega) > 0):
        Omega, Omega_se, n_shrink = utils.make_omega_pd(Omega, Omega_se)
        shrinkOmega = Omega.copy() * M
        shrinkOmega_se = Omega_se.copy() * M
        np.savez(os.path.join(out_dir, filename), oriOmega=oriOmega, oriOmega_se=oriOmega_se, 
                 shrinkOmega=shrinkOmega, shrinkOmega_se=shrinkOmega_se, n_shrink=n_shrink,
                 Sigma=Sigma, cell_types=cell_types)
    else:
        np.savez(os.path.join(out_dir, filename), oriOmega=oriOmega, oriOmega_se=oriOmega_se, 
                 Sigma=Sigma, cell_types=cell_types)
    
    ## ibsep
    betas_tensor, ses_tensor = core.run_ibsep(betas, ses, props, ld, Omega=Omega[:K,:K], Sigma=Sigma)
        
    zs_tensor = betas_tensor / ses_tensor
    pval_tensor = st.norm.sf(abs(zs_tensor)) * 2
    
    ## save results
    for i, cell_type in enumerate(cell_types):
        sumstat = sumstats[cell_type]
        sumstat['BETA_IBSEP'] = betas_tensor[:, i]
        sumstat['SE_IBSEP'] = ses_tensor[:, i]
        sumstat['PVAL_IBSEP'] = pval_tensor[:, i]
        sumstat.drop(columns=['Z', 'N'], inplace=True)
        filename = f'{gene_id}_{cell_type}_IBSEP'
        if trun_corr:
            filename = f'{gene_id}_{cell_type}_IBSEP_truncorr_pval{pval_thres}'
        if zero_corr:
            filename = f'{gene_id}_{cell_type}_IBSEP_zerocorr'
        if not reg_int_ident:
            filename = f'{filename}_scs'
        sumstat.to_csv(os.path.join(out_dir, f'{filename}.csv'), sep='\t', index=False)

    return K, cell_types, props
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', type=str, help='directory saving sumstats and ldscore', required=True)
    parser.add_argument('--out-dir', type=str, help='directory for saving results', required=True)
    parser.add_argument('--avg-props', type=str, help='average cell type proportion dataframe', required=True)
    parser.add_argument('--pval-thres', type=float, help='threshold of hij test', default=1e-10)
    parser.add_argument('--zero-corr', action='store_true', help='force hij to be zero')
    parser.add_argument('--trun-corr', action='store_true', help='truncate insignificant hij to be zero')
    parser.add_argument('--not-reg-int-ident', action='store_true', help="do not force intercept matrix to be identity")
    args = parser.parse_args()

    ## load cell type proportion
    cell_type_df = pd.read_csv(args.avg_props, sep='\t')
    cell_type_df.cell_type = cell_type_df.cell_type.apply(lambda x: x.replace(' ', '_'))
    cell_type_df = cell_type_df.sort_values(by='cell_type').reset_index(drop=True)
    cell_types = cell_type_df.cell_type.values.tolist()
    props = cell_type_df.prop.values
    print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] {len(cell_types)} cell types in total:")
    print(cell_type_df)

    ## gene IDs
    files = glob.glob(f'{args.data_dir}/ENSG*_gtex_withld.csv')
    gene_ids = [x.split('/')[-1].split('_')[0] for x in files]
    reg_int_ident = True
    if args.not_reg_int_ident:
        reg_int_ident = False
    print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] {len(gene_ids)} genes.")
    
    ## perform IBSEP
    t0 = time.time()
    stream = tqdm(gene_ids, ascii=' >=')
    for gene_id in stream:
        try:
            K, current_cell_types, current_props = run_experiment(args.data_dir, args.out_dir, gene_id, cell_type_df, args.pval_thres, args.trun_corr, args.zero_corr, reg_int_ident)
        except:
            print(f'!!! Error {gene_id} !!!')
            traceback.print_exc()
        abr_cell_props = [x[:4] for x, y in zip(current_cell_types, current_props)]
        stream.set_description(f'{gene_id}, {K} cell types: {abr_cell_props}')
    print(f"[{datetime.today().strftime('%Y-%m-%d %H:%M:%S')}] IBSEP running time: {time.time() - t0:.1f}s")
    


