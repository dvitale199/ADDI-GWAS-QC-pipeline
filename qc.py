import subprocess
import os
import sys
import shutil
import pandas as pd
import numpy as np


def shell_do(command, log=False, return_log=False):
    print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res=subprocess.run(command.split(), stdout=subprocess.PIPE)

    if log:
        print(res.stdout.decode('utf-8'))
    if return_log:
        return(res.stdout.decode('utf-8'))
        

def callrate_prune(geno_path, out_path, mind=0.02):
    
    # what step are we running?
    step = "callrate_prune"
    print()
    print(f"RUNNING: {step}")
    print()
    
    outliers_out = f'{out_path}.outliers'
    
    plink_cmd1 = f"plink2 --pfile {geno_path} --mind {mind} --make-pgen --out {out_path}"

    shell_do(plink_cmd1)
    
    if os.path.isfile(f'{out_path}.mindrem.id'):
        shutil.move(f'{out_path}.mindrem.id', outliers_out)

        outlier_count = sum(1 for line in open(f'{outliers_out}'))
        
    else:
        outlier_count = 0
        
    process_complete = True
    
    outfiles_dict = {
        'pruned_samples': f'{outliers_out}',
        'plink_out': f'{out_path}'
    }

    metrics_dict = {
        'outlier_count': outlier_count
    }

    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict




def het_prune(geno_path, out_path):
    
    # what step are we running?
    step = "het_prune"
    print()
    print(f"RUNNING: {step}")
    print()

    het_tmp = f"{out_path}_tmp"
    het_tmp2 = f"{out_path}_tmp2"
    het_tmp3 = f"{out_path}_tmp3"
    outliers_out = f"{out_path}.outliers"
    
    plink_cmd1 = f"plink2 --pfile {geno_path} --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out {het_tmp}"
    plink_cmd2 = f"plink2 --pfile {geno_path} --extract {het_tmp}.prune.in --make-pgen --out {het_tmp2}"
    plink_cmd3 = f"plink2 --pfile {het_tmp2} --het --out {het_tmp3}"

    cmds1 = [plink_cmd1, plink_cmd2, plink_cmd3]

    for cmd in cmds1:
        shell_do(cmd)
    
    # check if het_tmp3 is created. if not, skip this.
    # NOTE: there may be legitimate reasons for this, for example, if only one sample in genotype file (sometimes happens in ancestry split method)
    hetpath = f'{het_tmp3}.het'
    if os.path.isfile(hetpath):
        het = pd.read_csv(hetpath, sep='\s+')
        het_outliers = het[((het.F <= -0.25) | (het.F >= 0.25))]
        outlier_count = het_outliers.shape[0]
        het_outliers.to_csv(f'{outliers_out}', sep='\t', index=False)

        plink_cmd4 = f"plink2 --pfile {geno_path} --remove {outliers_out} --make-pgen --out {out_path}"

        shell_do(plink_cmd4)

        outfiles_dict = {
            'pruned_samples': outliers_out,
            'plink_out': out_path
        }

        metrics_dict = {
            'outlier_count': outlier_count
        }
        
        process_complete = True
    
    else:
        print(f'Heterozygosity pruning failed!')
        print(f'Check {het_tmp}.log, {het_tmp2}.log, or {het_tmp3}.log for more information')
        
        outfiles_dict = {
            'pruned_samples': 'Heterozygosity Pruning Failed!',
            'plink_out': [het_tmp, het_tmp2, het_tmp3]
        }

        metrics_dict = {
            'outlier_count': 0
        }
    
        process_complete = False
    
    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict


def merge_genos(geno_path1, geno_path2, out_name):
    # attempt 1 at merging genos
    bash1 = f"plink --bfile {geno_path1} --allow-no-sex --bmerge {geno_path2} --out {out_name} --make-bed"
    shell_do(bash1)

    # if {outname}-merge.missnp file created, snps need to be flipped and merge tried again
    if os.path.isfile(f'{out_name}-merge.missnp'):
        bash2 = f"plink --bfile {geno_path1} --allow-no-sex --flip {out_name}-merge.missnp --make-bed --out {geno_path1}_flip"
        bash3 = f"plink --bfile {geno_path1}_flip --allow-no-sex --bmerge {geno_path2} --out {out_name}_flip --make-bed"

        cmds1 = [bash2, bash3]

        for cmd in cmds1:
            shell_do(cmd)

        #if another -merge.missnp file is created, these are likely triallelic positions and must be excluded and then try merge again
        if os.path.isfile(f'{geno_path1}_flip-merge.missnp'):
            bash4 = f"plink --bfile {geno_path1}_flip --allow-no-sex --exclude {out_name}_flip-merge.missnp --out {geno_path1}_flip_pruned --make-bed"
            bash5 = f"plink --bfile {geno_path1}_flip_pruned --allow-no-sex --bmerge {geno_path2} --out  {out_name} --make-bed"

            cmds2 = [bash4, bash5]

            for cmd in cmds2:
                shell_do(cmd)

        # if second attempt at merge is successful, there are no triallelic snps and we can go ahead and move the _flip files to our _merged_ref_panel filenames 
        # for further processing
        else:
            suffix_list = ['bed','bim','fam']
            for suffix in suffix_list:
                shutil.copy(f'{out_name}_flip.{suffix}',f'{out_name}.{suffix}')

    # if first merge works, there are no alleles in need of flipping and no triallelic positions and we can proceed!
    else:
        pass

def get_outlier_ranges(pc_df, target_label):

    target_pc1 = pc_df.loc[pc_df.label==target_label,'PC1']
    target_pc2 = pc_df.loc[pc_df.label==target_label,'PC2']

    targetlowc1 = np.mean(target_pc1) - (6*np.std(target_pc1))
    targethighc1 = np.mean(target_pc1) + (6*np.std(target_pc1))
    targetlowc2 = np.mean(target_pc2) - (6*np.std(target_pc2))
    targethighc2 = np.mean(target_pc2) + (6*np.std(target_pc2))

    out_dict = {
        'lowc1': targetlowc1,
        'highc1': targethighc1,
        'lowc2': targetlowc2,
        'highc2': targethighc2
    }

    return out_dict

def ancestry_prune(geno_path, ref_path, ref_labels, out_path, target_label):
    merge_genos(geno_path, ref_path, out_path)

    plink_cmd1 = f"plink --bfile {out_path} --geno 0.01 --pca 4 --make-bed --out {out_path}.pca"
    shell_do(plink_cmd1)

    pcs = pd.read_csv(f'{out_path}.pca.eigenvec', sep='\s+', header=None, names=['FID','IID','PC1','PC2','PC3','PC4'], dtype={'FID':str,'IID':str})
    fam = pd.read_csv(f'{geno_path}.fam', sep='\s+', header=None, usecols=[0,1], names=['FID','IID'], dtype={'FID':str,'IID':str})
    fam.loc[:,'label'] = 'new'
    ref_lab_df = pd.read_csv(f'{ref_labels}', sep='\t', header=None, names=['FID','IID','label'])
    label_df = fam.append(ref_lab_df)
    ref_pcs_merge = pcs.merge(label_df, how='left', on=['FID','IID'])

    outlier_ranges = get_outlier_ranges(ref_pcs_merge, target_label)

    fam_pcs = fam.merge(pcs, how='left', on=['FID','IID'])
    fam_keep = fam_pcs.loc[(fam_pcs.PC1 >= outlier_ranges['lowc1']) & (fam_pcs.PC1 <= outlier_ranges['highc1']) & (fam_pcs.PC2 >= outlier_ranges['lowc2']) & (fam_pcs.PC2 <= outlier_ranges['highc2'])]
    fam_out_merge = fam_pcs.merge(fam_keep, how='left', on=['FID','IID'], indicator=True)
    fam_outliers = fam_out_merge.loc[fam_out_merge['_merge']=='left_only']
    fam_keep[['FID','IID']].to_csv(f'{out_path}.keep', sep='\t', header=True, index=False)
    fam_outliers[['FID','IID']].to_csv(f'{out_path}.outliers', sep='\t', header=True, index=False)
    outlier_count = fam_outliers.shape[0]

    plink_cmd2 = f'plink --bfile {geno_path} --keep {out_path}.keep --make-bed --out {out_path}'
    shell_do(plink_cmd2)

    process_complete = True
    step = 'ancestry_prune'

    outfiles_dict = {
        'pruned_samples': f'{out_path}.outliers',
        'plink_out': f'{out_path}',
        'keep_samples': f'{out_path}.keep'
    }

    metrics_dict = {
        'outlier_count': outlier_count
    }

    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict


def king_prune(geno_path, out_path, king_cutoff=0.0884):
    plink_cmd1 = f'plink2 --pfile {geno_path} --hwe 0.0001 --mac 2 --make-pgen --out {out_path}_tmp'
    
    plink_cmd2 = f'plink2 --pfile {out_path}_tmp --king-cutoff {king_cutoff} --out {out_path}' 
    
    plink_cmd3 = f'plink2 --pfile {out_path}_tmp --remove {out_path}.king.cutoff.out.id --make-pgen --out {out_path}'

    for cmd in [plink_cmd1, plink_cmd2, plink_cmd3]:
        shell_do(cmd)

    shutil.copy(f'{out_path}.king.cutoff.out.id',f'{out_path}.outliers')
    
    outlier_count = sum(1 for line in open(f'{out_path}.outliers'))

    process_complete = True
    step = 'king_prune'

    outfiles_dict = {
        'pruned_samples': f'{out_path}.outliers',
        'plink_out': f'{out_path}'
    }

    metrics_dict = {
        'outlier_count': outlier_count
    }

    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict


def report_outliers(steps, outpath):

    outliers_df = pd.DataFrame()
    
    for step in steps:
        outlier_file = step['output']['pruned_samples']

        if step['pass']:
            print(f"{step['step']} Passed!")
            if os.path.isfile(outlier_file):
                step_df = pd.read_csv(outlier_file, sep='\s+')
                step_df.rename(columns={'#FID':'FID','#IID':'IID'}, inplace=True)
                
                if step_df.shape[0] > 0:
                    if 'FID' not in step_df:
                        step_df.loc[:,'FID'] = '0'
                
                    outliers_df = outliers_df.append(step_df[['FID','IID']])
        else:
            print(f"{step['step']} failed!")
    
    outliers_df.to_csv(outpath, sep='\t', header=False, index=False)
    out_dict = {'outliers_df': outliers_df, 'outliers_path': outpath}

    return out_dict