import gzip
import os
import pandas as pd
import logging
import numpy as np
import sys
import argparse
import re
import glob
import logging
from collections import defaultdict
import matplotlib as mpl
import scanpy as sc

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s') 

def run_NMF(sample_name,matrix_10X,threads,K_range,K_selection,density_threshold,iteration,run_K):
    if not run_K:
        logging.info("reading matrix")
        adata = sc.read_10x_mtx(matrix_10X,var_names='gene_symbols',cache=False)
        adata.var_names_make_unique()
        outdir = "NMF_out/"
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        count_adat_fn = outdir + sample_name + '.h5ad'
        logging.info("writing h5 file")
        sc.write(count_adat_fn, adata)

        numiter=iteration # Number of NMF replicates. Set this to a larger value ~200 for real data. We set this to a relatively low value here for illustration at a faster speed
        numhvgenes=2000 ## Number of over-dispersed genes to use for running the actual factorizations
        ## Results will be saved to [output_directory]/[run_name] which in this example is example_PBMC/cNMF/pbmc_cNMF
        seed = 0 ## Specify a seed pseudorandom number generation for reproducibility

        numworkers= threads
        prepare_cmd = """python /SGRNJ/Database/script/soft/cNMF/cnmf.py \
        prepare --output-dir %s --name %s -c %s -k %s --n-iter %d \
        --total-workers %d --seed %d --numgenes %d --beta-loss frobenius""" % (outdir, 
            sample_name, count_adat_fn, K_range, numiter, numworkers, seed, numhvgenes)
        logging.info('Prepare command assuming parallelization with %d cores:\n%s' % (numworkers, prepare_cmd))
        os.system(prepare_cmd)

        ## Using GNU parallel
        worker_index = ' '.join([str(x) for x in range(numworkers)])
        factorize_cmd = """parallel python /SGRNJ/Database/script/soft/cNMF/cnmf.py \
            factorize --output-dir %s --name %s --worker-index {} ::: %s""" % (outdir, sample_name, worker_index)
        logging.info('Factorize command to simultaneously run factorization over %d cores using GNU parallel:\n%s' % (numworkers, factorize_cmd))
        os.system(factorize_cmd)

        # combine
        combine_cmd = 'python /SGRNJ/Database/script/soft/cNMF/cnmf.py \
            combine --output-dir %s --name %s' % (outdir, sample_name)
        logging.info(combine_cmd)
        os.system(combine_cmd)

        worker_index = ' '.join([str(x) for x in range(numworkers)])
        kselect_plot_cmd = 'python /SGRNJ/Database/script/soft/cNMF/cnmf.py \
            k_selection_plot --output-dir %s --name %s' % (outdir, sample_name)
        logging.info('K selection plot command: %s' % kselect_plot_cmd)
        os.system(kselect_plot_cmd)

    # run_K
    consensus_cmd = 'python /SGRNJ/Database/script/soft/cNMF/cnmf.py \
        consensus --output-dir %s --name %s --local-density-threshold %.2f \
        --components %d --show-clustering' % (outdir, sample_name, density_threshold, K_selection)
    logging.info('Consensus command for K=%d:\n%s' % (K_selection, consensus_cmd))
    os.system(consensus_cmd)

    ## Load the Z-scored GEPs which reflect how enriched a gene is in each GEP relative to all of the others
    density_threshold_str = ('%.2f' % density_threshold).replace('.', '_')
    gene_file = '{outdir}/{sample_name}/{sample_name}.gene_spectra_score.k_{K_selection}.dt_{density_threshold_str}.txt'.format(outdir=outdir,
        sample_name=sample_name,K_selection=K_selection,density_threshold_str=density_threshold_str)
    gene_scores = pd.read_csv(gene_file,sep='\t', index_col=0).T

    ## Obtain the top 100 genes for each GEP in sorted order and combine them into a single dataframe
    top_genes = []
    ngenes = 100
    for gep in gene_scores.columns:
        top_genes.append(list(gene_scores.sort_values(by=gep, ascending=False).index[:ngenes]))
        
    top_genes = pd.DataFrame(top_genes, index=gene_scores.columns).T
    top_genes_file = '{outdir}/{sample_name}_top100_genes.tsv'.format(outdir=outdir,sample_name=sample_name)
    top_genes.to_csv(top_genes_file,sep="\t")

    usage_file = '{outdir}/{sample_name}/{sample_name}.usages.k_{K_selection}.dt_{density_threshold_str}.consensus.txt'.format(outdir=outdir,
        sample_name=sample_name,K_selection=K_selection,density_threshold_str=density_threshold_str)
    logging.info("usage_file:"+usage_file)
    usage = pd.read_csv(usage_file, sep='\t', index_col=0)
    usage.columns = ['Usage_%s' % i for i in usage.columns]
    usage_norm = usage.div(usage.sum(axis=1), axis=0)
    usage_norm_file = '{outdir}/{sample_name}_usage_norm.tsv'.format(outdir=outdir,sample_name=sample_name)
    usage_norm.to_csv(usage_norm_file,sep="\t")

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NMF')
    parser.add_argument("--sample_name", help="sample name")
    parser.add_argument("--matrix_10X", help="10X matrix dir")
    parser.add_argument("--threads", default=5)
    parser.add_argument("--K_range", help="K range",default="7 8")
    parser.add_argument("--K_selection", help="K selection",default=7)
    parser.add_argument("-d","--density_threshold", help="--density threshold",default = 0.2)
    parser.add_argument("--iteration", help="iteration",default=100)
    parser.add_argument("--run_K", help="only run K selection",action="store_true")
    args = parser.parse_args()

    sample_name = args.sample_name
    matrix_10X = args.matrix_10X
    threads = int(args.threads)
    K_range = args.K_range
    K_selection = int(args.K_selection)
    density_threshold = float(args.density_threshold)
    iteration = int(args.iteration)
    run_K = args.run_K

    run_NMF(sample_name,matrix_10X,threads,K_range,K_selection,density_threshold, iteration, run_K)