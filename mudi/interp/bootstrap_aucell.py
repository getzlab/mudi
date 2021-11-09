import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
from pyscenic.aucell import aucell
from .aucell import create_gene_signatures
from .aucell import assign_bootstrap

def main():
    parser = argparse.ArgumentParser(description='AUcell Bootstrapping.')
    parser.add_argument(
        '-i', '--in_file',
        required=True,
        help='<Required> Path to input expression matrix.',
        type=str
    )
    parser.add_argument(
        '-d', '--de_genes',
        required=True,
        help='<Required> Differential expression results.',
        type=str
    )
    parser.add_argument('-o', '--out_file',
        help='<Required> Output .h5 file to save results.',
        required=True,
        type=str
    )
    parser.add_argument('-n', '--niter',
        help='Number of iterations.',
        required=False,
        default=100,
        type=int
    )
    parser.add_argument('-s', '--subset_n',
        help='Number of genes to subset.',
        required=False,
        default=150,
        type=int
    )
    parser.add_argument('-w', '--n_workers',
        help='Number of workers.',
        required=False,
        default=8,
        type=int
    )
    parser.add_argument('-k', '--weight',
        help='Enrichment weight. Default is "t" statistic form differential expression.',
        required=False,
        default="t",
        type=str
    )
    parser.add_argument('-r', '--random_seed',
        help='Random seed for bootstrapping.',
        required=False,
        default=None
    )
    args = parser.parse_args()

    # Set random seed
    if args.random_seed is None:
        np.random.seed()
    else:
        np.random.seed(int(args.random_seed))

    # Load
    exp_mtx = pd.read_parquet(args.in_file)

    print("   * {} cells loaded".format(exp_mtx.shape[0]))
    print("   * {} genes detected".format(exp_mtx.shape[1]))

    # Load DE Genes
    de_df = pd.read_csv(args.de_genes, sep='\t').set_index("gene_name")

    store = pd.HDFStore(args.out_file,'a')

    for n in tqdm(range(args.niter)):
        gene_sigs = create_gene_signatures(de_df, n=args.subset_n, weight_idx=args.weight)
        enrich_df = aucell(exp_mtx, gene_sigs, normalize=False, num_workers=args.n_workers)
        store["perm{}".format(n)] = enrich_df

    store.close()

    # Assign bootstrapped
    print("   * assigning bootstrap results")
    bootstrap_df = assign_bootstrap(args.out_file, n=args.niter, norm=True)
    bootstrap_df.to_csv(args.out_file.split(".h5")[0]+".tsv", sep="\t")

if __name__ == "__main__":
	main()
