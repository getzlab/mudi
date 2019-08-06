from anndata import AnnData
import scanpy as sc
import numpy as np
import pandas as pd
import sys
import os
import torch

# NMF Engine
sig_analyzer_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "repos/SignatureAnalyzer-GPU")
sys.path.append(sig_analyzer_path)
from ARD_NMF import ARD_NMF, run_method_engine

# Relative Imports
from .nmf_utils import compute_phi, nmf_normalize, nmf_scale, nmf_markers

# ---------------------------------
# NMF Wrapper
# ---------------------------------
def NMF(adata, input='raw', use_highly_variable=True, filter_ribo=True, filter_mito=True, \
        K0=None, objective='poisson', max_iter=10000, del_=1, \
        tolerance=1e-6, phi=None, a=10.0, b=None, prior_on_W='L1', prior_on_H='L1', \
        report_frequency=100, active_thresh=1e-5, cut_norm=0.5, cut_diff=1.0,
        inplace=True, verbose=False):
        """
        AnnData wrapper for ARD-NMF.
        ------------------------
        Inputs
            * adata: AnnData object
            * input: either the raw object
            * use_highly_variable: whether or not to use highly variable genes (recommended)
            * filter_ribo: exclude ribosomal genes
            * filter_mito: exclude mitochondrial genes
            * K0: starting number of latent components
            * objective: objective function for optimizaiton
            * max_iter: maximum number of iterations for algorithm
            * del_: n/a
            * tolerance: stop point for optimization
            * phi: dispersion parameter
            * a: shape parameter
            * b: shape parameter
            * prior_on_W: L1 or L2
            * prior_on_H: L1 or L2
            * report_frequency: how often to print stats
            * parameters: parameters file
            * cut_norm
            * cut_diff: difference between mean signature and rest of signatures
                for marker selction
            * active_thresh: threshold for a latent component's impact on signature
                if the latent factor is less than this, it does not contribute
            * inplace: whether or not to edit the AnnData object directly

        Outputs


        It is highly reccomended to use only highly variable genes for factorization.
        Default parameters are to use highly variable genes and filter out mitochondrial
        or ribosomal genes.
        """
        # ---------------------------------
        # Select Genes for NMF
        # ---------------------------------
        if use_highly_variable:
            assert 'highly_variable' in list(adata.var), 'Please compute highly variable genes.'
            if verbose: print("Using {} highly variable genes.".format(sum(adata.var.highly_variable)))
            genes_for_nmf = set(adata.var[adata.var['highly_variable']].index)
        else:
            genes_for_nmf = set(adata.var.index)

        if filter_mito:
            mito_genes = {x for x in genes_for_nmf if x.startswith('MT-')}
            if verbose: print("Filtering {} mito genes.".format(len(mito_genes)))
            genes_for_nmf -= mito_genes
        if filter_ribo:
            ribo_genes = {x for x in genes_for_nmf if x.startswith('RPS') or x.startswith('RPL')}
            if verbose: print("Filtering {} ribo genes.".format(len(ribo_genes)))
            genes_for_nmf -= ribo_genes

        # ---------------------------------
        # Inputs
        # ---------------------------------
        if input == 'raw':
            df = pd.DataFrame(data=adata.raw.X.toarray(), index=adata.obs_names, columns=adata.raw.var_names).T
        elif input == 'count':
            assert adata.layers['counts'], 'Please save raw counts in adata.lyers["counts"]'
            df = pd.DataFrame(data=adata.layers['counts'].toarray(), index=adata.obs_names, columns=adata.var_names).T

        df = df.loc[genes_for_nmf]

        if objective == 'poisson':
            Beta = 1
        elif objective == 'gaussian':
            Beta = 2
        else:
            ValueError("Objective should be either 'gaussian' or 'poisson'")

        data = ARD_NMF(df, objective)
        channel_names = data.channel_names
        sample_names = data.sample_names

        if phi is None:
            phi = compute_phi(np.mean(df.values), np.var(df.values), Beta)

        # ---------------------------------
        # Run NMF
        # ---------------------------------
        W, H, cost = run_method_engine(
            data, \
            a, \
            phi, \
            b, \
            Beta, \
            prior_on_W, \
            prior_on_H, \
            K0, \
            tolerance, \
            max_iter \
        )

        W, H, nsig = nmf_normalize(W, H, active_thresh=active_thresh)
        sig_names = [str(i) for i in range(1,nsig+1)]
        W = pd.DataFrame(data=W, index=channel_names, columns=sig_names)
        H = pd.DataFrame(data=H, index=sig_names, columns=sample_names)

        W,H = nmf_scale(W,H)
        markers, gene_signatures = nmf_markers(df,W,H,cut_norm=cut_norm,cut_diff=cut_diff)

        if inplace:
            adata.obs = adata.obs.join(H)
            adata.obs['max_id'] = adata.obs['max_id'].astype('category')
            adata.obs = adata.obs.rename(columns={'max_id':'nmf'})

            adata.uns['signatures'] = list(H)[:-3]
            adata.obsm['X_nmf'] = H.iloc[:,:-3].values

            _nmf_genes = {}
            _nmf_genes['cols'] = list(gene_signatures)
            _nmf_genes['rows'] = gene_signatures.index
            _nmf_genes['values'] = gene_signatures.values
            adata.uns['nmf_genes'] = _nmf_genes

            _nmf_markers = {}
            _nmf_markers['cols'] = list(markers)
            _nmf_markers['rows'] = markers.index
            _nmf_markers['values'] = markers.values
            adata.uns['nmf_markers'] = _nmf_markers
            return
        else:
            return H,W,markers,gene_signatures



# # ---------------------------------
# # NMF Wrapper Class
# # ---------------------------------
# class NMF(object):
#     """
#     AnnData wrapper for ARD-NMF.
#
#     Inputs
#         - adata: AnnData object
#         - K0: starting number of latent components
#         - objective: objective function for optimizaiton
#         - max_iter: maximum number of iterations for algorithm
#         - del_: n/a
#         - tolerance: stop point for optimization
#         - phi: dispersion parameter
#         - a
#         - b
#         - prior_on_W
#         - prior_on_H
#
#
#
#     It is highly reccomended to use only highly variable genes for factorization.
#     Default parameters are to use highly variable genes and filter out mitochondrial
#     or ribosomal genes.
#
#
#     """
#
#     def __init__(self, adata, use_highly_variable=True, filter_ribo=True, filter_mito=True, raw=False, verbose=False, K0=None, objective='poisson', max_iter=10000, del_=1, tolerance=1e-6, phi=None, a=10.0, b=None, prior_on_W='L1', prior_on_H='L1', report_frequency=100, parameters=None, dtype='Float32', active_thresh=1e-5, num_processes='auto'):
#         """
#         Parse input to run NMF.
#
#         """
#         if use_highly_variable:
#             if verbose: print("Using {} highly variable genes.".format(sum(adata.var.highly_variable)))
#             genes_for_nmf = set(adata.var[adata.var['highly_variable']].index)
#         else:
#             genes_for_nmf = set(adata.var.index)
#
#         if filter_mito:
#             mito_genes = {x for x in adata.var_names if x.startswith('MT-')}
#             if verbose: print("Filtering {} mito genes.".format(len(mito_genes)))
#             genes_for_nmf -= mito_genes
#         if filter_ribo:
#             ribo_genes = {x for x in adata.var_names if x.startswith('RPS') or x.startswith('RPL')}
#             if verbose: print("Filtering {} ribo genes.".format(len(ribo_genes)))
#             genes_for_nmf -= ribo_genes
#
#         #print(adata)
#         return
#
#
#
#
#
#         self.X = X
#         self.objective = objective
#         self.parameters = parameters
#         self.mu = np.mean(self.X.values)
#         self.var = np.var(self.X.values)
#
#         if dtype == 'Float32':
#             self.dtype = torch.float32
#         elif dtype == 'Float16':
#             self.dtype = torch.float16
#
#         if self.objective == 'poisson':
#             self.Beta = 1
#         elif self.objective == 'gaussian':
#             self.Beta = 2
#         else:
#             ValueError("Objective should be either 'gaussian' or 'poisson'")
#
#         start_alg_time = time.time()
#
#         data = ARD_NMF(self.X, self.objective)
#
#         self.channel_names = data.channel_names
#         self.sample_names = data.sample_names
#
#         if phi is None:
#             print("Computing dispersion parameter.")
#             self.phi = self.compute_phi()
#         else:
#             self.phi = phi
#
#         self.W, self.H, self.cost = run_method_engine(data, a, self.phi, b, self.Beta, prior_on_W, prior_on_H, K0, tolerance, max_iter)
#         self.alg_time = time.time() - start_alg_time
#
#         self.W, self.H, self.nsig = self.normalize_nmf_output()
#         self.W, self.H = self.nmf_to_pd()
#
#     def compute_nmf_result(self, cut_norm=0.5, cut_diff=1.0):
#         """Assigns NMF_Result object to result attribute."""
#         self.result = NMF_Result(self.X, self.W, self.H, cut_norm=cut_norm, cut_diff=cut_diff, objective=self.objective)
#
#     def compute_phi(self):
#         """
#         Compute dispersion parameter (phi).
#         """
#         return self.var / (self.mu ** (2 - self.Beta))
#
#     def normalize_nmf_output(self, active_thresh=1e-5):
#         """
#         Prunes output from ARD-NMF.
#         """
#         nonzero_idx = (np.sum(self.H, axis=1) * np.sum(self.W, axis=0)) > active_thresh
#         W_active = self.W[:, nonzero_idx]
#         H_active = self.H[nonzero_idx, :]
#         nsig = np.sum(nonzero_idx)
#
#         # Normalize W and transfer weight to H matrix
#         W_weight = np.sum(W_active, axis=0)
#         W_final = W_active / W_weight
#         H_final = W_weight[:, np.newaxis] * H_active
#
#         return W_final, H_final, nsig
#
#     def nmf_to_pd(self):
#         """
#         Collapse NMF output to pandas dataframe.
#         """
#         sig_names = [str(i) for i in range(1,self.nsig+1)]
#         return pd.DataFrame(data=self.W, index=self.channel_names, columns=sig_names), pd.DataFrame(data=self.H, index=sig_names, columns=self.sample_names)
