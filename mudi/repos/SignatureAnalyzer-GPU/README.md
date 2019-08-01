# SignatureAnalyzer-GPU

# Installation
```
git clone https://github.com/broadinstitute/SignatureAnalyzer-GPU.git
```
To install pytorch please use Anaconda (find more details @ https://pytorch.org/):
```
conda install pytorch torchvision cudatoolkit=9.0 -c pytorch
```

# Setup
For easy set up you can create a python virtual enviroment which matches our own: 
```
$ virtualenv venv

$ source venv/bin/activate .

(venv)$ pip install -r requirements-py3.txt
```

# How to run a single decomposition
SignatureAnalyzer runs on a count matrix (passed to the argument --data) and performs regularized NMF (Bayes NMF). You can specify the regularization you want on the resulting W and H matrices by using the arguments --prior_on_W and --prior_on_H . Passing "L1" is equivalent to an exponential prior and "L2" is half-normal.

For mathematical details see:

1. Tan, V. Y. F., Edric, C.  & Evotte, F. Automatic Relevance Determination in Nonnegative Matrix Factorization with the β-Divergence. (2012). (https://arxiv.org/pdf/1111.6085.pdf)

SignatureAnalyzer-CPU source publications:

1. Kim, J. et al. Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors. Nat. Genet. 48, 600–606 (2016). (https://www.nature.com/articles/ng.3557)

2. Kasar, S. et al. Whole-genome sequencing reveals activation-induced cytidine deaminase signatures during indolent chronic lymphocytic leukaemia evolution. Nat. Commun. 6, 8866 (2015). (https://www.nature.com/articles/ncomms9866)


Note that as part of this work we derived the form for a mixed prior (e.g. L1 on W and L2 on H) see the supplemental note in the repo. 

Example command line for a single run of SignatureAnalyzer-GPU:
```
python SignatureAnalyzer-GPU.py --data example_data/POLEMSI_counts_matrix.txt --max_iter=100000 --output_dir POLEMSI_EXAMPLE --prior_on_W L1 --prior_on_H L2 --labeled
```
Data should be formatted so that the rows are the categories and the columns are the samples. For a full description of inputs and outputs please see the repository wiki. 


# How to run an array of decompositions
The short run time of SignatureAnalyzer-GPU enables performing a parameter search or running the same parameter settings many times to find a maximum likely decomposition or characterize the modal number of clusters/signatures for some setting. To perform such an analysis simply save parameters you would like to run in a tsv and pass it to the --parameters_file argument. We provide the parameters file and count matrix used to generate Figure 1B from the manuscript in the example_data directory. 

NOTE this is automatically configured to run on a single or multiple GPUs just run as usual to perform parallel runs. 
```
python SignatureAnalyzer-GPU.py --data example_data/POLEMSI_counts_matrix.txt --prior_on_W L1 --prior_on_H L2 --output_dir example_data/POLEMSI_outputs/ --parameters_file example_data/POLEMSI_params.txt --max_iter 20000 --labeled --tolerance 1e-7 
```
For a full description of inputs and outputs related to parameters_file runs see the repository wiki. 
