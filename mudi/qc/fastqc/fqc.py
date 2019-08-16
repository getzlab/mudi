import os
from io import StringIO
import subprocess
import zipfile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class FastQC(object):
    """
    FastQC
    ------------------------
    Python class for interpreting FastQC output.
    Reads in .zip file output from fastqc.

    Example usage:
        fqc = FastQC(<FILE_NAME>.zip)
    """

    def __init__(self, dir_path, verbose=False):
        self.path = dir_path
        self.name = '.'.join(dir_path.split('/')[-1].split('.')[:-1])

        self.sample = '_'.join(self.name.split('_')[:2])
        self.lane = self.name.split('_')[-4]

        if verbose: print(self.name,self.sample,self.lane)

        with zipfile.ZipFile(self.path,'r') as z:
            with z.open(os.path.join(self.name,'fastqc_data.txt')) as f:
                self.modules = f.read().decode('utf-8').split('>>END_MODULE')[:-1]

            self.map = {list(pd.read_csv(StringIO(chunk),sep='\t', nrows=1))[0][2:]:idx for idx,chunk in enumerate(self.modules)}
            self.assays = list(self.map.keys())

    def df(self, assay, n=2):
        """
        Query dataframe.
        """
        return pd.read_csv(StringIO(self.modules[self.map[assay]]),sep='\t', skiprows=n)

    @property
    def metrics(self):
        return self.df('FastQC')

    @property
    def seq(self):
        return self.df('Per base sequence quality')

    @property
    def gc(self):
        '''
        Per sequence GC content
        ---------------------------------
        This module measures the GC content across the whole length of each sequence in a file
        and compares it to a modelled normal distribution of GC content.

        In a normal random library you would expect to see a roughly normal distribution of GC content
        where the central peak corresponds to the overall GC content of the underlying genome. Since we
        don't know the the GC content of the genome the modal GC content is calculated from the observed
        data and used to build a reference distribution.

        An unusually shaped distribution could indicate a contaminated library or some other kinds of
        biased subset. A normal distribution which is shifted indicates some systematic bias which is
        independent of base position. If there is a systematic bias which creates a shifted normal
        then this won't be flagged as an error by the module since it doesn't know what your genome's
        GC content should be.

        See here: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html

        '''
        _df = self.df('Per sequence GC content')
        _df['sample'] = self.sample
        _df['lane'] = self.lane
        return _df

    @property
    def dup(self):
        _df = self.df('Sequence Duplication Levels',n=3)
        _df['sample'] = self.sample
        _df['lane'] = self.lane
        return _df

    @property
    def adapt(self):
        _df = self.df('Adapter Content')
        _df['sample'] = self.sample
        _df['lane'] = self.lane
        return _df

    @property
    def overrep_seq(self):
        _df = self.df('Overrepresented sequences')
        _df['sample'] = self.sample
        _df['lane'] = self.lane
        return _df

def aggr_metric(list_of_fqcs, metric, n=2, cols=2):
    """
    Takes in a list of FastQC objects and aggregates them by a given modules.
    This may be used for easy downstream plotting.
    """
    _dfs = list()

    for f in list_of_fqcs:
        _df = f.df(metric,n=n).iloc[:,:cols]
        _df['sample'] = f.sample
        _df['lane'] = f.lane
        _dfs.append(_df)

    return pd.concat(_dfs,0)
