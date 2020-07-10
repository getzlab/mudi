from setuptools import setup
import re
import os
import sys

ver_info = sys.version_info
if ver_info < (3,6,0):
    raise RuntimeError("signatureanalyzer requires at least python 3.6.0")

with open(os.path.join(os.path.dirname(__file__), 'mudi', '__init__.py')) as r:
    version = re.search(r'__version__ = \'(\d+\.\d+\.\d+[-_a-zA-Z0-9]*)\'', r.read()).group(1)

setup(
    name='mudi',
    python_requires='>3.6.0',
    author='Shankara Anand - Broad Institute - Cancer Genome Computational Analysis',
    author_email='sanand@broadinstitute.org',
    url = 'https://github.com/broadinstitute/mudi',
    long_description = open("README.md", encoding="utf-8").read(),
    long_description_content_type = 'text/markdown',
    description='mudi: a collection of single-cell utilities. (Getz Lab).',
    packages = [
        'mudi',
        'mudi.plotting',
        'mudi.diffexp',
        'mudi.diffexp.glmm',
        'mudi.nmf',
        'mudi.norm',
        'mudi.qc',
        'mudi.qc.fastqc',
        'mudi.integration',
        'mudi.interp',
    ],
    install_requires=[
        "bbknn>=1.3.4",
        #"cellxgene>=0.10.1",
        "anndata>=0.7.1",
        "scanpy>=1.4.3",
        "rpy2>=3.0.5",
        "torch>=1.1.0",
        "scrublet",
        "openpyxl",
        "gprofiler",
        "tzlocal",
        "get_version",
        "signatureanalyzer",
        "anndata2ri",
        "ipywidgets", # Added for scanpy,
        "agutil",
        "tqdm"
    ],
    package_data = {
    "":[
        "ref/cell_cycle_genes/*",
        "ref/markers/*"
    ],
    },
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT"
)
