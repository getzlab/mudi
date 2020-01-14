from setuptools import setup

setup(
    name='mudi',
    python_requires='>3.6.0',
    author='Shankara Anand',
    author_email='sanand@broadinstitute.org',
    description='MUDI: a collection of single-cell utilities. (Getz Lab).',
    install_requires=[
        "bbknn>=1.3.4",
        "cellxgene>=0.10.1",
        "anndata>=0.6.21",
        "scanpy>=1.4.3",
        "rpy2>=3.0.5",
        "torch>=1.1.0"
        "scrublet",
        "multicoreTSNE",
        "openpyxl",
        "gprofiler",
        "tzlocal",
        "get_version",
        "signatureanalyzer"
    ]
)
