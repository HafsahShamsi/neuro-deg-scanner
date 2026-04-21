from setuptools import setup, find_packages

setup(
    name             = 'neuro-deg-scanner',
    version          = '1.0.0',
    author           = 'Hafsah Shamsi',
    description      = 'Differential gene expression pipeline for neuroimmune disease GEO datasets',
    long_description = open('README.md').read(),
    long_description_content_type = 'text/markdown',
    packages         = find_packages(),
    python_requires  = '>=3.8',
    install_requires = [
        'pandas>=2.0',
        'numpy>=1.26',
        'matplotlib>=3.7',
        'seaborn>=0.13',
        'scipy>=1.11',
        'scikit-learn>=1.3',
        'statsmodels>=0.14',
        'GEOparse>=2.0',
    ],
    extras_require = {
        'labels': ['adjustText>=1.0'],
    },
    classifiers = [
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
