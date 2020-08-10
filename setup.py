#!usr/bin/env python

#from distutils.core import setup
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

setup(
    name='dxm',
    version='0.1.1',
    author='Jerry Fong',
    author_email='fongj@wustl.edu',
    maintainer='John Edwards',
    maintainer_email='jredwards@wustl.edu',
    url='',
    license='GPLv3',
    description='Deconvolution of Methylation Sequencing Data',
    long_description='DXM will deconvolve the processed methylation sequencing data from a heterogeneous population into its subpopulation methylation profiles and their relative prevalence.',
    install_requires=['numpy','cython'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: POSIX',
        'Topic :: Science/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.5'
        ],
    packages=find_packages(),
    entry_points = {
        'console_scripts': [
            'dxm_solveMethylation = dxm.dxm_solveMethylation:main',
            'dxm_callIDMR = dxm.dxm_callIDMR:main',
            'dxm_estimateFracs = dxm.dxm_estimateFracs:main',
            ],
        },
    include_package_data=True,
    data_files=[('dxm',['dxm/Bins.txt','dxm/Freqs.txt','dxm/DXMfunctions.pyx'])],
    ext_modules = cythonize("dxm/DXMfunctions.pyx"), include_dirs=[numpy.get_include()],
    zip_safe=False,
    python_requires='>=3.5',
)
