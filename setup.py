from setuptools import setup, find_packages
import os, sys
#import setuptools 
#from numpy.distutils.core import setup, Extention 

with open("README.md", "r") as fh:
     long_description = fh.read()


if sys.version_info >= (3,7):
    install_requires_py_ver=[
    "numpy>=1.21",
    "numba==0.52",
    "scipy>=1.2.1",
    "matplotlib==3.3.1",
    "formlayout==1.2.0",
    "PyQt5>=5.15.4",
    "jupyter",
    "jupyter-client",
    "ipykernel",
    "qtconsole",
    "pathos>=0.2.5",
    "emcee>=3.0.2",
    "corner>=2.1.0",
    "celerite>=0.3.1",
    "batman-package>=2.4.8",
    "transitleastsquares==1.0.24",
    "dynesty>=1.1",
    "ttvfast>=0.3.0",
    "wotan>=1.7",
#    "wrapt>=1.12.1"
]
else:
    install_requires_py_ver=[
    "numpy>=1.16.6",
    "scipy>=1.2.1",
    "matplotlib==3.2.1",
    "formlayout==1.2.0",
    "PyQt5==5.9.2",
    "jupyter",
    "jupyter-client",
    "ipykernel",
    "qtconsole",
    "pathos>=0.2.5",
    "emcee>=3.0.2",
    "corner>=2.1.0",
    "celerite>=0.3.1",
    "batman-package>=2.4.8",
    "transitleastsquares==1.0.24",
    "dynesty>=1.0.1",
    "ttvfast>=0.3.0",
    "wotan>=1.7"]


setup(
name='exostriker',  
version='0.72',
scripts=['scripts/exostriker'],
author="Trifon Trifonov",
author_email="trifonov@mpia.de",
description="This is the 'Transit and Radial velocity Interactive Fitting tool for Orbital analysis and N-body simulations: The Exo-Striker'",
long_description=long_description,
long_description_content_type="text/markdown",
url="https://github.com/3fon3fonov/exostriker",
#packages=['exostriker'],
packages=find_packages(),
include_package_data = True,
classifiers=[
 "Programming Language :: Python :: 3",
 "License :: MIT License",
 "Operating System :: OS Independent",
],
install_requires=install_requires_py_ver,
extras_requires={
            'pexpect>=4.8.0':['pexpect']},    

license="MIT"
 )
 
 
