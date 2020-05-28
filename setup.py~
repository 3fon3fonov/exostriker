from setuptools import setup, find_packages
import os
import setuptools 
#from numpy.distutils.core import setup, Extention 

with open("README.md", "r") as fh:
     long_description = fh.read()

#os.chdir("exostriker")

#os.system("echo hallo world!")
#os.system("""
#gfortran -O3 ./source/latest_f/kepfit_LM.f -o ./lib/fr/chi2_kep; # chi2 keplerian
#gfortran -O3 ./source/latest_f/dynfit_LM.f -o ./lib/fr/chi2_dyn; # chi2 dynamical
#gfortran -O3 ./source/latest_f/kepfit_amoeba.f -o ./lib/fr/loglik_kep; # lnL keplerian
#gfortran -O3 ./source/latest_f/dynfit_amoeba.f -o ./lib/fr/loglik_dyn; # lnL dynamical
#gfortran -O3 ./source/latest_f/dynfit_amoeba+.f -o ./lib/fr/loglik_dyn+; # lnL dynamical/keplerian mixed
#""")

#os.chdir("../")
 

setup(
name='exostriker',  
version='0.27',
scripts=['scripts/exostriker'],
author="Trifon Trifonov",
author_email="trifonov@mpia.de",
description="This is The Exo-Striker",
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
install_requires=[
"numpy>=1.18.4",
"scipy>=1.2.1",
"matplotlib==3.2.1",
"PyQt5==5.9.2",
#"PyQt5.QtSvg",
"qtconsole>=4.5.5",
"jupyter==1.0.0",
"jupyter-client==5.3.4",
"ipykernel==5.1.3",
"pathos>=0.2.5",
"dill",
"emcee>=3.0.2",
"corner",
"celerite>=0.3.1",
"transitleastsquares>=1.0.24",
"dynesty>=1.0.1",
"ttvfast>=0.3.0",
"wotan>=1.7"],
license="MIT"
 )
 
 
