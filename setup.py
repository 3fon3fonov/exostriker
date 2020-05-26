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
packages=['exostriker'],
include_package_data = True,
classifiers=[
 "Programming Language :: Python :: 3",
 "License :: MIT License",
 "Operating System :: OS Independent",
],
install_requires=["numpy","scipy",
"matplotlib",
"PyQt5",
#"PyQt5.QtSvg",
"qtconsole",
"jupyter",
"pathos",
"dill",
"emcee",
"corner",
"celerite",
"transitleastsquares",
"dynesty",
#"rxvt",
"batman",
"ttvfast",
"wotan"],
license="MIT"
 )
