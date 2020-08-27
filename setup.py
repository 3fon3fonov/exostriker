from setuptools import setup, find_packages
import os
import setuptools 
#from numpy.distutils.core import setup, Extention 

with open("README.md", "r") as fh:
     long_description = fh.read()

setup(
name='exostriker',  
version='0.40',
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
"numpy>=1.16.6",
"scipy>=1.2.1",
"matplotlib==3.2.1",
"PyQt5==5.9.2",
"qtconsole",
"jupyter",
"jupyter-client",
"ipykernel",
"dill>=0.3.1",
"pathos>=0.2.5",
"emcee>=3.0.2",
"corner>=2.1.0",
"celerite>=0.3.1",
"transitleastsquares==1.0.24",
"dynesty>=1.0.1",
"ttvfast>=0.3.0",
"wotan>=1.7"],
extra_requires=['pexpect>=4.8.0'],
license="MIT"
 )
 
 
