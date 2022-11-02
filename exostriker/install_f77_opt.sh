#!/bin/bash
OPT=$1
echo " " 
echo " " 
echo "Installing the Fortran routines for a first time!"
echo " " 
echo " " 

 
gfortran $OPT  ./source/latest_f/kepfit_amoeba.f -o ./lib/fr/loglik_kep
gfortran $OPT  ./source/latest_f/kepfit_LM.f -o ./lib/fr/chi2_kep 
gfortran $OPT  ./source/latest_f/dynfit_amoeba.f -o ./lib/fr/loglik_dyn
gfortran $OPT  ./source/latest_f/dynfit_LM.f -o ./lib/fr/chi2_dyn
gfortran $OPT  ./source/latest_f/dynfit_amoeba+.f -o ./lib/fr/loglik_dyn+
