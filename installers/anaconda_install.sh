#!/bin/bash
## Check your Debian based system if is ready for RVmod/TRIFON. 


echo "For which anaconda version we want to check/install packages?"
select py in "anaconda2" "anaconda3"; do
   case $py in
       anaconda2 ) python="python2" ; break;;
       anaconda3 ) python="python3"; break;;
   esac
done 


conda config --add channels conda-forge

#system needed
arr=( "gfortran")

for i in "${arr[@]}";
do
   if type $i >/dev/null 2>&1; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Please install $i version <= 7 and try again"; 
       exit 1
   fi
done


#system needed
arr=( "csh")

for i in "${arr[@]}";
do
   if type $i >/dev/null 2>&1; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Please install $i and try again"; 
       exit 1
   fi
done


#system optional
arr=( "rxvt" )

for i in "${arr[@]}";
do
   if type $i >/dev/null 2>&1; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Please install the $i bash shell for better experience (xterm used by default)"; 
   fi
done



#python system install 
arr=( "PyQt5" )

for i in "${arr[@]}";
do
   if $python -c "import $i" &> /dev/null; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i in your $python?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) conda install pyqt; break;;
               No ) echo "WARNING: 'The Exo-Striker' will not work without $i!!!"; break;;
           esac
       done
   fi
done
 


#python system install 
arr=( "numpy" "scipy" "matplotlib" "qtconsole" "jupyter" "pathos" "dill" "emcee" "corner" "celerite" "dynesty" "ttvfast")

for i in "${arr[@]}";
do
   if $python -c "import $i" &> /dev/null; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i in your $python?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) conda install $i; break;;
               No ) echo "WARNING: 'The Exo-Striker' will not work without $i!!!"; break;;
           esac
       done
   fi
done


 


#python local install (under test)
arr=( "batman")

for i in "${arr[@]}";
do
   if $python -c "import $i" &> /dev/null; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Only a local instalation possible... Do you wish a install $i from source (from ./source directory)?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) if [ $i=="batman" ]; then 
                         echo "Installing $i from source ./source ---> ./lib  ...";
                         cd source/batman-package-2.4.6/;
                         $python setup.py install;
                         cp -r ./build/lib*/batman ../../lib/;
                         rm -r ./build
                         cd ../../;
                     else 
                         echo "TBD for another package!";
                     fi
                     break;;          

               No ) echo "WARNING: 'The Exo-Striker' will not work without $i!!!"; break;;

           esac
       done
   fi
done


echo " " 
echo " " 
echo "Installing the swift N-body lib, OK?  (you must, if you haven't done it already!)"
echo " " 
echo " " 

select yn in "Yes" "No"; do
   case $yn in
       Yes ) cd source/swift_j/;
             awk -v a="$PWD" '{ if (NR == 3) print "set SWIFT_DIR="a; else print $0}' @make_temp > @make
             csh @makeall;
             cp libswift.a ../../lib/;
             cd ../../;         
             break;;
       No ) echo "skipped..."; break;;
   esac
done


echo " " 
echo " " 
echo "Compiling the fortran fitting routines, OK? (you must, if you haven't done it already!)"
echo " " 
echo " " 

select yn in "Yes" "No"; do
   case $yn in
       Yes ) gfortran -O3 ./source/latest_f/kepfit_LM.f -o ./lib/fr/chi2_kep; # chi2 keplerian
             gfortran -O3 ./source/latest_f/dynfit_LM.f -o ./lib/fr/chi2_dyn; # chi2 dynamical
             gfortran -O3 ./source/latest_f/kepfit_amoeba.f -o ./lib/fr/loglik_kep; # lnL keplerian
             gfortran -O3 ./source/latest_f/dynfit_amoeba.f -o ./lib/fr/loglik_dyn; # lnL dynamical               
             gfortran -O3 ./source/latest_f/dynfit_amoeba+.f -o ./lib/fr/loglik_dyn+; # lnL dynamical/keplerian mixed                 
             break;;
       No ) echo "skipped..."; break;;
   esac
done


echo " " 
echo " " 
echo "Compiling Symba/mvs and other N-body routines, OK? (you must, if you haven't done it already!)"
echo " " 
echo " " 

select yn in "Yes" "No"; do
   case $yn in
       Yes ) gfortran -O3 ./source/latest_f/symba_f/swift_symba5_j.f -o  ./stability/symba/swift_symba5_j ./lib/libswift.a;           
             gfortran -O3 ./source/latest_f/mvs_f/swift_mvs_j.f -o ./stability/mvs/swift_mvs_j ./lib/libswift.a;     
             gfortran -O3 ./source/latest_f/mvs_f/swift_mvs_j_GR.f -o ./stability/mvs_gr/swift_mvs_j_GR ./lib/libswift.a;
             gfortran -O3 ./source/latest_f/symba_f/follow_symba2.f -o ./stability/symba/follow_symba2 ./lib/libswift.a;                    
             gfortran -O3 ./source/latest_f/mvs_f/follow2.f -o ./stability/mvs/follow2 ./lib/libswift.a;
             gfortran -O3 ./source/latest_f/mvs_f/follow2.f -o ./stability/mvs_gr/follow2 ./lib/libswift.a;                
             gfortran -O3 ./source/latest_f/symba_f/geninit_j3_in_days.f -o ./stability/symba/geninit_j3_in_days ./lib/libswift.a;   
             gfortran -O3 ./source/latest_f/mvs_f/geninit_j3_in_days.f -o ./stability/mvs/geninit_j3_in_days ./lib/libswift.a;              
             gfortran -O3 ./source/latest_f/mvs_f/geninit_j3_in_days.f -o ./stability/mvs_gr/geninit_j3_in_days ./lib/libswift.a;                                    
             break;;
       No ) echo "skipped..."; break;;
   esac
done





 
