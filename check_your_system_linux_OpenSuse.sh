#!/bin/bash
## Check your Debian based system if is ready for RVmod/TRIFON. 


echo "For which Python version we want to check/install packages?"
select py in "Python2" "Python3"; do
   case $py in
       Python2 ) python="python2" ; pip="pip2" ; break;;
       Python3 ) python="python3"; pip="pip3"; break;;
   esac
done  


#system needed
arr=( "gfortran" "csh")

for i in "${arr[@]}";
do
   if type $i >/dev/null 2>&1; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) sudo zypper install gcc-fortran; break;;
               No ) echo "WARNING: RVmod/TRIFON will not work without $i!!!"; break;;
           esac
       done     
   fi
done



#python system install 
arr=( "setuptools" "pip" "numpy" "scipy" "matplotlib")

for i in "${arr[@]}";
do
   if $python -c "import $i" &> /dev/null; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i?"
       select yn in "Yes" "No"; do
           case $yn in
#               Yes ) sudo zypper install $python-$i; break;;
               Yes ) $python -m install --user $i; break;;
               No ) echo "WARNING: RVmod/TRIFON will not work without $i!!!"; break;;
           esac
       done
   fi
done


#python system install 
arr=( "PyQt5" )

for i in "${arr[@]}";
do
   if $python -c "import $i" &> /dev/null; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) sudo zypper install $python-qt5; break;;
               No ) echo "WARNING: RVmod/TRIFON will not work without $i!!!"; break;;
           esac
       done
   fi
done

#python system install (sudo apt install)
arr=("PyQt5.QtSvg")

for i in "${arr[@]}";
do
   if $python -c "import $i" &> /dev/null; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) sudo zypper install $python-qt5.qtsvg; break;;
               No ) echo "WARNING: RVmod/TRIFON will not work without $i!!!"; break;;
           esac
       done
   fi
done




#python pip install
arr=( "qtconsole" "jupyter" "dill" "emcee" "corner" "celerite")

for i in "${arr[@]}";
do
   if $python -c "import $i" &> /dev/null; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) sudo $pip install $i; break;;
               No ) echo "WARNING: RVmod/TRIFON will not work without $i!!!"; break;;
           esac
       done
   fi
done


# odd python isnstall 




#system optional
arr=( "rxvt" )

for i in "${arr[@]}";
do
   if type $i >/dev/null 2>&1; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) sudo zypper install rxvt-unicode; break;;
               No ) echo "WARNING: RVmod/TRIFON may not work properly without $i!!!"; break;;
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
       echo "$i - not installed! Do you wish a install $i from source (from ./deps directory)?"
       select yn in "Yes" "Yes-Local" "No"; do
           case $yn in
               Yes ) if [ $i=="batman" ]; then 
                         sudo $pip install $i-package;
                     else 
                         echo "TBD for another package!";
                     fi
                     break;;   

               Yes-Local ) if [ $i=="batman" ]; then 
                         echo "We are installing $i from source ./deps ---> ./addons  ...";
                         cd deps/batman-package-2.4.6/;
                         sudo $python setup.py install;
                         cp -r ./build/lib.linux-x86_64-2.7/batman ../../addons/;
                         cd ../../;
                     else 
                         echo "TBD for another package!";
                     fi
                     break;;          

               No ) echo "WARNING: RVmod/TRIFON will not work without $i!!!"; break;;

           esac
       done
   fi
done


echo " " 
echo " " 
echo "Installing the swift N-body lib, OK?  (you must if you havent already!)"
echo " " 
echo " " 

select yn in "Yes" "No"; do
   case $yn in
       Yes ) cd deps/swift_j/;
             awk -v a="$PWD" '{ if (NR == 3) print "set SWIFT_DIR="a; else print $0}' @make_temp > @make
             csh @makeall;
             cp libswift.a ../../addons/;
             cd ../../;         
             break;;
       No ) echo "skiped..."; break;;
   esac
done


echo " " 
echo " " 
echo "Compiling the fortran fitting routines, OK? (you must if you havent already!)"
echo " " 
echo " " 

select yn in "Yes" "No"; do
   case $yn in
       Yes ) gfortran -O3 ./latest_f/kepfit_LM_v1b.f -o ./fitting_routines/chi2_kep ./addons/libswift.a; # chi2 keplerian
             gfortran -O3 ./latest_f/dynfit_LM_v1b.f -o ./fitting_routines/chi2_dyn ./addons/libswift.a; # chi2 dynamical
             gfortran -O3 ./latest_f/kepfit_amoeba_v1b.f -o ./fitting_routines/loglik_kep ./addons/libswift.a; # lnL keplerian
             gfortran -O3 ./latest_f/dynfit_amoeba_v1b.f -o ./fitting_routines/loglik_dyn ./addons/libswift.a; # lnL dynamical               
             break;;
       No ) echo "skiped..."; break;;
   esac
done


echo " " 
echo " " 
echo "Compiling Symba/mvs and other N-body routines, OK? (you must if you havent already!)"
echo " " 
echo " " 

select yn in "Yes" "No"; do
   case $yn in
       Yes ) gfortran -O3 ./latest_f/symba_f/follow_symba2.f -o ./stability/symba/follow_symba2 ./addons/libswift.a; # chi2 keplerian
             gfortran -O3 ./latest_f/symba_f/geninit_j3_in_days.f -o ./stability/symba/geninit_j3_in_days ./addons/libswift.a; # chi2 dynamical
             gfortran -O3 ./latest_f/symba_f/swift_symba5_j.f -o  ./stability/symba/swift_symba5_j ./addons/libswift.a; # lnL keplerian
             break;;
       No ) echo "skiped..."; break;;
   esac
done










# Not working, not sure why...


#python system install (sudo apt install)
#arr=( "setuptools" "pip" "numpy" "scipy" "matplotlib" "PyQt5" "PyQt5.QtSvg")

#for i in "${arr[@]}";
#do
#   if $python -c "import $i" &> /dev/null; then
#       echo "$i - yes!"
#   else
#       echo "$i - not installed! Do you wish to install $i?"
#       select yn in "Yes" "No"; do
#           case $yn in
#               Yes ) if [ $i=="PyQt5" ]; then 
#                         sudo apt install $python-pyqt5;
#                     elif [ $i=="PyQt5.QtSvg" ]; then 
#                         sudo apt install $python-pyqt5.qtsvg; 
#                     else
#                         sudo apt install $python-$i; 
#                     fi
#                     break;;
#               No ) echo "WARNING: RVmod/TRIFON will not work without $i!!!"; break;;
#           esac
#       done
#   fi
#done






