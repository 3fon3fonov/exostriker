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
               No ) echo "WARNING: RVmod/TRIFON will not work without $i!!!"; break;;
           esac
       done
   fi
done
 


#python system install 
arr=( "numpy" "scipy" "matplotlib" "qtconsole" "jupyter" "dill" "emcee" "corner" "celerite")

for i in "${arr[@]}";
do
   if $python -c "import $i" &> /dev/null; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i in your $python?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) conda install $i; break;;
               No ) echo "WARNING: RVmod/TRIFON will not work without $i!!!"; break;;
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
                         $python setup.py install;
                         cp -r ./build/lib*/batman ../../addons/;
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
echo "Installing the swift N-body lib, OK?  (you must if you haven't already!)"
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
       No ) echo "skipped..."; break;;
   esac
done


echo " " 
echo " " 
echo "Compiling the fortran fitting routines, OK? (you must if you haven't already!)"
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
echo "Compiling Symba/mvs and other N-body routines, OK? (you must if you haven't already!)"
echo " " 
echo " " 

select yn in "Yes" "No"; do
   case $yn in
       Yes ) gfortran -O3 ./latest_f/symba_f/swift_symba5_j.f -o  ./stability/symba/swift_symba5_j ./addons/libswift.a;           
             gfortran -O3 ./latest_f/mvs_f/swift_mvs_j.f -o ./stability/mvs/swift_mvs_j ./addons/libswift.a;     
             gfortran -O3 ./latest_f/mvs_f/swift_mvs_j_GR.f -o ./stability/mvs_gr/swift_mvs_j_GR ./addons/libswift.a;
             gfortran -O3 ./latest_f/symba_f/follow_symba2.f -o ./stability/symba/follow_symba2 ./addons/libswift.a;                    
             gfortran -O3 ./latest_f/mvs_f/follow2.f -o ./stability/mvs/follow2 ./addons/libswift.a;
             gfortran -O3 ./latest_f/mvs_f/follow2.f -o ./stability/mvs_gr/follow2 ./addons/libswift.a;                
             gfortran -O3 ./latest_f/symba_f/geninit_j3_in_days.f -o ./stability/symba/geninit_j3_in_days ./addons/libswift.a;   
             gfortran -O3 ./latest_f/mvs_f/geninit_j3_in_days.f -o ./stability/mvs/geninit_j3_in_days ./addons/libswift.a;              
             gfortran -O3 ./latest_f/mvs_f/geninit_j3_in_days.f -o ./stability/mvs_gr/geninit_j3_in_days ./addons/libswift.a;                                    
             break;;
       No ) echo "skipped..."; break;;
   esac
done





 
