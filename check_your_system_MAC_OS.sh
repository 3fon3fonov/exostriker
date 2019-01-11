#!/bin/bash
## Check your Mac OS based system if is ready for RVmod/TRIFON. 


echo "For which Python version we want to check/install packages?"
select py in "Python2" "Python3"; do
   case $py in
       Python2 ) python="python" ; pip="pip" ; break;;
       Python3 ) python="python3"; pip="pip3"; break;;
   esac
done  

#system needed
arr=( "gfortran")

for i in "${arr[@]}";
do
   if type $i >/dev/null 2>&1; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) brew install gcc; break;;
               No ) echo "WARNING: RVmod/TRIFON will not work without $i!!!"; break;;
           esac
       done     
   fi
done

arr=( "csh")

for i in "${arr[@]}";
do
   if type $i >/dev/null 2>&1; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) brew install tcsh; break;;
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
               Yes ) brew install pyqt5; break;;
               No ) echo "WARNING: RVmod/TRIFON will not work without $i!!!"; break;;
           esac
       done
   fi
done

#python system install  
arr=( "PyQt5.QtSvg" )

for i in "${arr[@]}";
do
   if $python -c "import $i" &> /dev/null; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) brew install pyqt5.qtsvg; break;;
               No ) echo "WARNING: RVmod/TRIFON will not work without $i!!!"; break;;
           esac
       done
   fi
done




#python system install 
arr=( "setuptools" "numpy" "scipy" "matplotlib")

for i in "${arr[@]}";
do
   if $python -c "import $i" &> /dev/null; then
       echo "$i - yes!"
   else
       echo "$i - not installed! Do you wish to install $i?"
       select yn in "Yes" "No"; do
           case $yn in
               Yes ) $pip install $i --user; break;;
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
               Yes ) $pip install $i --user; break;;
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
                         $pip install --user --install-option="--prefix=" $i-package;
                     else 
                         echo "TBD for another package!";
                     fi
                     break;;   

               Yes-Local ) if [ $i=="batman" ]; then 
                         echo "We are installing $i from source ./deps ---> ./addons  ...";
                         cd deps/batman-package-2.4.6/;
                         $python setup.py install --user --prefix=;
                         #sudo python setup.py install;
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
       No ) echo "skiped..."; break;;
   esac
done


#system optional
#arr=( "urxvt" )

#for i in "${arr[@]}";
#do
#   if type $i >/dev/null 2>&1; then
#       echo "$i - yes!"
#   else
#       echo "$i - not installed! Do you wish to install $i?"
#       select yn in "Yes" "No"; do
#           case $yn in
#               Yes ) brew cask install xquartz; brew install rxvt-unicode; break;;
#               No ) echo "WARNING: RVmod/TRIFON may not work properly without $i!!!"; break;;
#           esac
#       done  
#   fi
#done



