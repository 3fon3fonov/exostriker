#!/bin/bash
 
echo " " 
echo " " 
echo "Installing the swift N-body lib for a first time!"
echo " " 
echo " " 

cd source/swift_j/;
awk -v a="$PWD" '{ if (NR == 3) print "set SWIFT_DIR="a; else print $0}' @make_temp > @make
csh @makeall;
cp libswift.a ../../lib/;
cd ../../;


echo " " 
echo " " 
echo "Compiling Symba/mvs and other N-body routines for a first time!"
echo " " 
echo " " 

gfortran -O3 ./source/latest_f/symba_f/swift_symba5_j.f -o  ./stability/symba/swift_symba5_j ./lib/libswift.a;
gfortran -O3 ./source/latest_f/mvs_f/swift_mvs_j.f -o ./stability/mvs/swift_mvs_j ./lib/libswift.a;
gfortran -O3 ./source/latest_f/mvs_f/swift_mvs_j_GR.f -o ./stability/mvs_gr/swift_mvs_j_GR ./lib/libswift.a;
gfortran -O3 ./source/latest_f/symba_f/follow_symba2.f -o ./stability/symba/follow_symba2 ./lib/libswift.a;
gfortran -O3 ./source/latest_f/mvs_f/follow2.f -o ./stability/mvs/follow2 ./lib/libswift.a;
gfortran -O3 ./source/latest_f/mvs_f/follow2.f -o ./stability/mvs_gr/follow2 ./lib/libswift.a;
gfortran -O3 ./source/latest_f/symba_f/geninit_j3_in_days.f -o ./stability/symba/geninit_j3_in_days ./lib/libswift.a;
gfortran -O3 ./source/latest_f/mvs_f/geninit_j3_in_days.f -o ./stability/mvs/geninit_j3_in_days ./lib/libswift.a;
gfortran -O3 ./source/latest_f/mvs_f/geninit_j3_in_days.f -o ./stability/mvs_gr/geninit_j3_in_days ./lib/libswift.a;
 


 




