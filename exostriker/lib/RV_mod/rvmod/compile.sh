#!/bin/sh
set -eu
cd "$(dirname "$0")"
FFLAGS="${FFLAGS:-} -O3 -std=legacy -frecursive -I'$(pwd)'" \
    "${PYTHON:-python3}" -m numpy.f2py -c --backend meson \
    -m rvmod_for rvmod_for.f95 --build-dir bdir \
    only: kepfit_amoeba kepfit_lm dynfit_amoeba dynfit_lm rv_capacity_check rv_get_status :
#python -m numpy.f2py -c --opt='-O3 -frecursive' -m rvmod_for rvmod_for.f90
#f2py -c --opt='-O3 -frecursive' -m rvmod rvmod.f
#mv rvmod.cpython-311-x86_64-linux-gnu.so rvmod.so
