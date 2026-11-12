#!/bin/sh
set -eu
cd "$(dirname "$0")"
FFLAGS="${FFLAGS:-} -O3 -std=legacy -fmax-stack-var-size=32768 -I'$(pwd)'" \
    "${PYTHON:-python3}" -m numpy.f2py -c --backend meson \
    -m rvmod_for rvmod_for.f95 joint_fit.f95 shared_nbody.f95 model_only.f95 --build-dir bdir \
    only: kepfit_amoeba kepfit_lm dynfit_amoeba dynfit_lm rv_capacity_check rv_get_status rv_transit_model rv_transit_nbody rv_transit_model_parts rv_transit_nbody_parts rv_transit_phased_model rv_transit_density rv_transit_geometry timing_model astrometry_model astrometry_nbody transit_event_model :
#python -m numpy.f2py -c --opt='-O3 -frecursive' -m rvmod_for rvmod_for.f90
#f2py -c --opt='-O3 -frecursive' -m rvmod rvmod.f
#mv rvmod.cpython-311-x86_64-linux-gnu.so rvmod.so
