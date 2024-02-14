python3.11 -m numpy.f2py -c --opt='-O3 -frecursive --std=f95' -m rvmod_for rvmod_for.f95
#python -m numpy.f2py -c --opt='-O3 -frecursive' -m rvmod_for rvmod_for.f90
#f2py -c --opt='-O3 -frecursive' -m rvmod rvmod.f
#mv rvmod.cpython-311-x86_64-linux-gnu.so rvmod.so
