# The batman package: fast computation of exoplanet transit light curves
# Copyright (C) 2015 Laura Kreidberg	 
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
import numpy as np
import math
from .transitmodel import *
from . import _quadratic_ld
from . import _nonlinear_ld
import timeit
import batman

__all__ = ['make_plots']

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped


def make_plots():
	import matplotlib.pyplot as plt
	"""zs = np.linspace(0., 1., 1000)
	rp = 0.1
	wrapped = wrapper(_quadratic_ld._quadratic_ld, zs, rp, 0.1, 0.3, 1)
	t = timeit.timeit(wrapped,number=10000)
	print("time:", t)"""
	"""zs = np.linspace(0., 1., 1000)
	rp = 0.1
	u = [0., 0.7, 0.0, -0.3]
	f = _nonlinear_ld._nonlinear_ld(zs, rp, u[0], u[1], u[2], u[3], 1.0e-2, 4)
	fhi = _nonlinear_ld._nonlinear_ld(zs, rp, u[0], u[1], u[2], u[3], 1.0e-4, 4)
	fquad = occultquad.occultquad(zs, rp, 0.1, 0.3, 4)
	#for i in range(len(f)): print "z, fnl, fquad", zs[i], f[i], fquad[i]

	for i in range(1,16):
		wrapped = wrapper(occultquad.occultquad, zs, rp, 0.1, 0.3, i)
		t = timeit.timeit(wrapped,number=1)
		print i, t

	plt.plot(zs, (f - fhi)*1.0e6)
	plt.plot(zs, (fhi - fquad)*1.0e6, color='r')
	plt.axvline(0.9)
	plt.show()"""


	#generates Figure FIXME: max err as a function of function call time
	zs = np.linspace(0., 1., 1000)
	rp = 0.1
	u = [0., 0.7, 0.0, -0.3]
	n = 20
	ts = []
	errs = []
	f_ref = _nonlinear_ld._nonlinear_ld(zs, rp, u[0], u[1], u[2], u[3], 1.0e-4, 4)
	fac = np.logspace(-3, -1, n)
	for i in range(n):
		wrapped = wrapper(_nonlinear_ld._nonlinear_ld, zs, rp, u[0], u[1], u[2], u[3], fac[i], 1)
		t = timeit.timeit(wrapped,number=10)/10.
		ts.append(t)
		print(t)
		f= _nonlinear_ld._nonlinear_ld(zs, rp, u[0], u[1], u[2], u[3], fac[i], 12)
		err = np.max(np.abs(f - f_ref))
		errs.append(err)
	plt.plot(np.array(ts), np.array(errs)*1.0e6, color='k')
	plt.xlim((1.0e-3, 1.0e-1))
	plt.yscale('log')
	plt.xscale('log')
	plt.xlabel("Time (s)")
	plt.ylabel("Max Err (ppm)")
	plt.show()
	

	

