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
import timeit
from .transitmodel import *
from .openmp import detect

__all__ = ['test']

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped

def test():
	print("\nStarting tests...\n")
	failures = 0

	params = TransitParams()
	params.t0 = 0.
	params.per = 1.0
	params.rp = 0.1
	params.a = 15.0
	params.inc = 90.
	params.ecc = 0.
	params.w = 90. 
	t = np.linspace(-0.01, 0.01, 1000)

	u = [[], [0.1], [0.1,0.3], [0.1,-0.3], [0.1,-0.03], [0.0, 0.7, 0.0, -0.3]]
	ld = ["uniform", "linear", "quadratic", "exponential", "logarithmic", "nonlinear"]
	
	#Testing different limb darkening models
	for i in range(len(u)):
		params.u = u[i]
		params.limb_dark = ld[i]
		print("Testing "+ld[i]+" limb darkening model...\t")
		m = TransitModel(params, t)
		lc = m.light_curve(params)
		if np.min(lc) <= 1.-params.rp**2: print("\ttest passed")
		else: 
			print("test failed")
			failures += 1

	#Testing truncation error tolerance
	params.u = [0.1, 0.3]
	params.limb_dark = "quadratic"
	m = TransitModel(params, t)
	f_ref = m.light_curve(params)

	params.u = np.array([0.0, 0.7, 0.0, -0.3])
	params.limb_dark = "nonlinear"

	print("\nTesting error tolerance...")
	max_err = np.logspace(1, -2, 4) 
	for i in range(len(max_err)):
		m = TransitModel(params, t, max_err = max_err[i])
		f = m.light_curve(params)
		if np.max(1.0e6*(abs(f-f_ref)) < max_err[i]): print("\t"+str(max_err[i])+ " ppm passed")
		else: 
			print(str(max_err[i]) + " failed")
			failures += 1

	print("\nTesting multithreading...")
	if detect():
		params.u = [0.1, 0.3]
		params.limb_dark = "quadratic"
		m = TransitModel(params, t, nthreads = 1)
		wrapped = wrapper(m.light_curve, params)
		t1 = timeit.timeit(wrapped,number=100)
		m = TransitModel(params, t, nthreads = 2)
		wrapped = wrapper(m.light_curve, params)
		t2 = timeit.timeit(wrapped,number=100)
		if(t1 > 0.9*t2): print("\ttest passed")
		else:
			print("\ttest failed")
			failures +=1
	else: print("\tOpenMP not supported; test ignored")

	if failures == 0: 
		print("\nCongratulations! all tests passed")
		print("""
  ____       ____
  )   \     /   (
   )_  \_V_/  _(
     )__   __(
        `-'""")
	else: 
		print("Uh oh; tests finished with " + "{0}".format(failures) + " failures")
	
		print("""
         *
    __  /|
   /_ \/ |__
  *" \.--._ \           
     ( ^^ )\/
      \__/ *""")
	         
