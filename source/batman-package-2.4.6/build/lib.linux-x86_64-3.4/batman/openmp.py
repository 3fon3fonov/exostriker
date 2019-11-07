from __future__ import print_function
from distutils.ccompiler import new_compiler
import os
import sys
import tempfile

"""
Check for OpenMP based on
https://github.com/MDAnalysis/mdanalysis/tree/develop/package/setup.py
retrieved 06/15/15
"""

__all__ = ['detect']

def detect():
	"""Does this compiler support OpenMP parallelization?"""
	compiler = new_compiler()
	hasopenmp = hasfunction(compiler, 'omp_get_num_threads()')
	needs_gomp = hasopenmp
	if not hasopenmp:
		compiler.add_library('gomp')
	hasopenmp = hasfunction(compiler, 'omp_get_num_threads()')
	needs_gomp = hasopenmp
	return hasopenmp

def hasfunction(cc, funcname, include=None, extra_postargs=None):
	# From http://stackoverflow.com/questions/
	#            7018879/disabling-output-when-compiling-with-distutils
	tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
	devnull = oldstderr = None
	try:
		try:
			fname = os.path.join(tmpdir, 'funcname.c')
			f = open(fname, 'w')
			if include is not None:
				f.write('#include %s\n' % include)
			f.write('int main(void) {\n')
			f.write('    %s;\n' % funcname)
			f.write('}\n')
			f.close()
			# Redirect stderr to /dev/null to hide any error messages
			# from the compiler.
			devnull = open(os.devnull, 'w')
			oldstderr = sys.stderr.fileno()
			os.dup2(devnull.fileno(), sys.stderr.fileno())
			objects = cc.compile([fname], output_dir=tmpdir, extra_postargs=extra_postargs)
			cc.link_executable(objects, os.path.join(tmpdir, "a.out"))
		except Exception as e:
			return False
		return True
	finally:
		if oldstderr is not None:
			os.dup2(oldstderr, sys.stderr.fileno())
		if devnull is not None:
			devnull.close()

#checks whether OpenMP is supported
