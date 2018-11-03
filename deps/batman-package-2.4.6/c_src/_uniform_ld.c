/* The batman package: fast computation of exoplanet transit light curves
 * Copyright (C) 2015 Laura Kreidberg	 
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include "numpy/arrayobject.h"

#if defined (_OPENACC) && defined(__PGI)
#  include <accelmath.h>
#else
#  include <math.h>
#endif

#if defined (_OPENMP) && !defined(_OPENACC)
#  include <omp.h>
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


static PyObject *_uniform_ld(PyObject *self, PyObject *args)
{
	int nthreads;
	double p;

	PyArrayObject *ds, *flux;
	npy_intp dims[1];
	
  	if(!PyArg_ParseTuple(args, "Odi", &ds, &p, &nthreads)) return NULL;		//parses function input

	dims[0] = PyArray_DIMS(ds)[0]; 
	flux = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(ds));	//creates numpy array to store return flux values
	
	double *f_array = PyArray_DATA(flux);
	double *d_array = PyArray_DATA(ds);

	if(fabs(p - 0.5) < 1.e-3) p = 0.5;

	#if defined (_OPENMP) && !defined(_OPENACC)
	omp_set_num_threads(nthreads);	//specifies number of threads (if OpenMP is supported)
	#endif

	#if defined (_OPENACC)
	#pragma acc parallel loop copyout(f_array[:dims[0]])
	#elif defined (_OPENMP)
	#pragma omp parallel for
	#endif
	for(int i=0; i<dims[0]; i++)
	{
		double d = d_array[i]; 				// separation of centers
		
		if(d >= 1. + p) f_array[i] = 1.;		//no overlap
		if(p >= 1. && d <= p - 1.) f_array[i] = 0.;	//total eclipse of the star
		else if(d <= 1. - p) f_array[i] = 1. - p*p;	//planet is fully in transit
		else						//planet is crossing the limb
		{
			double kap1=acos(fmin((1. - p*p + d*d)/2./d, 1.));
			double kap0=acos(fmin((p*p + d*d - 1.)/2./p/d, 1.));
			f_array[i] = 1. - (p*p*kap0 + kap1 - 0.5*sqrt(fmax(4.*d*d - pow(1. + d*d - p*p, 2.), 0.)))/M_PI;
		}
	}

	return PyArray_Return((PyArrayObject *)flux);
}


static char _uniform_ld_doc[] = "This extension module returns a limb darkened light curve for a uniform stellar intensity profile.";

static PyMethodDef _uniform_ld_methods[] = {
  {"_uniform_ld", _uniform_ld, METH_VARARGS, _uniform_ld_doc},{NULL}};

#if PY_MAJOR_VERSION >= 3
	static struct PyModuleDef _uniform_ld_module = {
		PyModuleDef_HEAD_INIT,
		"_uniform_ld",
		_uniform_ld_doc,
		-1, 
		_uniform_ld_methods
	};

	PyMODINIT_FUNC
	PyInit__uniform_ld(void)
	{
		PyObject* module = PyModule_Create(&_uniform_ld_module);
		if(!module)
		{
			return NULL;
		}
		import_array(); 
		return module;
	}
#else

	void init_uniform_ld(void)
	{
	  	Py_InitModule("_uniform_ld", _uniform_ld_methods);
		import_array(); 
	}
#endif

