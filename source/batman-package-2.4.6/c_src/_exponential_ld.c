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
#include "common.h"


static PyObject *_exponential_ld(PyObject *self, PyObject *args);

inline double intensity(double x, double* args)
{
	double c1=args[0], c2=args[1], norm=args[2];
	if(x > 0.99995) x = 0.99995;
	if(x < 0.00005) x = 0.00005;
	double mu = sqrt(1. - x*x);
	return (1. - c1*(1. - mu) - c2/(1. - exp(mu)))/norm; 
}


static PyObject *_exponential_ld(PyObject *self, PyObject *args)
{
	double rprs, fac;
	int nthreads;
	npy_intp dims[1];
	double c1, c2;
	
	PyArrayObject *ds, *flux;
  	if(!PyArg_ParseTuple(args,"Oddddi", &ds, &rprs, &c1, &c2, &fac, &nthreads)) return NULL;

	dims[0] = PyArray_DIMS(ds)[0]; 
	flux = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(ds));	//creates numpy array to store return flux values

	double *f_array = PyArray_DATA(flux);
	double *d_array = PyArray_DATA(ds);

	/*
		NOTE:  the safest way to access numpy arrays is to use the PyArray_GETITEM and PyArray_SETITEM functions.
		Here we use a trick for faster access and more convenient access, where we set a pointer to the 
		beginning of the array with the PyArray_DATA (e.g., f_array) and access elements with e.g., f_array[i].
		Success of this operation depends on the numpy array storing data in blocks equal in size to a C double.
		If you run into trouble along these lines, I recommend changing the array access to something like:
			d = PyFloat_AsDouble(PyArray_GETITEM(ds, PyArray_GetPtr(ds, &i))); 
		where ds is a numpy array object.


		Laura Kreidberg 07/2015
	*/
	
	double norm = 2.*M_PI*(0.5 - 0.1666666667*c1 + 0.77750463*c2); 	//normalization for intensity profile (faster to calculate it once, rather than every time intensity is called)		
	double intensity_args[] = {c1, c2, norm};
	#pragma acc data copyin(intensity_args)
	calc_limb_darkening(f_array, d_array, dims[0], rprs, fac, nthreads, intensity_args);
	
	return PyArray_Return((PyArrayObject *)flux);

} 

static char _exponential_ld_doc[] = "This extension module returns a limb darkened light curve for a exponential stellar intensity profile.";

static PyMethodDef _exponential_ld_methods[] = {
  {"_exponential_ld", _exponential_ld, METH_VARARGS, _exponential_ld_doc},{NULL}};

#if PY_MAJOR_VERSION >= 3
	static struct PyModuleDef _exponential_ld_module = {
		PyModuleDef_HEAD_INIT,
		"_exponential_ld",
		_exponential_ld_doc,
		-1, 
		_exponential_ld_methods
	};

	PyMODINIT_FUNC
	PyInit__exponential_ld(void)
	{
		PyObject* module = PyModule_Create(&_exponential_ld_module);
		if(!module)
		{
			return NULL;
		}
		import_array(); 
		return module;
	}
#else

	void init_exponential_ld(void)
	{
	  	Py_InitModule("_exponential_ld", _exponential_ld_methods);
		import_array(); 
	}
#endif

