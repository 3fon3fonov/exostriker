.. _tutorials2:

Fiting the RV data
..................

If you want to use the library on the Python shell/script

In [1]: import exostriker

or e.g. to load the RV routines:

In [1]: import exostriker.lib.RV_mod as rv   

In [2]: fit = rv.signal_fit(name="hip5364") #creates the "fit" object that contains everything.    

In [3]: fit.add_dataset("hip5364-Lick","./datafiles/hip5364.vels",0.0.10.0) # add the data file, initial offset and jitter   

In [4]: fit.add_planet(K=50,P=400,e=0,w=0,M0=0,i=90,cap=0)   # planet 1    

In [5]: fit.add_planet(K=50,P=700,e=0,w=0,M0=180,i=90,cap=0) # planet 2    

In [6]: fit.fitting() #optimize the parameters    

In [7]: fit.run_mcmc() # run MCMC, etc...    
 
(However, one must be familiar with the functions... A manual on RVmod is planned, but not available at the moment.)



* All Fortran and python codes in this version need serious clean-up from junk.

* All Fortran codes are planned to be translated with f2py into python-importable libraries.

* Don't ever run MCMC runs (or Nest. Samp.) directly in the embedded Jupyter shell! This will freeze the GUI until is the MCMC done!
  This is not a bug, simply the jupyter shell needs its thread and this is not done, yet. Hopefully, this will be fixed soon.
  Please use the GUI navigation.
