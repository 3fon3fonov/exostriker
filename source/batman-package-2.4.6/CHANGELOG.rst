.. :changelog:
2.4.6 (2017-11-25)
~~~~~~~~~~~~~~~~~~
- ensure numerical stability for pathological transit depth case

2.4.5 (2017-05-11)
~~~~~~~~~~~~~~~~~~
- correct tab spacing inconsistency

2.4.4 (2017-05-08)
~~~~~~~~~~~~~~~~~~
- ok now actually added .h files to MANIFEST.in :)

2.4.3 (2017-05-08)
~~~~~~~~~~~~~~~~~~
- added .h files to MANIFEST.in 

2.4.2 (2017-05-08)
~~~~~~~~~~~~~~~~~~
- added -std=c99 compile flags

2.4.1 (2017-05-07)
~~~~~~~~~~~~~~~~~~
- fix calculation of eccentric anomaly to handle diabolical inputs (following Eastmane et al. 2013)
- Optimized all files for GPU
- Restructured code so that _power2, _exponential, _logarithmic, _nonlinear, and _custom_ld all use a common code path
- Fixed bug where openmp.py duplicated a file handle without closing, causing operating system to run out of file handles after repeated re-runs

2.4.0 (2017-05-03)
~~~~~~~~~~~~~~~~~~
- add power2 limb darkening model from Morello et al. 2017
- bug fix for secondary eclipse calculation (make variables private to avoid race condition)
- OpenACC to take advantage of GPU acceleration (thanks to Michael Zhang)

2.3.0 (2015-03-11)
~~~~~~~~~~~~~~~~~~
- add get_true_anomaly() method
- remove redundant arrays
- improve accuracy for special cases (e.g. rp near 0.5)

2.2.0 (2015-12-19)
~~~~~~~~~~~~~~~~~~
- add inverse transit capability (can now handle negative rp)
- speed up super-sampling


2.1.0 (2015-08-06)
~~~~~~~~~~~~~~~~~~
- add get_t_conjunction() method 
- change eclipse model normalization so that stellar flux is unity

2.0.0 (2015-08-04)
~~~~~~~~~~~~~~~~~~
- add secondary eclipse model
- change model parameterization to time of inferior conjunction from time of periastron (backwards-incompatible change in response to referee)


1.0.0 (2015-07-29)
~~~~~~~~~~~~~~~~~~
- first stable release


0.9.1 (2015-06-24)
~~~~~~~~~~~~~~~~~~

- fixing bug in call to _rsky


0.9.0 (2015-06-24)
~~~~~~~~~~~~~~~~~~

- Beta version 
