
Welcome to the complete beginner's guide to Exostriker!
=======================================================

Review
......

Transit and Radial velocity Interactive Fitting tool for Orbital analysis and N-body simulations: 
The Exo-Striker.

The Exo-Striker analyzes exoplanet orbitals, performs N-body simulations, and models the RV stellar reflex motion
caused by dynamically interacting planets in multi-planetary systems. It offers a broad range of tools for detailed 
analysis of transit and Doppler data, including power spectrum analysis for Doppler and transit data, Keplerian and 
dynamical modeling of multi-planet systems, MCMC and nested sampling, Gaussian Processes modeling, and a long-term 
stability check of multi-planet systems. The Exo-Striker can also perform Mean Motion Resonance (MMR) analysis, create
fast fully interactive plots, and export ready-to-use LaTeX tables with best-fit parameters, errors, and statistics. 
It combines Fortran efficiency and Python flexibility and is cross-platform compatible (MAC OS, Linux, Windows).

Features
........

* Keplerian and Dynamical modeling of RV & Transit photometry exoplanet data.
* RV signal and alias search: via GLS periodogram & maximum lnL periodogram (MLP).
* Activity index signal search via GLS periodogram.
* RVs vs. Activity time-series correlation analysis/plots.
* RV auto-fit (RV automated planet-finder algorithm).
* Transit signal search (via "TLS").
* Interactive transit photometry detrending (via "wotan"), interactive outlier removal, and more.
* Linear models for detrending ground-based transit photometry.
* GP modeling (via "celerite").
* Joint RVs + Transit + GPs best-fit optimization (internal Fortran Simplex and L-M minimizers, or many more via "SciPyOp").
* Joint RVs + Transit + GPs MCMC/Nested Sampling (via "emcee" & "dynesty").
* TTVs extraction.
* TTVs and/or joint TTVs + RVs analysis.
* Fit for apsidal orbital precession, or apply General Relativity (GR) precession.
* Long-term stability check of multiple planetary systems using SyMBA, MVS, and MVS with a GR precession.
* Instant AMD stability check for multiple planetary systems (including during optimization or MCMC/Nested Sampling).
* Multi-platform: It works on MAC OS (10.6+), Linux (Suse, Mint, Ubuntu, etc.) and Windows 10.
* Import/Export of work sessions and multi-sessions.
* Export plots to a matplotlib window for further customization.
* Export ready to use LaTeX tables with best-fit parameters, errors, and statistics.
* Print the GUI screen into a .jpeg/.png image (useful for sharing quick results).
* Fully interactive, super-fast, high-quality, exportable plots.
* Integrated Bash-shell (Linux only).
* Integrated Jupyter shell.
* Handy "cornerplot" GUI control.
* Instant online access to the "RVBank" database (over 212 000 RVs and activity indices of about 3000 HARPS stars & over 64 000 RVs and activity indices of about 1700 HIRES stars).
* Handy text-editor and calculator tools.
* Direct import of TESS lc.fits and CHEOPS SCI_COR.fits files.
* Importable as a standard python library (i.e., "import exostriker").

What is to be implemented
.........................

* Larger arsenal of N-body/dynamical simulation/analysis tools (+ "REBOUND" is planned to be included).
* Internal TTV and photo-dynamical modeling (i.e. the external "TTVFast" will become a secondary option).
* Swap "celerite" with "celerite2" (the dSHO kernel from "celerite2" is available).
* Combined modeling with Astrometry (work in progress).

Developers
..........

* **Trifon Trifonov**, MPIA Heidelberg.
* with contributions by **Mathias Zechmeister**, **Jakub Morawski**, **Man Hoi Lee**, **Stefan Dreizler**, **Grigorii Smirnov-Pinchukov**, **Stephan Stock**, **Jonas Kemmer** and **Harry Psarakis**.

.. toctree::
   :hidden:
   :maxdepth: 5
   :caption: User guide

   userguide

.. toctree::
   :hidden:
   :maxdepth: 6
   :caption: tutorials based on gui
   
   gui
   rvs
   transit
   rvtran
   otbfpu
   stability

.. toctree::
   :hidden:
   :maxdepth: 5
   :caption: tutorials based on Jupyter shell

   tutorials2

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Credits

   credits

