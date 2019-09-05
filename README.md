
**T**ransit and **R**adial velocity **I**nteractive **F**itting tool for **O**rbital analysis and **N**-body simulations : **The Exo-Striker** 

<p align="center">
  <img width="400" src="https://github.com/3fon3fonov/trifon/blob/master/lib/33_striker.png">
</p>
 
The Exo-Striker analyzes exoplanet orbitals, performs N-body simulations, and models the RV stellar reflex motion caused by dynamically interacting planets in multi-planetary systems. It offers a broad range of tools for detailed analysis of transit and Doppler data, including power spectrum analysis for Doppler and transit data; Keplerian and dynamical modeling of multi-planet systems; MCMC and nested sampling; Gaussian Processes modeling; and a long-term stability check of multi-planet systems. The Exo-Striker can also perform Mean Motion Resonance (MMR) analysis, create fast fully interactive plots, and export ready-to-use LaTeX tables with best-fit parameters, errors, and statistics. It combines Fortran efficiency and Python flexibility and is cross-platform compatible (MAC OS, Linux, Windows). The tool relies on a number of open-source packages, including RVmod engine (Trifonov et al. in prep.), emcee (ascl:1303.002), batman (ascl:1510.002), celerite (ascl:1709.008), and dynesty (ascl:1809.013).

![new_es](https://user-images.githubusercontent.com/44244057/62121145-48b99700-b2c3-11e9-9411-7fc3d1a4467d.png)


**WARNING!** This tool is under active development and its functionality is enhanced on a daily basis! Therefore, although very unlikely, the version you download today may not be fully compatible with the version uploaded tomorrow! Use at your own risk!

Also, please keep in mind that this software is developed mostly for my needs and for fun. I hope, however, that you may find 
it capable to solve your scientific problems, too. If you made the use of The Exo-Striker for your paper please 
cite the tool with its ASCL ID ascl:1906.004 (see https://ascl.net/1906.004).


What works:

* Period search: GLS periodograms (RVs, act. data) & TLS (transit data).
* Keplerian and Dynamical RV modeling. 
* Transit photometry modeling.
* GP modeling (only SHO and Rot. GP kernels integrated so far).
* Joint RVs + Transit + GPs best-fit optimization.
* Joint RVs + Transit + GPs MCMC/Nested Sampling (via "emcee" & "dynesty").
* RV auto-fit (RV automated planet-finder algorithm).
* Possible to fit for Precession and/or GR precession. 
* Long-term stability check of multi-planet systems using SyMBA, MVS and MVS with an GR precession.
* Various of minimization methods (via "SciPyOp").
* Fully interactive, super-fast, high-quality, exportable plots.
* RV vs. Activity time series correlation analysis/plots.
* Import/Export of work sessions and multi-sessions. 
* Export plots to a matplotlib window for further customization.
* Export ready to use LaTeX tables with best-fit parameters, errors, and statistics. 
* Handy text editor and calculator tools.
* Multi-platform: It works on MAC OS (10.6+), Linux (Suse, Mint, Ubuntu, etc.) and Windows 10.
* Integrated Bash-shell (Linux only).
* Integrated Jupyter shell.
* Importable as a standard python library (i.e. ideal for scripting and notebooks, see "Examples").
* Print the GUI screen into a .jpeg image (useful for sharing quick results, just like the image above)

What is to be implemented:

* More GP kernels (work in progress). 
* Larger arsenal of N-body/dynamical simulation/analysis tools (REBOUND to be included). 
* A pip installer, and a ready-to-use pre-installed binary of the tool (work in progress). 
* Binary/Triple star modeling mode.
* Combined modeling with Astrometry.
* Documentation, Instructions and Video tutorials.
* For more "TBD" list see: "focus_matrix_TBFixed.doc".

If you use this tool and you find a bug or a problem, please report it!
A wish-list with your favorite tools and methods to be implemented is also welcome!
Just open an "Issue" on the github, or send a PM to trifonov@mpia.de.



