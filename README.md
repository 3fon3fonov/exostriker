
**T**ransit and **R**adial velocity **I**nteractive **F**itting tool for **O**rbital analysis and **N**-body simulations: **The Exo-Striker** 

<p align="center">
  <img width="400" src="https://github.com/3fon3fonov/trifon/blob/master/lib/UI/33_striker.png">
</p>
 
The Exo-Striker analyzes exoplanet orbitals, performs N-body simulations, and models the RV stellar reflex motion caused by dynamically interacting planets in multi-planetary systems. It offers a broad range of tools for detailed analysis of transit and Doppler data, including power spectrum analysis for Doppler and transit data; Keplerian and dynamical modeling of multi-planet systems; MCMC and nested sampling; Gaussian Processes modeling; and a long-term stability check of multi-planet systems. The Exo-Striker can also perform Mean Motion Resonance (MMR) analysis, create fast fully interactive plots, and export ready-to-use LaTeX tables with best-fit parameters, errors, and statistics. It combines Fortran efficiency and Python flexibility and is cross-platform compatible (MAC OS, Linux, Windows). 

![new_es](/lib/png/Exo_striker_demo_image.png)
 

What works:

* RV signal and alias search: via GLS periodogram & maximum lnL periodogram (MLP).
* Activity index signal search via GLS periodogram.
* Keplerian and Dynamical RV modeling. 
* Instant online access to the HARPS RVBank (over 212 000 RVs and activity indices of about 3000 stars!!!).
* Transit signal search via TLS.
* Transit photometry detrending (via "wotan").
* Transit photometry modeling.
* GP modeling (only SHO and Rot. GP "celerite" kernels integrated so far).
* Joint RVs + Transit + GPs best-fit optimization.
* Joint RVs + Transit + GPs MCMC/Nested Sampling (via "emcee" & "dynesty").
* TTVs and joint TTVs + RVs analysis.
* RV auto-fit (RV automated planet-finder algorithm).
* Fit for apsidal orbital precession, or apply General Relativity (GR) precession. 
* Instant AMD stability check for multiple planetary systems (including during optimization or MCMC/Nested Sampling).
* Long-term stability check of multiple planetary systems using SyMBA, MVS, and MVS with a GR precession.
* Various of minimization methods (via "SciPyOp").
* Fully interactive, super-fast, high-quality, exportable plots.
* RV vs. Activity time-series correlation analysis/plots.
* Import/Export of work sessions and multi-sessions. 
* Export plots to a matplotlib window for further customization.
* Export ready to use LaTeX tables with best-fit parameters, errors, and statistics. 
* Handy text-editor and calculator tools.
* Multi-platform: It works on MAC OS (10.6+), Linux (Suse, Mint, Ubuntu, etc.) and Windows 10.
* Integrated Bash-shell (Linux only).
* Integrated Jupyter shell.
* Importable as a standard python library (i.e., ideal for scripting and notebooks, see "Notebook_and_script_examples").
* Print the GUI screen into a .jpeg/.png image (useful for sharing quick results, just like the image above).

What is to be implemented:

* More GP kernels (work in progress). 
* Larger arsenal of N-body/dynamical simulation/analysis tools (+ REBOUND is planned to be included). 
* A photo-dynamical transit model.
* Internal TTV modeling (i.e. the external "TTVFast" will become a secondary option).
* A pip installer, and a ready-to-use pre-installed binary of the tool (work in progress). 
* Combined modeling with Astrometry.
* Documentation, Instructions, and Video tutorials (work in progress here: https://exostriker.readthedocs.io)


**Installation:**

Since the Exo-Striker tool is under active development, and its functionality is enhanced on a daily basis,
I intentionally did not make a "pip" installer. This will be done close to the official Ver.1 release. The only way to install the tool is simply to "git clone" and then run "pyhton es_gui.py" inside the root directory.
Generally, you do not need to install anything if you already have all the dependencies.
The Exo-Striker will automatically compile the Fortran code for you at the first start of the program and will keep you updated if the source code was updated. 
Yet, to avoid problems, I recommend to follow the steps listed below:


**for a first download, just do:**

$ git clone https://github.com/3fon3fonov/exostriker  

**then, e.g. for an Ubuntu (Debian) installation:**

$ cd exostriker  
$ bash installers/linux_debian_install.sh  

**Successful installation should look like this:**

For which Python version we want to check/install packages?  
1) Python2  
2) Python3  
*#? 2*  

**(now it will check if you have the needed dependencies. If packages are missing it will ask you to install it")**

gfortran - yes!  
csh - yes!  
setuptools - yes!  
pip - yes!  
numpy - yes!  
scipy - yes!  
matplotlib - yes!  
PyQt5 - yes!  
PyQt5.QtSvg - yes!   
qtconsole - yes!   
jupyter - yes!   
pathos - yes!  
dill - yes!  
emcee - yes!  
corner - yes!  
celerite - yes!  
transitleastsquares - yes!  
dynesty - yes!  
rxvt - yes!  
batman - yes!   
ttvfast - yes!   
wotan - yes!   

*Installing the swift N-body lib, OK?  (you must if you haven't done it already!)*
  
1) Yes  
2) No   
*#? 1*  

**(some output will be printed, it takes ~1 min.)**  
DONE
 

*Compiling the Fortran fitting routines, OK? (you must if you haven't done it already!)*
 
1) Yes  
2) No   
*#? 1*  
 
*Compiling Symba/mvs and other N-body routines, OK? (you must if you haven't done it already!)*
 
1) Yes  
2) No   
*#? 1*   

**That's it! Then simply:**

$ python es_gui.py

**to update the tool (inside the "exostriker" directory) just do:**

$ git pull  



**Also please read "README_for_installation"!**

I believe the installation instructions are very clear and easy to run.

If you still cannot boot the tool after a 'successful' installation, please try:

$ python es_gui.py -debug 

Then, copy the output error and please open a 'GitHub' issue. Otherwise, all possible problems/bugs/crashes will be displayed on the 
'stdout/stderr' tab of the tool. If you use this tool and you find a bug or a problem, please report it!

A wish-list with your favorite tools and methods to be implemented is also welcome!
Just open an "Issue" on the GitHub, or send a PM to trifonov@mpia.de.



**Credit**

If you made the use of The Exo-Striker for your paper, I would appreciate if you give credit to it.
As it is unlikely that I will find time to write a refereed paper on the Exo-Striker soon, please cite the tool with its ASCL ID ascl:1906.004 (see https://ascl.net/1906.004).
 
The Exo-Striker relies on many open-source packages, which if you had made the use of them while working with the tool, 
you should acknowledge too. (It is your responsibility to find the correct references in the literature):    


* The interactive plotting is done with a custom version of the "pyqtgraph": 

http://www.pyqtgraph.org/

* "GLS" and "MLP" periodograms are taken from Mathias Zechmeister's repo: 

https://github.com/mzechmeister/python

* "TLS" and "wotan" are taken from: 

https://github.com/hippke/tls

https://github.com/hippke/wotan

* The transit modeling is done with "batman":
 
https://github.com/lkreidberg/batman

* MCMC sampling is done with "emcee": 

https://github.com/dfm/emcee

* Nested Sampling is done with "dynesty": 

https://github.com/joshspeagle/dynesty

* TTV models are adopted from "TTVfast-python":

https://github.com/mindriot101/ttvfast-python

* The "Text editor" used in the tool is a hack between "Megasolid Idiom" 
and "PyEdit2":

https://github.com/mfitzp/15-minute-apps/tree/master/wordprocessor

https://github.com/Axel-Erfurt/PyEdit2

* N-body tests are performed using a custom version of the "Swift" N-body library,
modified by Man Hoi Lee (HKU) and Trifon Trifonov (MPIA).

https://www.boulder.swri.edu/~hal/swift.html

* Additionally, the Exo-Striker uses many "standard" Python libraries like 
"PyQt5", "matplotlib", "numpy", "scipy", "dill", "Jupyter", "qtconsole",
and more.


**Papers used The Exo-Striker (to my knowledge)**

* Stock, S., J. Kemmer, S. Reffert, et al. (2020). The CARMENES search for exoplanets around M dwarfs. Characterization of the nearby ultra-compact multiplanetary system YZ Ceti. (A&A, in press) arXiv e-prints arXiv:2002.01772.
 
* Trifonov, T., M. H. Lee, M. Kürster, et al. (2020). The CARMENES search for exoplanets around M dwarfs. Dynamical characterization of the multiple planet system GJ 1148 and prospects of habitable exomoons around GJ 1148 b. (A&A, in press) arXiv e-prints arXiv:2002.00906.

* Trifon Trifonov, Lev Tal-Or, Mathias Zechmeister, et al. (2020). A public HARPS radial velocity database corrected for systematic errors. A&A, 636, A74 
 
* Luque, R., T. Trifonov, S. Reffert, et al. (2019). Precise radial velocities of giant stars. XIII. A second Jupiter orbiting in 4:3 resonance in the 7 CMa system. ApJ, 631, A136.

* Trifonov, Trifon, Stephan Stock, Thomas Henning, et al. (2019). Two Jovian Planets around the Giant Star HD 202696: A Growing Population of Packed Massive Planetary Pairs around Massive Stars?. AJ, 157, 93.

* Trifonov, Trifon, Jan Rybizki, and Martin Kürster. (2019). TESS exoplanet candidates validated with HARPS archival data. A massive Neptune around GJ 143 and two Neptunes around HD 23472. A&A, 622, L7. 

(And as of April 2020, at least 7 more papers in preperation !!!)


