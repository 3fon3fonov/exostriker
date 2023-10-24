
# Transit and Radial velocity Interactive Fitting tool for Orbital analysis and N-body simulations: The Exo-Striker

<p align="center">
  <img width="400" src="./exostriker/source/png/33_striker.png">
</p>
 
The Exo-Striker analyzes exoplanet orbitals, performs N-body simulations, and models the RV stellar reflex motion caused by dynamically interacting planets in multi-planetary systems. It offers a broad range of tools for detailed analysis of transit and Doppler data, including power spectrum analysis for Doppler and transit data; Keplerian and dynamical modeling of multi-planet systems; MCMC and nested sampling; Gaussian Processes modeling; and a long-term stability check of multi-planet systems. The Exo-Striker can also perform Mean Motion Resonance (MMR) analysis, create fast fully interactive plots, and export ready-to-use LaTeX tables with best-fit parameters, errors, and statistics. It combines Fortran efficiency and Python flexibility and is cross-platform compatible (MAC OS, Linux, Windows).

 <p align="center">
  <img width="800" src="./exostriker/source/png/ES_new.gif">
</p>


### Documentation, Instructions, and Tutorials

* Documentation is available at [https://exostriker-manual.readthedocs.io](https://exostriker-manual.readthedocs.io/en/latest/). (Work in progress)

### Developer

* Trifon Trifonov, MPIA Heidelberg.
* With contributions by: Mathias Zechmeister, Jakub Morawski, Man Hoi Lee, Stefan Dreizler, Grigorii Smirnov-Pinchukov, Stephan Stock, Jonas Kemmer, Harry Psarakis, Tom Schiwy, and Desislava Antonova.

### Features

- [x] RV signal and alias search: via GLS periodogram & maximum lnL periodogram (MLP).
- [x] Transit signal search (via "TLS").
- [x] Interactive transit photometry detrending (via "wotan"), interactive outlier removal, and more.
- [x] Keplerian and Dynamical modeling of RV & Transit photometry exoplanet data.
- [x] Joint RVs + Transit + GPs best-fit optimization (internal Fortran Simplex and L-M minimizers, or many more via "SciPyOp").
- [x] Joint RVs + Transit + GPs MCMC/Nested Sampling (via "emcee" & "dynesty") 
- [x] TTVs extraction.
- [x] TTVs and/or joint TTVs + RVs analysis.
- [x] Relative astrometry fitting.
- [x] GP modeling (via "celerite").
- [x] Linear models for detrending ground-based transit photometry.
- [x] Activity index signal search via GLS periodogram.
- [x] RVs vs. Activity time-series correlation analysis/plots.
- [x] RV auto-fit (RV automated planet-finder algorithm).
- [x] Fit for apsidal orbital precession, or apply General Relativity (GR) precession. 
- [x] Instant online access to the "RVBank" database (over 212 000 RVs and activity indices of about 3000 HARPS stars & over 64 000 RVs and activity indices of about 1700 HIRES stars !!!).
- [x] Instant AMD stability check for multiple planetary systems (including during optimization or MCMC/Nested Sampling).
- [x] Long-term stability check of multiple planetary systems using SyMBA, MVS, and MVS with a GR precession.
- [x] Fully interactive, super-fast, high-quality, exportable plots.
- [x] Handy "cornerplot" GUI control.
- [x] Import/Export of work sessions and multi-sessions. 
- [x] Export plots to a matplotlib window for further customization.
- [x] Export ready-to-use LaTeX tables with best-fit parameters, errors, and statistics. 
- [x] Handy text-editor and calculator tools.
- [x] Multi-platform: It works on MAC OS (10.6+), Linux (Suse, Mint, Ubuntu, etc.) and Windows 10.
- [x] Integrated Bash-shell (Linux only).
- [x] Integrated Jupyter shell.
- [x] Integrated AI Boot (via openai GPT-3 -- ChatGPT).
- [x] Importable as a standard python library (i.e., "import exostriker").
- [x] Print the GUI screen into a .jpeg/.png image (useful for sharing quick results, just like the image above).
- [x] Direct import of TESS & K2 *lc.fits, and CHEOPS *SCI_COR*.fits files.
- [ ] Larger arsenal of N-body/dynamical simulation/analysis tools (+ "REBOUND" is planned to be included). 
- [ ] Internal TTV and photo-dynamical modeling (i.e. the external "TTVFast" will become a secondary option).
- [ ] Combined modeling with Astrometry (As of Ver 0.75 this is possible, but is still work in progress).


Please keep in mind that this software is developed mostly for my needs and fun. I hope, however, that you may find it capable of solving your scientific problems too. At the moment, there is NO documentation,
but as you will find, the GUI is self-explanatory. The import as a Python library indeed needs documentation, but I plan to write such after approx. Ver. 0.8, when
most of the functions and classes will be firmly defined and will start optimizing the software core. 
If you need the Exo-Striker as a library, you can check the [Notebook_and_script_examples](https://github.com/3fon3fonov/exostriker/tree/exostriker-ready/Notebook_and_script_examples),
or you can drop me a line at trifonov@mpia.de for further instructions. 

Feedback and help in further development will be highly appreciated!
A wish-list with your favorite tools and methods to be implemented is also welcome!    

Just open an "Issue" on the GitHub, or send a PM to trifonov@mpia.de.    


### Installation
 
Python 3.8+ is strongly recommended!!! The Exo-Striker works with Python 3.6, 3.7 and even Python 2, but you will likely have problems with some dependencies, which you may have to solve.
If your Python version is 3.6 or 3.7, it is recommended to not upgrade the `python3`, but to install `python3.8` alongside your system `python3` (see the instructions for [Python 3.8](https://linuxize.com/post/how-to-install-python-3-8-on-ubuntu-18-04/), or for [Python 3.9](https://linuxize.com/post/how-to-install-python-3-9-on-ubuntu-20-04/)). 

Python 3.10 is now supported! Please note that the ES with Python 3.8 and Python 3.9 will not be supported starting from July 1st 2023. 

For Ubuntu 20.04+ you may also need `python-dev` to install some of the dependencies. For example, if during installation you get something like
```sh
$ 20 | #include "Python.h"   
$    |          ^~~~~~~~~~   
$  compilation terminated.
```
, you may have to do
```sh
$ sudo apt-get install python3.8-dev
```
and try again.

Apart from Python 3.8/3.9, make sure you have `gcc`, `gfortran`, and `csh`, which will be needed for compiling some of the important binaries!

#### Installation through pip

Assuming you have Python 3.8, then it is recommended to start by installing `pip3.8`:
```sh
$ curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py   
$ python3.8 get-pip.py    
```
Then, to install the Exo-Striker:
```sh
$ pip3.8 install git+https://github.com/3fon3fonov/exostriker --user
# (for user installation -- recommended!)
```
Then just run
```sh
$ exostriker
```
**WARNING!** if you install with
```sh
$ sudo pip3.8 install git+https://github.com/3fon3fonov/exostriker
```
, then every time you should run
```sh
$ sudo exostriker 
```

#### Installation through git

Alternatively, you can clone the repository through `git`:
```sh
$ git clone https://github.com/3fon3fonov/exostriker   
$ cd exostriker   
$ python3.8 setup.py install   
```
However, please read the [Installation instructions](README_for_installation),
because some problems may occur depending on your OS system. Finally, you can run the Exo-striker by
```sh
$ python3.8 exostriker_gui.py  
```
The last command requires that all the dependencies, specified in `setup.py`, are installed. This is somewhat more flexible IMHO, since one can have multiple "exostriker" directories dedicated to different projects. I usually go with: "exostriker-GJ436", "exostriker-GJ876" etc.; and I store all datafiles, sessions, plots, drafts etc.
inside these directories.
 
### Usage

To load the GUI, type the following command in bash:
```sh
$ exostriker
# (in case of pip3.8 install, see above)
```
or just do
```sh
$ python3.8 exostriker_gui.py
# (inside of the git clone directory, see above)
``` 
If you want to use the library on the Python shell/script:
```py
import exostriker
```
or, for example, to load the RV routines:
```py
import exostriker.lib.RV_mod as rv    
``` 
Remember! Every button/action of the GUI is a wrapper of a convenient Python routine. This makes scripting very easy:    

<p align="center">
  <img width="400" src="./exostriker/source/png/ES_terminal.png">
</p>

 
(However, one must be familiar with the functions and the 'fit' object... 
A manual is planned but not available at the moment.)


### Credit

If you made use of The Exo-Striker for your paper, I would appreciate it if you give credit to it.
As it is unlikely that I will find time to write a refereed paper on The Exo-Striker soon, please cite the tool with its ASCL ID ascl:1906.004 (see https://ascl.net/1906.004).    
 
The Exo-Striker relies on many open-source packages. If you had made use of some of them while working with the tool, you should acknowledge them too (it is your responsibility to find the correct references in the literature):    


* The interactive plotting is done with a custom version of pyqtgraph: http://www.pyqtgraph.org/

* "GLS" and "MLP" periodograms are taken from Mathias Zechmeister's repo: https://github.com/mzechmeister/python

* "TLS" and "wotan" are taken from: https://github.com/hippke/tls , https://github.com/hippke/wotan

* The transit modeling is done with batman: https://github.com/lkreidberg/batman

* MCMC sampling is done with emcee: https://github.com/dfm/emcee

* Nested Sampling is done with dynesty: https://github.com/joshspeagle/dynesty

* TTV models are adopted from TTVfast-python: https://github.com/mindriot101/ttvfast-python

* The "Text editor" used in the tool is a hack between "Megasolid Idiom" and "PyEdit2": https://github.com/mfitzp/15-minute-apps/tree/master/wordprocessor , https://github.com/Axel-Erfurt/PyEdit2

* N-body tests are performed using a custom version of the Swift N-body library,
modified by Man Hoi Lee (HKU) and Trifon Trifonov (MPIA): https://www.boulder.swri.edu/~hal/swift.html

* Additionally, the Exo-Striker uses many "standard" Python libraries like:
PyQt5,
matplotlib,
numpy,
scipy,
dill,
Jupyter,
qtconsole
and more.

* The Exo-Striker project was inspired by the Systemic project: http://www.stefanom.org/systemic/

### Scientific papers which one way or another made the use of the Exo-Striker (to my knowledge):
 
[Check in ADS.](https://ui.adsabs.harvard.edu/abs/2019ascl.soft06004T/citations)
