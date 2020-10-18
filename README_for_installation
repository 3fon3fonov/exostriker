To install the Exo-Striker, first please make sure that you have:

* csh
* gfortran
* pip


The instructions below will install the following dependencies: 

* numpy
* scipy
* matplotlib
* PyQt5
* jupyter
* pathos
* emcee  
* celerite
* qtconsole
* dynesty
* wotan 
* corner
* transitleastsquares
* ttvfast-python


Also, please only use Python 3 for installing the Exo-Striker! The tool works with Python 2, 
but since Python 2 is no longer maintained (since Jan 1, 2020), I do not assist in case of problems.    

(see below comments on Linux, Mac OS, and Windows 10 installations)




Currently there are three ways to install/run the Exo-Striker:    

#######################################################    
*  The simplest way to "git clone":    

$ git clone https://github.com/3fon3fonov/exostriker  

and then:    

$ cd exostriker  
$ python3 exostriker_gui.py    

Generally, you do not need to install anything if you already have all the required dependencies for the tool to run. For the dependency list, see the "setup.py" file. The Exo-Striker will automatically compile the Fortran77 code for you at the first start of the program and will keep you updated if the source code was updated (if you regularly "git pull").    
 
#######################################################    
*  The second way to install the tool from the source:    

$ git clone https://github.com/3fon3fonov/exostriker     

and then:    

$ cd exostriker    
$ python3 setup.py install    

This will install the tool in your system.     
Then, open a terminal and:     

$ exostriker    

Should start the Exo-Striker  

#######################################################    
*  and last, you can try pip install:    

$ pip install git+https://github.com/3fon3fonov/exostriker    

This will install the tool in your system.    
Then, open a terminal and:     

$ exostriker

Should start the Exo-Striker  

#######################################################     

If you still cannot boot the tool after a 'successful' installation, please try:

$ python3 exostriker_gui.py -debug 

or 

$ exostriker -debug 

(depending on how you use the tool)

Then, copy the output error, and please open a 'GitHub' issue. Otherwise, all possible problems/bugs/crashes will be displayed on the 
'stdout/stderr' tab of the tool. If you use this tool, and you find a bug or a problem, please report it!    

#######################################################




Some OS related comments:


##############################  LINUX  #########################################

The above instructions usually work without any problem on Linux OS.
 
 
For full functionality on Linux, you will also need to install:

* rxvt (optional, better bash shell)
 
 

################################  MAC OS  ######################################
 
You will need to install "homebrew", which requires ``sudo'' privileges. 
According to 

https://brew.sh/

To install "homebrew", it should be as simple as this:


/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

try it out in bash shell terminal.


If you have MAC OS 10.14+ it is likely to have problems with the numpy headers.

if you see an error like this during the installations (mostly during the "batman-package" build)

...
...
    c_src/_nonlinear_ld.c:21:10: fatal error: 'numpy/arrayobject.h' file not found
    #include "numpy/arrayobject.h"
             ^~~~~~~~~~~~~~~~~~~~~
    1 error generated.
    error: command 'clang' failed with exit status 1
...
...
Please do the following:

export CFLAGS="-I /usr/local/lib/pythonX.X/site-packages/numpy/core/include $CFLAGS"

(Where X is your Python version, e.g., 3.6, 3.7 to 3.8)

and then re-run your Exo-Striker installation (via pip or setup.py).
 

################################  WINDOWS 10  ######################################


Installation on Windows 10 works troughs the "Windows Subsystem for Linux".
Please follow this guide:

https://docs.microsoft.com/en-us/windows/wsl/install-win10

This way you will be able to run all Linux native programs on your WINDOWS 10 
bash shell, which is very useful in general!

To make The Exo-Striker work, however, you also will need an XServer installed.
Follow these instructions:

https://seanthegeek.net/234/graphical-linux-applications-bash-ubuntu-windows/


These two tutorials worked for me, but there might be other options too. 


In case there is a problem of the appearance of the GUI the problem could be 
your DPI setup. On high DPI displays with Windows display scaling activated (>100%),
the X-server DPI has to be set, otherwise, the Exo-Striker will not display correctly (text clipping).
Edit the file `.Xresources` in the home directory of the WSL installation,
for example with `nano ~/.Xresources`.
Add the line `Xft.dpi: 100` and save & close with Ctrl+O Ctrl+X,
then reload the X configuration with `xrdb ~/.Xresources`.
On Ubuntu WSL you might need to install `x11-xserver-utils`
with `sudo apt install x11-xserver-utils` for xrdb to be available.

Launch the Exo-Striker and check the scaling.
If text is clipping, the DPI needs to be set lower, if everything is too small,
the dpi needs to be higher. Remember always to reload the configuration with `xrdb ~/.Xresources`.
For the configuration to automatically load at startup,
add the xrdb command to your ~/.bashrc, after the `export DISPLAY=:0.0`.


Running the tool via the official Windows 10 python3 installation should generally work too,
it was never tried! If you want to experiment, and you successfully install the tool under the official
Windows python path, I would appreciate it if you share your experience and some instructions.


For now, the recommended WINDOWS 10 installation option of the Exo-Striker is via the "Windows 
Subsystem for Linux" as pointed above.
 
 
 
################################################################################
########################### Some known problems   ##############################
################################################################################
 
 
 See issues:
 
 https://github.com/3fon3fonov/exostriker/issues
 
################################## Finally #####################################

To load the GUI 

$ python3 exostriker_gui.py 
 
 
 
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




Some comments:

* All Fortran and python codes in this version need serious clean-up from junk.

* All Fortran codes are planned to be translated with f2py into python-importable libraries.

* Don't ever run MCMC runs (or Nest. Samp.) directly in the embedded Jupyter shell! This will freeze the GUI until is the MCMC done!
This is not a bug, simply the jupyter shell needs its thread and this is not done, yet. Hopefully, this will be fixed soon.
Please use the GUI navigation. 

* AT this point other bugs may occur! Please see: https://github.com/3fon3fonov/exostriker/issues
for work in progress issues.

* There is NO manual page at the moment.... This takes time, and at the moment I do
not have such! The GUI however, is, I believe, very intuitive, and with a trial and error one can figure out how everything works. 




 

Some credits:

* Some of the dynamical RV fitting routines are done by 
Xianyu Tan (the University of Hong Kong, now in Oxford) during
his Master studies with Man Hoi Lee (HKU). 

* Some of the Keplerian RV fitting routines and other N-body codes are initially written by Man Hoi Lee. 

* Jakub Morawski (Warsaw University) worked on the "RVmod" library as my intern student in the summer of 2018. His numerous contributions 
during the early stages of this project are highly appreciated.

* Grigorii Smirnov-Pinchukov helped to make the setup.py installer

* The AMD stability check function was donated by Stefan Dreizler (IAG, Germany). 


* The interactive plotting is done with a custom version of the "pyqtgraph": 

http://www.pyqtgraph.org/

* "GLS" and "MLP" periodograms are taken from Mathias Zechmeister's repo: 

https://github.com/mzechmeister/python

* "TLS" is taken from: 

https://github.com/hippke/tls

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

https://github.com/Axel-Erfurt/PyEdit2'

* N-body tests are performed using a custom version of the "Swift" N-body library,
modified by Man Hoi Lee (HKU) and Trifon Trifonov (MPIA).

https://www.boulder.swri.edu/~hal/swift.html

* Additionally, the Exo-Striker uses many "standard" Python libraries like 
"PyQt5", "matplotlib", "numpy", "scipy", "dill", "Jupyter", "qtconsole",
and more.



