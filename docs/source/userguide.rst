.. _userguide:

Installation / Loading the GUI
..............................

To install the Exo-Striker, first please make sure that you have:

* csh
* gfortran
* pip

Also, please only use **Python 3** for installing the Exo-Striker! The tool works with Python 2, 
but since Python 2 is no longer maintained (since Jan 1, 2020), I do not assist in case of problems.    

Currently there are **three ways to install/run** the Exo-Striker:
   
* The **simplest way** to "git clone":    

$ git clone https://github.com/3fon3fonov/exostriker  

$ cd exostriker

and then **load the gui**:    

$ python3 exostriker_gui.py 

or 

$ ./exostriker_gui.py

Generally, you do not need to install anything if you already have all the required dependencies for the tool to run.
For the dependency list, see the "setup.py" file. The Exo-Striker will automatically compile the Fortran77 code for you
at the first start of the program and will keep you updated if the source code was updated (if you regularly "git pull").    

--------------------------------------------------------------

* The **second way** to install the tool from the source:    

$ git clone https://github.com/3fon3fonov/exostriker and then:    

$ cd exostriker

$ python3 setup.py install    

This will install the tool in your system.     
Then,      

$ python3 exostriker_gui.py

This should start the Exo-Striker.  

----------------------------------------------------------------

* and **last**, you can try pip install:    

$ pip install git+https://github.com/3fon3fonov/exostriker    

This will install the tool in your system.    
Then,     

$ python3 exostriker_gui.py


This should start the Exo-Striker.

---------------------------------------------------------------

If you **still cannot boot** the tool after a 'successful' installation, please try:

$ python3 exostriker_gui.py -debug 

or 

$ exostriker -debug 

(depending on how you use the tool)

Then, copy the output error, and please open a 'GitHub' issue. Otherwise, all possible problems/bugs/crashes
will be displayed on the 'stdout/stderr' tab of the tool. If you use this tool, and you find a bug or a problem,
please report it!    

The instructions above will install the following **dependencies**: 

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

if not, install them manually. e.g pip install jupyter, ..., etc.

-----------------------------------------------------------------------------------------------------

OS related comments
...................

**LINUX**

The above instructions usually work without any problem on Linux OS.

Though, the installation might be easier (no problem with missing libraries) if you have the Anaconda package already
installed on your system. See anaconda documentation : https://docs.anaconda.com/anaconda/install/
  
For full functionality on Linux, you will also need to install:

* rxvt (optional, better bash shell)

-----------------------------------------------------------------

**MAC OS** 

You will need to install "homebrew", which requires "sudo" privileges. 
According to 

https://brew.sh/

To install "homebrew", it should be as simple as this:

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

try it out in bash shell terminal.


If you have MAC OS 10.14+ it is likely to have problems with the numpy headers.

if you see an error like this during the installations (mostly during the "batman-package" build):



    c_src/_nonlinear_ld.c:21:10: fatal error: 'numpy/arrayobject.h' file not found
    #include "numpy/arrayobject.h"
             
    1 error generated.
    error: command 'clang' failed with exit status 1



Please do the following:

export CFLAGS="-I /usr/local/lib/pythonX.X/site-packages/numpy/core/include $CFLAGS"

(Where X is your Python version, e.g., 3.6, 3.7 to 3.8)

and then re-run your Exo-Striker installation (via pip or setup.py).
 

----------------------------------------------------------------------------------

**WINDOWS 10**

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
 
------------------------------------------------------------------------------------------------------------

Some known problems
...................  
 
For work in progress issues see:
 
https://github.com/3fon3fonov/exostriker/issues
 
------------------------------------------------

Reporting an issue
..................

If you run into issues or bugs, do not hesitate to report it (**New issue**) on the GitHub repository (https://github.com/3fon3fonov/exostriker/issues) 
or send a PM to trifonov@mpia.de.

Feedback and help in further development will be highly appreciated! A wish-list with your favorite tools and methods to be implemented is also welcome!

