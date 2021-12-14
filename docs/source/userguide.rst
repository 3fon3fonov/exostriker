.. _userguide:

Installing / Loading the GUI
............................

To **install** the Exo-Striker, first please make sure that you have:

* csh
* gfortran
* pip

Also, please only use **Python 3** for installing the Exo-Striker! The tool works with Python 2, 
but since Python 2 is no longer maintained (since Jan 1, 2020), I do not assist in case of problems.    

-------------------------------------------------------------------

Currently there are **three ways to install/run** the Exo-Striker:
   
* The **simplest way** is to **git clone**:    

**$** git clone https://github.com/3fon3fonov/exostriker  

**$** cd exostriker

and then **load the gui**:    

**$** python3 exostriker_gui.py 

**or** 

**$** ./exostriker_gui.py

--------------------------------------------------------------

* The **second way** to install the tool from the source:    

**$** git clone https://github.com/3fon3fonov/exostriker and then:    

**$** cd exostriker

**$** python3 setup.py install    

This will install the tool in your system.     
Then,      

**$** python3 exostriker_gui.py

This should start the Exo-Striker !!!  

----------------------------------------------------------------

* and **last**, you can try **pip install**:    

**$** pip install git+https://github.com/3fon3fonov/exostriker    

This will install the tool in your system.    
Then,     

**$** python3 exostriker_gui.py

This should start the Exo-Striker !!!

----------------------------------------------------------------------------

Generally, **you do not need to install anything if you already have all the required dependencies** for the tool to run.
The dependency list, including the version of each dependency, can also be seen at the "setup.py" file.

.. IMPORTANT::
   List of required **dependencies** : 

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

You can also install these manually, e.g **pip install jupyter**, ..., etc.

---------------------------------------------------------------

If you **still cannot boot** the tool after a 'successful' installation, please try:

**$** python3 exostriker_gui.py -debug 

or 

**$** exostriker -debug 

(depending on how you use the tool)

Then, copy the output error, and please **open a GitHub issue**. Otherwise, all possible problems/bugs/crashes
will be displayed on the 'stdout/stderr' tab of the tool. If you use this tool, and you find a bug or a problem,
please report it!    

-----------------------------------------------------------------------------------------------------

OS related comments
...................

**LINUX**

The above instructions usually work without any problem on Linux OS.

Though, the installation might be easier (no problem with missing libraries) if you have the Anaconda package already
installed on your system. See anaconda documentation : https://docs.anaconda.com/anaconda/install/
  
For full functionality on Linux, **you will also need to install**:

* rxvt (optional, better bash shell)

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

Reporting an issue
..................

If you run into issues or bugs, do not hesitate to report it **(New issue)** on the `GitHub repository`_ or send a PM to trifonov@mpia.de.

.. _GitHub repository: https://github.com/3fon3fonov/exostriker/issues 

**Feedback** and **help** in further development will be highly appreciated! A wish-list with your favorite tools and methods to be implemented is also welcome!

