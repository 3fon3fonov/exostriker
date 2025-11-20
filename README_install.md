# Pre-installation - Anaconda

It is recommended to use Anaconda, especially for Windows users. For this:

- First, download anaconda [here](https://www.anaconda.com/products/individual#download-section) and execute the script.

    - For Linux users, if needed, export anaconda to your path by adding the line `export PATH="/home/USERNAME/anaconda3/bin:$PATH"` (changing the `USERNAME` to your user) to the `~/.bashrc` file. Close and open a new terminal to update the path.

    - If desirable, disable auto activation of anaconda with `conda config --set auto_activate_base false`.

- Create Exostriker environment with `conda create python -n exostriker` (you may change the environment name after the `-n` flag);

- Then, to activate the environment, type `conda activate exostriker`;

## Additional step

This is especially necessary for Windows and can also be needed for some new M-chip macOS users.

Depending on the type of GFortran installation, anaconda might not be able to find it. In order to satisfactorily compile the Fortran code of Exostriker when using Anaconda, install the Gfortran directly from the Anaconda repository:

```bash
conda install conda-forge::gfortran
```

Obs: this installation needs to be done with the exostriker environment acitvated.

# Instalation

There are many options to run Exo-Striker. One may try an Anaconda environment, install it from pip, or run directly the source code from the repository.  

All the methods listed below work in Linux, Windows, or MacOS if you have python and gfortran installed.

## Getting the pip version (recommended)
You can easily get the pip version by typing in the terminal:

```
pip install exostriker
```
or for e.g., python3.10

```
python3.10 -m pip install exostriker
```

**This is all.** Below other installation methods are described.

## Cloning from GitHub

You can directly clone from the GitHub page with:

```bash
git clone https://github.com/3fon3fonov/exostriker.git
```

After cloning, you can install it with pip or run it without any installation (see the Execution section).  

To build the package and install it with pip from the source, type:  

```bash
pip install build
cd exostriker
python3 -m build
pip install .
```

## Adding the PATH (for non-Anaconda installations)
If you do not use Anaconda, you need to add the path of the pip scripts into the system path variable if you have not already done so.  
**This step is not necessary if you are running without installation.**

### Linux
For Linux users, you can add the following line in the end of the file `~/.bashrc`, changing `USER` to your own user.

```bash
export PATH="${PATH}:/home/USER/.local/bin/"
```

The locale of the executables can slightly change for different Linux distributions and Python installations. If this does not work, you can try different paths, such as `/usr/local/bin` and others.


### Windows
For Windows users, you need to open the menu and search for *path*, click in *Edit the system environment variables*, at the bottom right click in *Environment Variables...*, in the tab *System variables* (attention: not the *User variables for Username*), look for the variable *Path* and click on *Edit*. Add a new line with one of the following (check the Python location first):

```bash
C:\Users\Windows\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.11_qbz5n2kfra8p0\LocalCache\local-packages\Python311\Scripts
```

Or:

```bash
C:\Users\USERNAME\AppData\Local\Programs\Python\Python311\Scripts
```

**Be aware of different Python versions, the path will change also. Always verify if the current path exists.**

## Step-by-Step for Windows (without Anaconda)
Obs: this step-by-step guide was made using a fresh new install of Windows 10.

Python: Go to [Python Windows Releases](https://www.python.org/downloads/windows/) and download a stable release by clicking on it and downloading the Windows installer (64-bit) Recommended. After downloading, execute the files and follow the installation process.
    - Add the Python folder to the path.
        - Open the menu and search for *path*, click in *Edit the system environment variables*, at the bottom right click in *Environment Variables...*, in the tab *System variables* (attention: not the *User variables for Username*), look for the variable *Path* and click on *Edit*. Add two new lines:
            - `C:\Users\USERNAME\AppData\Local\Programs\Python\Python311`
            - `C:\Users\USERNAME\AppData\Local\Programs\Python\Python311\Scripts`
            - Change `USERNAME` with your username and `Python311` to the actual version that you installed.
    - Go to the python folder (ex: `C:\Users\USERNAME\AppData\Local\Programs\Python\Python311`) and create a copy of the `python.exe` file in the same location with the name `python3.11.exe`
- Gfortran: Download the *x64* version of [Winlibs](https://winlibs.com).
    - Choose the x64 POSIX **with** all the libraries;
    - After downloading, decompress the file and move the `mingw64` folder to `C:\Program Files\`;
    - Now add the folder `bin` of the `mingw64` folder to the path:
         - Use the same steps as before, but now add the line: `C:\Program Files\mingw64\bin`.
- Microsoft Visual C++ for the package TTVFast: [Download](https://visualstudio.microsoft.com/visual-cpp-build-tools/) it and install the *Microsft Build Tools*, in the installation process, check the box for the *Desktop development with C++* (6.52 GB) and then click to install.
- Exostriker: In CMD type `pip install exostriker` or follow the above steps to install it.

**WARNING** if your USERNAME contains a "space", e.g., "Trifon Trifonov", the installation of RVMod will fail!!!! No space allowed as of Version 0.93. TBFixed!

# Execution

## With Installation
You can create menu entries for Exo-Striker by typing in the terminal or cmd:  

```bash
exostriker-desktop-create
```

Then, search in the menu for *Exostriker*, and in Windows, a Desktop link will also be created.  
However, if you do not want to create menu entries, you can launch it by typing `exostriker` or `python—m exostriker` in the terminal or cmd.

## Without Installation
It is possible to execute the Exo-Striker without doing any installation! E.g.: 

```bash
git clone https://github.com/3fon3fonov/exostriker.git
cd exostriker
python3 -m exostriker
```


Simply execute the file `gui.py` and it will power the GUI.

# Uninstall
If you did not install with pip (both directly or by cloning), you could just erase the *exostriker* folder.  

Otherwise, first, remove the menu entries (if you have added it) and then uninstall:

```bash
exostriker-desktop-remove
pip uninstall exostriker
```
