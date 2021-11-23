.. _gui:

GUI Layout
..........

Before starting the tutorials, it is kind to familiarize yourself with the
basic parameters that exostriker uses to determine the goodness of a fit but also
the ones that describe the planets orbits. All of these parameters can be seen on the 
homepage of exostriker.

.. figure:: /images/homepage.png
   

   *Home Page of the GUI*

**Statistical parameters**

* **rms**: root-mean-square.
* **wrms**: weighted root-mean-square.
* **χ**\ :sup:`2`: chi-squared.
* **χ**\ :sup:`2` :sub:`reduced`\ : chi-squared reduced.
* **lnL** : likelihood function.
* **BIC** : Bayesian information criterion.
* **AIC** : Akaike information criterion.
* **N data**: number of data/observations.
* **DOF**: degrees of freedom. 

**Control parameters**

* **Simplex** : fitting curves using the Simplex algorithm.
* **L-M** : fitting curves using the Leveberg-Marquardt algorithm.
* **Keplerian** : perform a Keplerian analysis.
* **Dynamical**: perform a Dynamical analysis.
* **Initialize**: ?
* **Fit**: optimization parameter.
* **Run MCMC** : generates samples using the Markov chain Monte Carlo algorithm.
* **Run Nest.samp** : generates samples using the Nested sampling alogithm.
* **Run orbital evolution**: perform parameter evolution.
* **RV auto fit**: apply a curve to data.

**Input/Output parameters**

* **P [d]**: planets period.
* **K [m/s]**: planets RV amplitude. 
* **e**: eccentricity of the orbit.
* **ω [deg]**: argument of periastron.
* **Ma [deg]**: mean anomaly at the first observational epoch.
* **inc**: inclination of the orbit.
* **Ω [deg]**: longitude of the ascending node.  
* **ώ [deg/yr]**: rate of argument of periastron.
* **t**\ :sub:`0`\ **[d]**: time of the first transit.
* **R**\ :sub:`pl`\ /**R**\ :sub:`*`\ : planet radius in units of stellar radius.
* **a**\ :sub:`pl`\ /**R**\ :sub:`*`\ : planet semimajor axis in units of stellar radius.
* **a [au]**: semimajor axis.
* **m** [**M**\ :sub:`jup`\ ] : planets mass. 
* **t**\ :sub:`ω`\ [**d**]: ? 

----------------------------------------------------------------------------------------------------------

Data area
---------

Depending on the type of data that you are trying to fit, you can choose
between Radial Velocities (RV data), Transits (Transit data) and TTVs (Transit timing variations).

* RV data
   Load RVs, include Offset/Jitter parameters, choose a RV trend, bin data.

* Transit data
   Load Transits, include Offset/Jitter parameters, add Limb-darkening parameters to the model,
   detrend the data.

* TTVs (Transit-Timing Variations)
   ?

* Activity
   ?

* Limits and Priors
   Set limits to the parameters range.

Help widgets area
-----------------

* Shells 
   *Exostriker* offers 3 command-line interpreters.

   **Jupyter**: A Qt-based console for working with Jupyter kernels. It provides a number of enhancements 
   only possible in a GUI, such as inline figures, proper multi-line editing with syntax highlighting, graphical
   calltips, and much more. For more information visit `qtconsole documentation`_.

   .. _qtconsole documentation : https://qtconsole.readthedocs.io/en/stable/

   **Bash shell**: Work in progress.

   **pqg shell**: PyQtGraph is a graphics and user interface library for Python. For more information visit `pyqtgraph documentation`_.

   .. _pyqtgraph documentation : https://pyqtgraph.readthedocs.io/en/latest/index.html

* Extra plots 
   In this section plots of the most prominent peaks of the RV data are displayed phase folded (phase diagrams).
   Additionally, periodograms of the RV data are included.   
   
* Data inspector
   Inspect the data on your local machine through the options *This computer* or *RVBank* and load them to exostriker. 

* Text editor
   Through the *text editor* you can inspect and edit the data files. (Works for .dat, .tran, .vels extensions)

* **Calculator**

* Stdout/Stderr
   This section provides information about the version of the GUI you are using. Also
   informs about the progress of all processes.

   .. WARNING::
      Before starting any project make sure that you run the latest version of *exostriker*. You can
      be updated about the latest version of *exostriker* on exostriker's `github page`_ .
      
      .. _github page: https://github.com/3fon3fonov/exostriker

Plotting widgets area
---------------------

* RV
   **RVs**: Radial velocity graph.
   
   **RVs o-c**: Radial velocity residuals graph.
   
   **GLS**: Generalized Lomb-Scargle periodogram of the initial signal. Options including Cross hair & jitter to graph.
   
   **GLS o-c**: Generalized Lomb-Scargle periodogram of the residual signal. Adopt best parameter option ??
   
   **MLP**: Maximum Likelihood Periodogram. 
   
   **Window (DFT)**: ?

    
* Transit
   **Tran.**: Transit graph.
   
   **Tran. o-c**: Transit residuals graph.
   
   **TLS**: Transit Least Squares of the initial signal.
   
   **TLS o-c**: Transit Least Squares of the residual signal.

* TTV
   **TTVs**: TTVs graph.
   
   **TTVs o-c**: TTVs residuals graph.

* Activity
   **Time series**: ?
   
   **GLS**: ?
   
   **Correlations**: ?

* Sample corr.
   ?

* Orb. Evol.
   Orbital parameters evolution graphs.
   
Input/Output parameters area
----------------------------

* Planet param.
   Planetary parameters. The values change whenever a model is fitted. The values can also be fixed.

* GP param.
   Gaussian processes parameters. ?

* Stellar param.
   Edit stellar parameters depending on your system.

* **Models param.**


   .. image:: images/modelsparam.gif


   Edit RV model parameters.
   
   Choose between different minimizers in SciPy param. section.
   
   Configure the GLS/MLP/TLS options and MCMC/NS simulation parameters.
   
   Set the maximum number of planets Auto fit to look for. 

* Limits and Priors
   Set bounds to planetary parameters before the simulations.

* N-body
   Perform long-term stability check of multi-planet systems by setting the maximum time 
   of evolution. 
   
   Evolution of arbitrary planetary values can also be performed. 

* **Plot opt.**


   .. image:: images/plotopt.gif


   Customize the RV/Transit/TTVs graph (Change the size of the data points, their transparency (Alpha).
    
   Enable cross hair.
   
   Shift the planets phase signal.
   
   Configure the model.
   
   Configure GLS/MLP/TLS/DFT graphs (Select the number of peaks that 
   will be visible in the graphs).
   
   Show aliases in cross hair.
   
   Customize MCMC/NS sampling cornerplots and generate them.
