.. _gui:

GUI Layout
..........

First it is kind to familiarize yourself with the
basic parameters that exostriker uses to determine the goodness of a fit but also
the ones that describe the planets orbitals. All of these parameters can be seen on the 
frontend of the GUI.

.. figure:: images/homepage.png
   :target: _images/homepage.png
   
   *Frontend of the GUI.*

------------------------------------------------------------------------------------------------------------

**Action menu bar**

* **File**: New session, open session/multi session, save session/multi session, open RVmod init file, open RVBank file, Quit.

* **Edit**: Push Jupyter var. to GUI, RV Auto fit, Reset mid-pannel buttons (dev-mode!).?

* **View**:
   **Print window GUI to file.**

   **Print f-test FAP**: ?

   .. Note::
      Each of the following "Get" actions is a command line to the Jupyter shell, except "Get All plots". You can
      modify its parameters, depending on the output.

   **Get LaTeX table with parameters**: Best fit parameters as a LaTeX table (.tex file on your local exostriker folder). 

   **Get LaTex table with priors**: Priors as a LaTeX table (.tex file on your local exostriker folder).

   **Get RV model**: Models information, BJD[days] / RV model[m/s] (.txt file on your local exostriker folder).

   **Get RV data**: Datapoints information, BJD[days]/ RVs[m/s]/ RVs errors[m/s]/ Number of dataset starting from 0 (.txt file on your local exostriker folder).

   **Get Orb evol.**: Orbital parameters, Time[years], Semimajor axis[au]/ Eccentricity/ Arg. of periastron[deg]/ Mean anomaly[deg]/ Inclination[deg]/ Longitude of the ascending node[deg] (.txt file on your local exostriker folder).

   **Get All plots**: All plots from the analysis are exported to "exported_plots" folder on your local exostriker folder.

   **Confidence interval table.**

* **Settings**: Change widget style, GUI & plots font.

* **Help**: Exostriker page on Github & Credits.

* **Control Sessions**: Navigate through sessions (new/copy/remove session).

**Statistical parameters**

* **rms**: Root-mean-square.
* **wrms**: Weighted root-mean-square.
* **χ**\ :sup:`2`: Chi-squared.
* **χ**\ :sup:`2` :sub:`reduced`\ : Chi-squared reduced.
* **lnL** : Log-likelihood function.
* **BIC** : Bayesian information criterion.
* **AIC** : Akaike information criterion.
* **N data**: Number of data/observations.
* **DOF**: Degrees of freedom.
* **AMD stable**: Checking the stability of a system (Green/Red).
* **More stat.info**: Provides information about the fit quality & RV data rms/wrms.

**Control parameters**

* **Simplex** : Fitting curves using the Simplex algorithm.
* **L-M** : Fitting curves using the Leveberg-Marquardt algorithm.
* **Keplerian** : Perform a Keplerian analysis.
* **Dynamical**: Perform a Dynamical analysis.
* **Initialize**: Fitting any change without optimizing (pressing Enter).
* **Fit**: Optimization parameter.
* **Run MCMC** : Triggers samples using the Markov chain Monte Carlo algorithm.
* **Run Nest.samp** : Triggers samples using the Nested sampling alogithm.
* **Run orbital evolution**: Perform parameter evolution.
* **RV auto fit**: RV automated planet-finder algorithm.

**Input/Output parameters**

* **P [d]**: Planets period.
* **K [m/s]**: RV amplitude. 
* **e**: Eccentricity.
* **ω [deg]**: Argument of periastron.
* **Ma [deg]**: Mean anomaly at the first observational epoch.
* **inc [deg]**: Inclination.
* **Ω [deg]**: Longitude of the ascending node.  
* **ώ [deg/yr]**: Rate of argument of periastron.
* **t**\ :sub:`0`\ **[d]**: Time of the first transit.
* **R**\ :sub:`pl`\ /**R**\ :sub:`*`\ : Planets radius in units of stellar radius.
* **a**\ :sub:`pl`\ /**R**\ :sub:`*`\ : Planets semimajor axis in units of stellar radius.
* **a [au]**: Semimajor axis.
* **m** [**M**\ :sub:`jup`\ ] : Mass. 
* **t**\ :sub:`ω`\ [**d**]: Time of periastron passage. 

----------------------------------------------------------------------------------------------------------

Data area
---------

Depending on the type of data that you are trying to fit, you can choose
between **Radial Velocities** (RV data), **Transits** (Transit data) and **TTVs** (Transit timing variations).

* RV data
   Load RVs, Include Offset/Jitter parameters, Choose a RV trend, Data options.

* Transit data
   Load Transits, Include Offset/Jitter parameters, Add trends, Add Limb-darkening parameters to model, Data options (detrend options).

* TTVs (Transit-Timing Variations)
   Load TTVs, Set the Epoch [BJD] of the first transit and the end of the model, Choose a time step in dyn. model.

* Activity
   Load Activity indicators from your local folder, Modify them and apply the changes.

* Limits and Priors
   Set limits to the parameters range.

------------------------------------------------------------------------------------------------------------

Help widgets area
-----------------

* Shells 
   *Exostriker* offers 3 command-line interpreters.

   **Jupyter**: A Qt-based console for working with Jupyter kernels. It provides a number of enhancements 
   only possible in a GUI, such as inline figures, proper multi-line editing with syntax highlighting, graphical
   calltips, and much more. For more information visit `qtconsole documentation`_.

   .. _qtconsole documentation : https://qtconsole.readthedocs.io/en/stable/

   **Bash shell**.

   **pqg shell**: PyQtGraph is a graphics and user interface library for Python. For more information visit `pyqtgraph documentation`_.

   .. _pyqtgraph documentation : https://pyqtgraph.readthedocs.io/en/latest/index.html

* Extra plots 
   In this section plots of the most prominent peaks of the RV data are displayed phase folded (phase diagrams).
   Additionally, periodograms of the RV data are included.   
   

      .. figure:: images/extraplots.gif
         :target: _images/extraplots.gif
         
         *Extra plots.*


* Data inspector
   Inspect the data on your local machine through the options **This computer** or **RVBank** and load them to exostriker. 


      .. figure:: images/datainspector.gif
         :target: _images/datainspector.gif
         
         *Data inspector.*


   The **RVBank** option offers data sets from *HARPS RVBank* and *HIRES NZP*. Choose between different types of **RV data** sets
   (RVs SERVAL + NZP correction etc.) and **Activity indicators** (CRX, dLW, .., etc.) 


      .. figure:: images/datainspector1.gif
         :target: _images/datainspector1.gif

         *RVBank.*

   Activity indicators can also be modified.

      .. figure:: images/modactivity.gif
         :target: _images/modactivity.gif

         *Activity indicators.*


* Text editor
   Through the *text editor* you can inspect and edit the data files. (Works for .dat, .tran, .vels extensions)

      .. figure:: images/texteditor.gif
         :target: _images/texteditor.gif

         *Text editor.*

* **Calculator**

* Stdout/Stderr
   This section provides information about the **version of the GUI** you are using. Also
   informs about the progress of all processes.

   .. WARNING::
      Before starting any project make sure that you run the latest version of *exostriker*. You can
      be updated about the latest version/updates of *exostriker* on exostriker's `github page`_ .
      
      .. _github page: https://github.com/3fon3fonov/exostriker

-------------------------------------------------------------------------------------------------------------------------------

Plotting widgets area
---------------------

* RV
   **RVs**: Radial velocities time series.
   
   **RVs o-c**: Radial velocities residuals.
   
   **GLS**: Generalized Lomb-Scargle periodogram of the initial signal.
   
   **GLS o-c**: Generalized Lomb-Scargle periodogram of the residual signal.
   
   **MLP**: Maximum Likelihood Periodogram. 
   
   **Window (DFT)**: Window function.
   
   For more information check the *Radial Velocity data* section.

* Transit
   **Tran.**: Transits time series.
   
   **Tran. o-c**: Transits residuals.
   
   **TLS**: Transit Least Squares of the initial signal.
   
   **TLS o-c**: Transit Least Squares of the residual signal.

   For more information check the *Transit data* section.

* TTV
   **TTVs**: TTVs.
   
   **TTVs o-c**: TTVs residuals.

* Activity
   **Time series**: Activity time series.
   
   **GLS**: Generalized Lomb-Scargle periodogram of the activity indicators.
   
   **Correlations**: Check the correlation between the RV data and the RV indicators.

* Sample correlation
   In this section graphs correlating the parameters samples that are generated through the MCMC or Nested Sampling
   methods are displayed. 
   
   (For more information check *Obtaining the best fit parameters uncertainties* section.)


   .. figure:: images/samplecor.gif
      :target: _images/samplecor.gif

      *Sample correlation.*


* Orb. Evol.
   Orbital parameters evolution time series. 
   
   (For more information check *Stability analysis* session.)
   
-----------------------------------------------------------------------------------------------------------------

Input/Output parameters area
----------------------------

* Planet param.
   Planetary parameters. The values change whenever a new model is applied. They can also get fixed.

* GP param.
   Gaussian processes parameters. ?

* Stellar param.
   Edit stellar parameters depending on your system.

* **Models param.**


   .. figure:: images/modelsparam.gif
      :target: _images/modelsparam.gif

      *Models parameters.*


   Edit RV model parameters.
   
   Choose between different minimizers in SciPy param. section.
   
   Configure the GLS/MLP/TLS options and MCMC/NS simulation parameters.
   
   Set the maximum number of planets Auto fit algorithm to look for. 

* Limits and Priors
   Set bounds to planetary parameters before the simulations.

* N-body
   Perform long-term stability check of multi-planet systems by setting the maximum time 
   of evolution. 
   
   Evolution of arbitrary planetary values can also be performed. 

   (For more information check the *Stability analysis* section.)

* **Plot opt.**


   .. figure:: images/plotopt.gif
      :target: _images/plotopt.gif

      *Plot options.*


   Customize the RV/Transit/TTVs graph, Change the size of the data points, their transparency (Alpha).
    
   Enable cross hair.
   
   Shift the planets phase signal.
   
   Configure the model.
   
   Configure GLS/MLP/TLS/DFT graphs (Select the number of peaks that 
   will be visible in the graphs).
   
   Show aliases in cross hair.
   
   Customize MCMC/NS sampling cornerplots and generate them.
