.. _gui:

GUI Layout
..........

Before starting the tutorials, it is kind to familiarize yourself with the
basic parameters that exostriker uses to determine the goodness of a fit but also
the ones that describe the planet orbits. All of these parameters can be seen on the 
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

Now, depending on the type of data that you are trying to fit, you have to choose
between Radial Velocities (RV data), Transits (Transit data) and TTVs (Transit timing variations) on the 
**Data area** and then add the data files.
