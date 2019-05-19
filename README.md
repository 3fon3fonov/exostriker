
**T**ransit and **R**adial velocity **I**nteractive **F**itting tool for **O**rbital analysis and **N**-body simulations : **The Exo-Striker** 

<p align="center">
  <img width="400" src="https://github.com/3fon3fonov/trifon/blob/master/lib/33_striker.png">
</p>
 

Very powerful and fast GUI tool for exoplanet orbital analysis. It uses a brand new RV fitting library called "RVmod", which can model the Stellar reflex motion caused by dynamically interacting planets in multi planetary systems. 

![new_es](https://user-images.githubusercontent.com/44244057/57973277-ccdbce00-79a6-11e9-9930-60a40beb0b04.png)

**WARNING!** This tool is under active development and its functionality is enhanced on a daily basis! Therefore, although very unlikely, the version you download today may not be fully compatible with the version uploaded tomorrow! Use at your own risk!

Also, please keep in mind that this software is developed mostly for my needs and for fun. I hope, however, that you may find 
it capable to solve your scientific problems, too. For updates, follow my Twitter account https://twitter.com/3fon3fonov 

What works:

* Period search: GLS periodograms (RVs, act. data) & TLS (transit data).
* Keplerian and Dynamical RV modeling. 
* Transit modeling (not really tested for 2+ transiting planets and more datasets)
* GP modeling (only SHO and Rot. GP kernels integrated so far).
* Joint RVs + RV GP + Transit best-fit optimization.
* Joint RVs + RV GP + Transit MCMC/Nested Sampling sampling (via "emcee" & "dynesty").
* RV auto-fit (RV automated planet-finder algorithm).
* Long-term stability check of multi-planet systems using SyMBA, MVS and MVS with an GR precession.
* Various of minimization methods (via "SciPyOp").
* Interactive plots.
* RV vs. Activity time series correlation analysis/plots.
* Import/Export of work sessions and multi-sessions. 
* Export plots to a matplotlib window for further customization.
* Export ready to use Latex tables with best-fit parameters, errors and statistics. 
* Text editor.
* Bash-shell (Linux only).
* Integrated Jupyter widget shell.
* Print the GUI screen into a .jpeg image (useful for sharing quick results, just like the image above)

What is to be implemented:

* More GP kernels (work in progress). 
* GP on the transit photometry  (work in progress)
* Transit auto-fit (transit automated planet-finder algorithm).
* Binary/Triple star modeling mode.
* Combined modeling with Astrometry.
* A larger arsenal of N-body/dynamical simulation/analysis tools (REBOUND to be included). 
* Documentation, Instructions and Video tutorials.
* For more "TBD" list see: "focus_matrix_TBFixed.doc".

If you use this tool and you find a bug or a problem, please report it!
A wish-list with your favorite tools and methods to be implemented is also welcome.



