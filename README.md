
**T**ransit and **R**adial velocity **I**nteractive **F**itting tool for **O**rbital analysis and **N**-body simulations : **The Exo-Striker** 

<p align="center">
  <img width="400" src="https://user-images.githubusercontent.com/44244057/51602645-d1680c80-1f07-11e9-8d5e-8ec9916552d5.png">
</p>
 

Very powerful and fast GUI tool for exoplanet orbital analysis. It uses a brang new RV fitting library called "RVmod", which can model the Stellar reflex motion caused by dynamically interacting planets in multi planetary systems. 

![new_es](https://user-images.githubusercontent.com/44244057/53755942-20667180-3eb8-11e9-9802-530618db7e7d.png)

**WARNING!** This tool is under active development and its functionality is enhanced on a daily basis! Therefore, although very unlikely, the version you download today may not be fully compatible with the version uploaded tomorrow! Use at your own risk!

Also, please keep in mind that this software is developed mostly for my needs and for fun. I hope, however, that you may find 
it capable to solve your scientific problems, too. For updates, follow my Twitter account https://twitter.com/3fon3fonov 

What works:

* Period search: GLS periodograms (RVs, act. data) & TLS (transit data).
* RVs Keplerian and Dynamical modeling + GP (only one GP kernel integrated so far).
* Transit modeling (so far only one dataset and not tested for 2+ transiting pl.)
* RVs + GP + Transit (in principle should work).
* RVs/Transit/GP MCMC sampling/modeling.
* RV auto-fit (RV automated planet-finder algortm).
* Long-term stability check of multiplanet systems using SyMBA, MVS, MVS with a GR precession.
* Interactive plots.
* RV vs. Activity time series correlation analysis/plots.
* Import/Export of work sessions and multi-sessions. 
* Export plots to a matplotlib window for further customization.
* Text editor/Calculator/Bash-shell widgets.
* Integrated Jupyter widget shell.
* Print the GUI screen into a .jpeg image (useful for sharing quick results, just like the image above)

What is to be implemented:

* Combined modeling with Astrometry.
* Variouse of minimization/sampling methods (SciPyOp implemented, but with only limited options, work in progress)
* Print all results and plots into a .tex ready environment.
* A large arsenal of N-body/dynamical simulation/analysis tools. 
* Documentation, Instructions and Video tutorials.
* For more "TBD" list see: "focus_matrix_TBFixed.doc".

If you use this tool and you find a bug or a problem, please report it!
