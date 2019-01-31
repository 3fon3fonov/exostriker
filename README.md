
**T**ransit and **R**adial velocity **I**nteractive **F**itting tool for **O**rbital analysis and **N**-body simulations : **The Exo-Striker** 

<p align="center">
  <img width="400" src="https://user-images.githubusercontent.com/44244057/51602645-d1680c80-1f07-11e9-8d5e-8ec9916552d5.png">
</p>
 

Very powerful and fast GUI tool for exoplanet orbital analysis. It uses a brang new RV fitting library called "RVmod", which can model the Stellar reflex motion caused by dynamically interacting planets in multi planetary systems. 

![exo_striker_new](https://user-images.githubusercontent.com/44244057/52046917-38a94080-2548-11e9-87f6-1873e167eb42.png)

**WARNING!** This tool is under active development and its functionality is enhanced on a daily basis! Therefore, although very unlikely, the version you download today may not be fully compatible with the version uploaded tomorrow! Use at your own risk!

Also, please keep in mind that this software is developed mostly for my needs and for fun. I hope, however, that you may find 
it capable to solve your scientific problems, too. For updates, follow my Twitter account https://twitter.com/3fon3fonov 

What works:

* Periodograms.
* RVs Keplerian and Dynamical modeling + GP (only one GP kernel integrated so far).
* RVs MCMC sampling/modeling.
* RV auto-fit (automated planet-finder algortm).
* Transit modeling (so far only one planet and one dataset)
* Long-term stability check of multiplanet systems using SyMBA, MVS, MVS with a GR precession.
* Interactive plots.
* Import/Export of work sessions and multi-sessions. 
* Export of plots to a matplotlib window for further customization.
* Text editor/Calculator/Bash-shell widgets.
* Integrated Jupyter widget shell.
* Print the GUI screen into a .jpeg image (useful for sharing quick results, just like the image above)

What is to be implemented:

* Combined modeling of data (Transit + RV + Astrometry +  GP/Moving average modeling, at once). 
* Variouse of minimization/sampling methods.
* Some more GUI plots and plot functionality.
* GUI accsess to parameter limits and priors (there, but not connected with RVmod, yet)
* Print all results and plots in a .tex ready environment
* Large arsenal of N-body/dynamical simulation/analysis tools. 
* Documentation, Instructions and Video tutorials.
* For more detailed TBD list see: "focus_matrix_TBFixed.doc".

If you use this tool and you find a bug or a problem, please report it!
