.. _rvs:

Adding the RV data
..................


*There are 3 ways to load RV data on exostriker*.

*  the **first** is through the **Data area**.


.. image:: images/1addrv.gif


* the **second** is through the **Data inspector** on the *Help widgets area*.


.. image:: images/2addrv.gif   


* the **third** is also through the **Data inspector** but from the RVBank option.


.. image:: images/3addrv.gif


---------------------------------------------------------------------------------------

Fitting the RV data
===================

Now lets perform a **simple fit** to data to showcase how exostriker works and we will
explain the remaining GUI options as a part of each tutorial.

*  Inspecting for **periodicities** on the RV graph.

In this tutorial we will fit the RV dataset with the **Keplerian model**. This model uses the radial 
velocity equation to fit the data. It's a model of motion of the star as expected from Kepler's laws. 


.. image:: images/checkperiodrv.gif
  

By *fitting* the data (Fit option), exostriker fixes the Offset parameter and also applies the jitter (jitter is only included on Simplex).
A curve is applied on the data by the *RV Auto fit* option. It is proper to fit between the 
Simplex & L-M algorithms and notice how the Statistical parameters change in order to acquire the best fit.
Now we notice that the RV graph is a periodic function. On the *GLS* (*generalized Lomb-Scargle periodogram*) section,
the periodogram of our initial RV periodic function is presented, including the main peaks that it is 
consisted of. Every significant peak needs investigation of whether it is a planet or stellar activity.
In this tutorial the initial signal is consisted of two planets.
A phase diagram of each significant peak is displayed by selecting *Extra Plots*. 

* Checking the residual signal. 

The *RVs o-c* & *GLS o-c* panels represent the residual RV signal and periodogram each time
a planet is included in the analysis. 

If the statistical parameters of a model are maximized/minimized then, probably 
the *best Keplerian fit* is adopted. 

----------------------------------------------------------------------------------------

Multiplanetary systems
======================

In case of *multiplanetary systems* it is kind to consider the planet masses and distances from each other (close orbits).
Massive planets with close distances from the host star, will surely interact with each other due to gravity. 
Then a further investigation using the **Dynamical model** is necessary. That will take into account the
gravitational interactions between the massive bodies by intergrating the equations of motion using the 
*Gragg-Bulirsch-Stoer* method.


.. image:: images/dynamicalrv.gif


Before enabling the *Dynamical option* make sure that the orbital parameters that are acquired so far 
correspond to the *best Keplerian fit*, because they will be used as a first guess for this *Dynamical fit*.
The next thing to be noticed is that the orbital parameters inclination (i) and the longitude 
of the ascending node (Ω) become available. The dynamical model has the advantage of being able to fit for 
mutually inclined orbits. For the purposes of this tutorial we assume edge-on coplanar
orbits (i=90, Ω=0) for consistency with the unperturbed Keplerian frame and in order to work with minimum
dynamical masses.

The difference in the phase diagrams between the two models is significant! 

------------------------------------------------------------------------------------------

Performing an orbital evolution
===============================

At this point we can perform an orbital evolution, using the *SyMBA N-body
symplectic integrator*, in order to notice how the orbital parameters develop with time.


.. figure:: images/1orbitalevo.gif

   Running an *orbital evolution* with Keplerian model.


First the stellar parameters need to be distinguished, by changing the values
of *Stellar param.* on the I/O parameters section. Then add the maximum time of evolution
by clicking at *N-body*. Running orbital evolution (*Run orb. evol*) automatically redirects
to the *Orb. Evol* section, where the orbital parameters evolution is revealed.


.. figure:: images/dynamicalorb.gif

   Running an *orbital evolution* with Dynamical model.


----------------------------------------------------------------------------------------------------

Obtaining the best fit parameters uncertainties I
=================================================

In order to perform a parameter distribution analysis and estimate the best fits parameters uncertainties,
they need to coupled with a *Markov chain Monte Carlo* (**Run MCMC**) sampling scheme using the *emcee
sampler*. 


.. image:: images/mcmc.gif


First, fill the parameters in *MCMC param.* (Models param.) on I/O parameters panel.

* **Burning phase samp.** : First steps in the MCMC chain to let the walkers explore the parameter space (At least 1000).
* **MCMC phase samp.** : Represents the total amount of samples (At least 5000).
* **N threads/CPUs** : Number of CPUs from your local machine that will be used for this process.
* **Init. Gaussian Ball** : How far from the starting point the sampler will start producing samples.
* **N walkers factor** : Each walker will produce a different chain (N walkers factor * DOF). 

Then, select between the options in *Adopt MCMC param. as*, in this tutorial we select the *best
maximum likelihood*.

Now, everything is ready for the MCMC process to start. A different set of parameters are being tested, **starting
from the best ones that are already obtained**, and for each set the likelihood is calculated. If a better lnL is
found then its parameters are counted as errors on the initials (from the best fit). 

For the purposes of this tutorial the values of *Burning phase & MCMC phase samples* are low, for the process to finish
at a reasonable time. Sometimes the MCMC process can take from hours to days, depending on the amount of samples and the
dimensions of the system (N data minus DOF). You can always check the MCMC progress at the **Stdout/Stderr** on the Help 
widgets area.

**Keep in mind** that setting bounds on each parameter on the *Limits and Priors*, on the I/O parameters panel, will save
computational time. In this tutorial we set the boundaries empirically. 


.. image:: images/mcmc1.gif


When the process is over, the samples are ready! *Go to "Make Cornerplot"* option redirectes to the *Plot options*. 
There a customization of the cornerplot is offered but also an option to include/exclude parameters from
being printed. By pressing *Make cornerplot* the final results are extracted on your local exostriker folder as a pdf
format. 


.. figure:: images/histograms.png

   Final histograms.


In our case the final histograms show the posterior MCMC distribution of the ﬁtted parameters with a dynamical modeling scheme whose orbital
conﬁguration is edge-on and coplanar. The histogram panels in the figure provide a comparison between the probability density distribution
of the overall MCMC samples for each fitted parameter. The two-dimensional parameter distribution panels represent all possible parameter 
correlations with respect to the best dynamical ﬁt (starting point), whose position is marked with blue lines. In color red all the samples
are represented. Also The red 2D contours are constructed from the overall MCMC samples and indicate the 68.3%, 95.5%, and 99.7% conﬁdence
levels (i.e., 1σ, 2σ, and 3σ).

Further information about *emcee sampler* and its modes can be found on `emcee documentation`_.

.. _emcee documentation: https://emcee.readthedocs.io/en/stable/

Obtaining the best fit parameters uncertainties II
==================================================

A modern way to perform parameter distribution analysis is by the *Nested Sampling technique* (**Run Nest. samp**)
using the *DYNESTY* sampler. Nested Sampling algorithm is a computational approach to the Bayesian statistics
problems of comparing models and generating samples from the posterior distributions.


.. image:: images/ns1.gif


First, fill the parameters in *Nest. Samp. param.* (Models param.) on the I/O parameters panel aswell.

* **Static** : Static sampling, the number of live points remain constant.
* **Dynamic** : Dynamic sampling, the number of live points vary during runtime.
* **Dynamic samp. opt.** : Some sampling options. 
* **N threads/CPUs** : Number of CPUs from your local machine that will be used for the process.
* **dlogz stop** : Stopping criterion.
* **Live points factor** : Number of live points used (live points factor * DOF).

Then, a selection between the options in *Adopt MCMC param. as* is required, in this tutorial we select the *best
maximum likelihood*. *Dynesty* also supports a number of options for bounding the target distribution (*Dynesty bound opt.*).

Unlike MCMC, Nested Sampling starts by randomly sampling from the entire parameter space specified by the prior.
This is not possible unless the priors are “proper” (i.e. that they integrate prior to 1). So while normal
priors spanning (−∞, +∞) are fine, uniform priors spanning the same range are not and must be bounded 
(*Limits and Priors*).

In this tutorial we run *dynesty* using the default settings whenever possible to highlight performance in a “typical” 
use case. You can always check the progress of the sampling at *Stdout/Stderr* on Help widgets area. 

.. image:: images/ns2.gif

When the process is over, the samples are ready! *Go to "Make Cornerplot"* option redirectes to the *Plot options*. 
There a customization of the cornerplot is offered but also an option to include/exclude parameters from
being printed. By pressing *Make cornerplot* the final results are extracted on your local exostriker folder as a pdf
format. 


.. figure:: images/ns1.png

   Final histograms.


The final histograms show the posterior distribution of the ﬁtted parameters with a dynamical modeling scheme whose orbital
conﬁguration is edge-on and coplanar. The histogram panels in the figure provide a comparison between the probability density distribution
of the overall NS samples for each fitted parameter. The two-dimensional parameter distribution panels represent all possible parameter 
correlations with respect to the best dynamical ﬁt (starting point), whose position is marked with red lines. In color blue all the samples
are represented. Also The blue 2D contours are constructed from the overall NS samples and indicate the 68.3%, 95.5%, and 99.7% conﬁdence
levels (i.e., 1σ, 2σ, and 3σ).

Further information about *DYNESTY* and its modes can be found on `dynestys documentation`_.

.. _dynestys documentation: https://dynesty.readthedocs.io/en/latest/index.html
