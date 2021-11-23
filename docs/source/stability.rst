.. _stability:

Stability analysis
..................

Performing an orbital evolution
===============================

At this point we can perform an orbital evolution, using the *SyMBA N-body
symplectic integrator*, in order to notice how the orbital parameters develop with time.
In the following tutorial the system consists of two planets. 


.. figure:: images/1orbitalevo.gif

   Running an *orbital evolution* with Keplerian model.


First the stellar parameters need to be distinguished, by changing the values
of *Stellar param.* on the I/O parameters section. Then add the maximum time of evolution
by clicking at *N-body* and the appropriate time step. Running orbital evolution (*Run orb. evol*) automatically redirects
to the *Orb. Evol* section, where the orbital parameters evolution is revealed.


.. figure:: images/dynamicalorb.gif

   Running an *orbital evolution* with Dynamical model.


Evaluating the *stability of a system* means that the orbital parameters have to be examined long-term (e.g 1Myr) and with a time
step of at least 100 points per orbit. For example, if the inner planet has a period of 200 days, then a time step of 2 
days is required. In case of planet–planet close encounters *SyMBA* automatically reduces the time step to ensure
an accurate simulation with high orbital resolution. *SyMBA* also checks for planet–planet or planet–star collisions or
planetary ejections and interrupts the integration if they occur. 

A planet is considered lost and the system unstable if, at any time:

* the mutual planet–planet separation is below the sum of their physical radii (assuming Jupiter mean density), i.e., the planets undergo collision.
* the star–planet separation exceeds two times the initial semimajor axis of the outermost planet, which we deﬁne as planetary ejection.
* the star–planet separation is below the physical stellar radius (R ≈ 0.03 au), which we consider a collision with the star.

All of these events are associated with large planetary eccentricities leading to crossing orbits, close planetary
encounters, rapid exchange of energy and angular momentum, and eventually instability. Therefore, these somewhat arbitrary
stability criteria are efﬁcient to detect unstable conﬁgurations and save CPU time.

!! Say some words about MVS & MVS_GR integrators !!

.. NOTE::
   In the **test arbitrary configuration** section tests with fixed values of planetary
   parameters can be run and their evolution can be examined.

----------------------------------------------------------------------------------------------------

Integrate mcmc/Nest. Samp. (Work in progress)





