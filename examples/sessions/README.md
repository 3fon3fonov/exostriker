# GetModel synthetic session

`GetModel_giant_planets.ses` is a fitted synthetic two-planet tutorial session,
with adjacent portable RV data `GetModel_giant_planets.vels`. The notebook
`../GetModel_synthetic_data.ipynb` repairs its relative data path after loading.

- Solar mass/radius star; injected planets 1.0 and 1.1 Mjup, P=320/680 days.
- Eccentricities .08/.12, omega=40/110 deg, inclination=89.98 deg, node=0,
  Rp/Rstar=.10045. These are fixed tutorial assumptions.
- 180 synthetic Keplerian RV observations, sigma=1.5 m/s, offset=12 m/s.
  Observing schedule seed=81; GetModel measurement-noise seed=42.
- Native LM maximum-likelihood fit: K/P/M0 per planet plus instrument offset;
  Gaussian Fisher covariance retained. Best lnL approximately -334.768.
- Real MCMC: 28 walkers, 200 burn-in steps discarded, 600 retained steps per
  walker, yielding 16,800 retained 7-parameter rows. This is a short tutorial
  chain, not a publication convergence claim.
- Real static nested sampling: 112 live points, rwalk, stopping dlogz=.2,
  3,803 aligned raw rows retained with log weights. A known-period/phase local
  tutorial search box avoids distant aliases. Raw maximum lnL about -335.093.

MCMC used the initial wider bounds: K=1..80 m/s, P=0.8..1.2 times each injected
period, M0=0..360 deg, instrument offset=-100..100 m/s. NS used:

| Planet | K [m/s] | P [day] | M0 [deg] |
| --- | --- | --- | --- |
| b | 15..45 | 300..340 | 5..35 |
| c | 15..40 | 640..720 | 280..320 |

The priors differ, so evidence/credible intervals are conditional on their
respective priors; this is not a controlled evidence comparison. Session
metadata `getmodel_demo` records the synthetic origin and NS bounds. The
notebook shows how to regenerate a fresh LM session and optionally run longer
MCMC/NS calculations without overwriting this fixture. Samplers are stochastic;
the bundled rows allow reproducible propagation without rerunning inference.

Transit/astrometry were not fitted. Their models/bands reuse RV orbital samples
with assumed fixed radii, inclinations and stellar/astrometric properties. One
shared stellar parallax is used for both standard motion and angular reflex
scaling. No real Gaia scanning law or real observations are implied.

Posterior storage uses NumPy arrays in standard-library `SimpleNamespace`
containers. Live sampler callbacks/classes are intentionally excluded so a
Python 3.10 fixture does not execute saved Python 3.10 bytecode in another
interpreter. Retained chains, log probabilities, acceptance fractions, raw NS
weights and posterior rows remain available. New sampler runs replace these
stored summaries with their normal runtime sampler objects.
