"""Offline guides for the application's help and READ ME buttons."""

HELP_PAGES = {
    'magnetic_cycle': ('Magnetic cycle — help', """
<h1>Magnetic cycle (shared sine model)</h1>
<p>Open <b>GP param. → Joint act.+ RV Magnetic cycle proxy (sine)</b>, after the
joint activity GP tab. These tabs can be dragged into another order.</p>
<p>This is a deterministic Exo-Striker model, outside celerite, with no planet
slot or transit GP component. It can be used alone or alongside any existing RV
GP. Enable it separately from <b>Use RV GP</b>.</p>
<pre>C(t) = sin[2π(t − tref)/P + phase]
RV(t) = orbital model + RV amplitude × C(t) + optional rotational GP
Index_j(t) = offset_j + amplitude_j × sin[2π(t − tref)/P + RV phase + phase offset_j]
             + optional GP response + noise</pre>
<p><b>P [d]</b> is shared across RVs and selected activity indices.
<b>RV phase [deg]</b> defines the RV cycle. Each index has its own
<b>Phase offset [deg]</b> relative to the RV phase; check its <b>Fit</b> box to
infer an independent phase. Zero offsets with Fit unchecked preserve the previous
shared-phase model. Offsets default to zero/fixed with bounds −180 to +180 degrees. <b>Phase reference [BJD]</b> is fixed, initially the first RV time; phase
is measured at that reference. Changing it changes the phase convention, so
refit/adjust phase accordingly. The reference is independent of the orbital epoch.</p>
<p><b>RV amplitude [m/s]</b> is a semi-amplitude (half the peak-to-peak range),
nonnegative by default. Each selected index has a signed <b>amplitude [U]</b>,
where U is its uploaded unit. Negative values allow anticorrelation. Its offset
and jitter also have units U. Measurement noise variance is error² + jitter².</p>
<h2>Independent selections</h2>
<p>Select a loaded file, check <b>Include in magnetic cycle</b>, then repeat for
other indices. The summary lists the cycle selection. This is separate from
<b>Joint activity GP</b>. For example, select FWHM here, and select BIS and Hα in
Joint activity GP with dSHO enabled in RV GP. Both processes contribute to RVs.</p>
<p>An index can participate in either or both components. If both are selected,
its cycle mean is subtracted before joint GP conditioning; its offset and jitter
are shared, and its observations contribute once to the likelihood. The offset
and jitter edits on the two forms refer to the same parameters.</p>
<h2>Fitting and priors</h2>
<p>Rows place Fit before the value and show an uncertainty label beside it
(three decimals, +/- the larger absolute lower/upper error, as in the planet fields; tiny nonzero errors use scientific notation). Activity
fields retain small values such as RHK-prime amplitudes, offsets, jitters and priors.
These display stored inference
errors from MCMC, nested sampling or SciPy when available. Rows supply bounds, prior family and Gaussian mean/sigma or Jeffreys
endpoints. These use the row's units. Period bounds must be positive; jitter
bounds nonnegative. Signed activity amplitudes normally need uniform bounds or
Gaussian priors. Jeffreys support must be strictly positive. Only checked Fit
boxes add dimensions. The fixed reference epoch is not a fitted parameter.</p>
<p>The initial period is max(3650 d, twice the RV time span); amplitudes are
scatter-based guesses. Set a plausible cycle period and sensible bounds before
fitting. Defaults leave the component disabled, with period, phase and amplitudes
free when enabled; activity offsets start at their medians and jitter at zero/fixed.</p>
<p>Initialize, Fit, MCMC and nested sampling evaluate the cycle jointly with the
orbital model and selected GP. Active cycles use the Python likelihood/SciPy fit
path, including when the main Fit button normally selects Fortran. No Fortran
optimizer change is required. Priors enter inference normally; reported lnL
is the data likelihood, excluding prior density. Statistics count the union of
cycle and GP activity observations, not two copies of overlapping series.</p>
<p>The RV model curve includes the sine contribution; RV residuals subtract it.
Activity time-series plots show the cycle or cycle-plus-GP mean after model
initialization/fitting. Save the session to retain values, fit flags, bounds,
priors, selections and uncertainties. Removing an index removes its selection.</p>
<h2>Limitations</h2>
<p>The model assumes one coherent sinusoid per series, sharing only its period
when the activity phase offsets are free. It does not model evolving cycle shapes.
A negative activity amplitude is equivalent to a positive amplitude and a 180-degree
phase change: use nonnegative activity-amplitude bounds when fitting a full-circle
phase offset to avoid duplicate solutions. Phases are periodic; −180 and +180
degrees are equivalent. An offset is a relative phase, not a unique physical time delay. Less than a full observed
cycle gives a poorly constrained period and can be degenerate with RV trends
or long-period planets. Compare alternative models and inspect posterior
correlations. Cycle-only indices use independent Gaussian errors; selected GP
indices use the existing covariance. The dense joint GP cost is unchanged.</p>
"""),
    'joint_activity_gp': ('Joint RV and activity Gaussian process — help', """
<h1>Joint RV and activity Gaussian process</h1>
<p>Fit RVs and several activity indicators (BIS, FWHM, Hα, NaD1, dLW, etc.)
simultaneously with one shared process G(t). <b>Time remains the regressor.</b>
Activity measurements constrain the same process as the RV residuals through
cross-correlations. This is not an activity-only fit subsequently copied into RV.</p>
<h2>1. Select the observations</h2>
<ol><li>Load RV data and activity files. Each activity file supplies time, value
and measurement uncertainty. Use finite values, positive errors and the same
time system. Unequal timestamps are supported: no interpolation or matching is required.</li>
<li>Under <b>GP param. → RV GP</b>, enable <b>Use RV GP</b> and choose the kernel,
its values and Fit flags. Existing kernel bounds and priors remain in use.</li>
<li>Open <b>Joint activity GP</b> (after Transit GP by default) or click
<b>Select activity indicators…</b>. The selector lists loaded files only.</li>
<li>Select a file, check <b>Include in joint RV GP</b>, then repeat for other
indicators. The summary lists all included series. Changing the selector only
changes which series' parameters are displayed; it does not deselect others.</li>
<li>Initialize to evaluate the model, then fit or sample. Inspect predictions,
residuals and posterior correlations. Drag the GP tabs to change their order.</li></ol>
<p>No included indicators uses the existing RV-only time GP. Switching off
Use RV GP disables joint conditioning but retains selections. <b>Include</b>
controls dataset participation; <b>Fit</b> controls whether a parameter varies.
Values and choices belong to the saved fit/session, not global plot settings.</p>
<h2>2. Shared process and response parameters</h2>
<pre>RV activity(t) = G(t) + b_RV dG/dt
Index_j(t)     = offset_j + a_j G(t) + b_j dG/dt + noise_j
noise variance = measurement error² + activity jitter²</pre>
<p>The existing kernel parameters describe G's amplitude and time dependence
(period, coherence, decay, etc.). The responses convert this process into each
indicator's units and allow different signs and shapes. Identical numerical
curves cannot generally describe observables with different units or baselines.</p>
<table cellspacing="8" cellpadding="4">
<tr><th align="left">Row</th><th align="left">Meaning and initial setting</th></tr>
<tr><td>Activity × G</td><td>Signed direct response a_j. Initially estimated from
scatter; Fit on. A negative value allows anticorrelation. For indicator units U,
its units are U/(m/s). The guess is not a physical calibration.</td></tr>
<tr><td>Activity × dG/dt</td><td>Optional derivative response b_j: sensitivity to
how quickly activity changes, allowing a different shape or apparent phase.
Units U day/(m/s). Initially zero with Fit off.</td></tr>
<tr><td>Activity offset</td><td>Individual baseline, initially the median of the
series with Fit on. G itself has zero mean. Units U.</td></tr>
<tr><td>Activity jitter</td><td>Additional independent scatter, added in quadrature
to measurement errors. Units U. Initially zero with Fit off.</td></tr>
<tr><td>RV × dG/dt [d]</td><td>One common RV derivative response for all RV
instruments, in days. It is the same parameter whichever activity file is
selected. Initially zero with Fit off.</td></tr></table>
<p><b>RV × G is fixed to one.</b> This anchors the amplitude to the existing RV
kernel convention and avoids redundant scaling between its amplitude and all
responses. Existing instrument-specific RV offsets and jitters are separate.
Kernel amplitude conventions differ: SHO S0 is a spectral-density scale;
Rotation/DRW amplitudes are covariance scales; sigma is a standard deviation.</p>
<h2>3. Which kernels allow derivatives?</h2>
<p><b>All five existing RV kernels support joint direct responses:</b> Rotation,
DRW/RealTerm, SHO, double SHO and the approximate Matérn 3/2.</p>
<p>Derivative responses require <b>SHO, double SHO or Matérn 3/2</b>. Rotation and
DRW contain exponential cusps, so the required process derivative is not defined
in the mean-square sense. Their derivative controls are grayed out and the
model uses zero fixed derivatives. Previous derivative values, Fit flags,
bounds and priors are retained and restored when you return to a supported kernel.</p>
<h2>4. Bounds, priors and Fit flags</h2>
<p>Each row supplies a value, Fit checkbox, lower/upper bounds and prior family.
The last two fields are Gaussian mean/standard deviation or Jeffreys lower/upper
support. Signed coefficients usually use uniform bounds or Gaussian priors;
Jeffreys support must be strictly positive. Review all starting guesses and
bounds. Fixing a parameter retains its value's contribution; it does not remove it.</p>
<p>Start with direct responses and offsets. Add derivatives only when justified
by the model and a compatible kernel, and jitter when the reported errors do not
explain independent scatter. Independently known responses or offsets can be fixed.</p>
<h2>5. How many extra free parameters?</h2>
<p>Only checked Fit boxes add dimensions. <b>Three indicators with direct response
and offset fitted add six parameters</b> (two per series), in addition to the
existing kernel, orbital and instrument parameters. Fixed derivative/jitter
rows add no free parameters. With M indicators, fitting all four per-series rows
plus the one shared RV derivative adds <b>4M + 1</b> parameters (13 for three).</p>
<h2>6. Inference, likelihood and plots</h2>
<p>The joint likelihood runs in Python with SciPy optimization, MCMC and nested
sampling; it requires no additional Fortran GP implementation. With enabled
priors, SciPy finds a MAP solution; MCMC samples the posterior; nested sampling
uses the data likelihood and its prior transform. <b>Reported lnL excludes prior
density</b> and includes the selected activity observations and other selected
components. The RV statistics row represents the joint RV–activity block,
including activity counts and covariance-weighted chi-square, rather than an
independent RV-only likelihood.</p>
<p>Activity uploads populate time-series and GLS views automatically. After
initialization/fitting, included indicators show their smooth conditional GP
mean in red; RV GP predictions also condition on all selected indicators.
These are predictions at the current parameters, not posterior-averaged means
or uncertainty bands. The activity GLS shows measurements, not GP-subtracted
residuals. Unselected or uninitialized indicators have no joint GP curve.</p>
<h2>7. Why can it be more expensive?</h2>
<p>There are two costs. More free parameters usually require more evaluations
and may introduce posterior degeneracies. Separately, the current joint solver
constructs a <b>dense covariance</b> from the existing celerite kernel coefficients:
N = RV observations + all included activity observations. Measurements of different
observables at the same epoch still count separately. Memory grows as O(N²)
and factorization as O(N³), even if every added response parameter is fixed.
Doubling N means roughly four times the covariance memory and eight times the
factorization work. This joint implementation does not use celerite's fast solver;
the fast RV-only path remains available when no activity series is included.</p>
<p>A common process is a modeling assumption. Check sensitivity to kernel choice
and priors and inspect correlations with planetary signals. A synthetic recovery
alone does not demonstrate that real stellar activity and planets are separated.</p>
"""),
    'ttv': ('Transit timing variations — help', """
<h1>Transit timing variations</h1>
<p>Fit measured mid-transit times with the dynamical TTVFast model, either alone
or together with radial velocities. Each timing dataset must be assigned to
an included planet and enabled for the fit.</p>
<h2>1. Prepare a timing file</h2>
<p>Use whitespace-separated numeric columns. Comment lines start with #.
The first three columns are transit number, mid-transit time in BJD, and its
uncertainty in days. Supply finite values and positive uncertainties.</p>
<pre># N   t0 [BJD]       sigma_t0 [days]
1     2458000.5      0.0022
2     2458020.5      0.0023
4     2458060.5      0.0021
8     2458140.5      0.0026</pre>
<p>Transit numbers start at 1 for the first modeled transit of that planet.
Keep missing transits as gaps in N; do not renumber the example's 4 and 8 as 3
and 4. Use full Julian dates consistently across the model and input files,
with a consistent time standard. Convert uncertainties in minutes to days by
dividing by 1440.</p>
<h2>2. Set the epoch and integration span</h2>
<p>Place the epoch slightly before the first transit that should be numbered 1,
and set <b>End of Model</b> later than the last observed transit. For this example,
an epoch of 2458000.0 precedes the first timing. Starting too late can shift
the model's transit numbering or omit an observation.</p>
<p>For a joint RV + TTV fit, the RV model epoch is used. Open
<b>Models param. → Models → RV Model</b>, uncheck <b>first RV</b>, and enter a
suitable epoch. Check the integration end time as well as the epoch.</p>
<h2>3. Check the model before sampling</h2>
<ol><li>Assign each timing file to the correct planet and enable that dataset.</li>
<li>Include the planets that perturb the transiting planet, even if they have
no timing dataset of their own. Check stellar mass and orbital geometry.</li>
<li>Run a fit and inspect the calculated times and O − C residuals.</li>
<li>Reduce the <b>time step in dynamical model</b> and compare the timings.
The numerical change should be small relative to the timing uncertainties.</li></ol>
<h2>Optional depth and duration columns</h2>
<p>The importer also accepts columns 4–7: depth, depth uncertainty, duration
in days, and duration uncertainty in days. Depth is a fractional flux loss
(e.g. 0.01 for 1%). These columns add depth/duration constraints, so omit them
for a timing-only file. Supply the stellar density and planet-to-star radius
ratio required by those constraints. The depth model uses a uniform stellar
disk; it is not a replacement for a limb-darkened light-curve fit.</p>
<h2>If the model does not cover the observations</h2>
<p>Check the planet assignment, dataset enable checkbox, transit numbering,
epoch, end time, and time units first. A missing calculated transit can also
indicate a non-transiting geometry. For assistance, include a small input
example and the relevant model settings in an
<a href="https://github.com/3fon3fonov/exostriker/issues">Exo-Striker issue</a>.</p>
"""),
    'gp': ('Gaussian process modeling — help', """
<h1>Gaussian process modeling</h1>
<p>A Gaussian process (GP) models time-correlated residuals alongside the
orbital or transit signal. The RV and transit GP controls are separate;
select a kernel, enable the GP, and choose which hyperparameters to fit.</p>
<h2>Parameter conventions</h2>
<p>Enter the physical parameter values displayed in the GUI. The backend
converts positive parameters to logarithms where required by celerite;
do not enter their logarithms yourself. Periods and time scales use days.
RV amplitudes use the RV data units, while transit amplitudes use the flux
units supplied to the fit. A covariance amplitude has squared data units.</p>
<p>Labels include units in the value, bounds and prior panels. <b>[1]</b> means
 dimensionless; <b>flux</b> means the units of the supplied photometric values
 (relative flux is dimensionless). With signal units Y = m/s for RV or flux for
 photometry, the conventions are:</p>
<table cellspacing="6" cellpadding="3">
<tr><th align="left">Parameter</th><th align="left">Units</th></tr>
<tr><td>SHO S₀</td><td>Y² d; variance = S₀ ω₀ Q</td></tr>
<tr><td>SHO ω₀</td><td>rad/d</td></tr>
<tr><td>SHO Q; double-SHO Q₀, ΔQ, f; Rotation factor</td><td>1 (dimensionless)</td></tr>
<tr><td>Double-SHO σ; Matérn σ</td><td>Y</td></tr>
<tr><td>Double-SHO P; Rotation period/decay time; Matérn ρ</td><td>d</td></tr>
<tr><td>Rotation amplitude; DRW a</td><td>Y²</td></tr>
<tr><td>DRW c; Matérn eps</td><td>d⁻¹</td></tr></table>
<p>In the implemented celerite approximation, eps multiplies the time lag in
 an oscillatory term, so its unit is inverse days. It controls approximation
 quality, rather than an additional physical variability amplitude. Bounds,
 Gaussian means/standard deviations and Jeffreys endpoints use the associated
 parameter's units.</p>
<h2>Available kernels</h2>
<table cellspacing="8" cellpadding="4">
<tr><th align="left">Kernel</th><th align="left">Parameters and interpretation</th></tr>
<tr><td>Rotation</td><td><b>Amp.</b>: covariance amplitude; <b>timescale</b>:
exponential decay time; <b>Period</b>: repeating time scale; <b>fact.</b>:
relative weight of the non-oscillating part. All four are positive.</td></tr>
<tr><td>SHO</td><td><b>S</b> (S0): spectral-density scale; <b>Q</b>: quality
factor; <b>ω</b>: angular frequency in radians/day. All are positive.
2π/ω is the undamped characteristic period; S is not the process standard
deviation.</td></tr>
<tr><td>Matérn 3/2</td><td><b>σ</b>: process standard deviation; <b>ρ</b>:
correlation time; <b>eps</b>: approximation parameter. Keep σ and ρ positive.
eps controls the celerite approximation, rather than a physical noise amplitude.</td></tr>
<tr><td>Double SHO</td><td><b>σ</b>: process standard deviation; <b>P</b>:
primary period; <b>Q0</b> and <b>dQ</b>: quality-factor controls; <b>f</b>:
relative power of the secondary mode. The model combines variability near P
and P/2. Use positive σ, P, Q0, and dQ and a positive f (usually below 1).
The additional dummy slot is not a kernel parameter.</td></tr>
<tr><td>DRW / RealTerm</td><td><b>a</b>: covariance amplitude; <b>c</b>:
inverse decay time. The covariance is a exp(−c |Δt|), so the correlation
time is 1/c. Both parameters are positive.</td></tr></table>
<h2>A practical workflow</h2>
<ol><li>Fit the deterministic model first and inspect its residuals.</li>
<li>Choose a kernel that can describe the variability you see. Set starting
values and bounds consistent with the data units and observing baseline.</li>
<li>Enable fitting only for the hyperparameters you intend to vary. Avoid
zero or negative bounds for parameters that are internally logged.</li>
<li>Refit, inspect the GP prediction and residuals, and then explore
uncertainties with MCMC or nested sampling.</li></ol>
<p>White-noise jitter and a correlated GP describe different components.
Inspect their posterior correlations with the planetary parameters. If you
link transit GP parameters to the RV GP, check which values are shared;
flux and RV amplitude scales are different.</p>
<h2>Further reading and citation</h2>
<p><a href="https://celerite.readthedocs.io/en/stable/python/kernel/">celerite kernel documentation</a>
and <a href="https://ui.adsabs.harvard.edu/abs/2017AJ....154..220F/abstract">
Foreman-Mackey et al. (2017)</a> describe the GP backend. Cite the backend
and the relevant kernel references when using it in a publication.</p>
"""),
    'mcmc': ('MCMC — help', """
<h1>MCMC posterior sampling</h1>
<p>Exo-Striker uses emcee to explore the posterior distribution of the selected
free parameters. Check the model, dataset selection, parameter bounds, and
priors before starting a run.</p>
<h2>Walkers, burn-in, and production</h2>
<p>Walkers are parallel chains in the ensemble. For the usual stretch move,
use at least twice as many walkers as free dimensions. Initialize them at
distinct valid positions. Burn-in lets the ensemble settle; exclude that
stage from the posterior summaries. Production steps are the samples used
for inference.</p>
<h2>Check convergence</h2>
<p>Inspect traces for drift and poor mixing, estimate autocorrelation times,
and compare independent runs. An acceptance fraction alone does not establish
convergence. Run long enough relative to the autocorrelation time before
reporting uncertainties. Thinning reduces stored samples; it does not repair
a chain that has not converged.</p>
<h2>Read the results</h2>
<p>Use a cornerplot to inspect parameter correlations and multiple modes.
The median, mean, and maximum-likelihood sample can differ. A narrow-looking
histogram does not establish that all posterior modes were explored.</p>
<p><a href="https://emcee.readthedocs.io/en/stable/tutorials/quickstart/">emcee quickstart</a>
· <a href="https://emcee.readthedocs.io/en/stable/tutorials/autocorr/">Convergence guide</a>
· <a href="https://ui.adsabs.harvard.edu/abs/2013PASP..125..306F/abstract">Foreman-Mackey et al. (2013)</a>.</p>
"""),
    'ns': ('Nested sampling — help', """
<h1>Nested sampling</h1>
<p>Exo-Striker uses dynesty to estimate Bayesian evidence and obtain posterior
samples. Define the model and priors before running; evidence depends on the
prior ranges as well as on the fit.</p>
<h2>Main controls</h2>
<p><b>Live points</b> control how well the evolving distribution is represented.
More live points generally increase cost and improve exploration.
<b>Static</b> sampling uses a fixed live-point allocation; <b>dynamic</b>
sampling allocates additional effort during the run. Bounding and sampling
options affect how new constrained points are proposed.</p>
<p>The evidence stopping tolerance (dlogz) controls the remaining estimated
evidence contribution. Smaller tolerances require more work. Iteration or
likelihood-call limits can stop a run before that tolerance is reached.</p>
<h2>Interpret the output</h2>
<p>Check the termination reason, evidence uncertainty, and posterior stability
across independent runs. Compare log-evidence only for models fitted to the
same observations with consistent likelihood conventions. Document the priors;
widening them can change the evidence. Use weighted samples or the application's
resampled posterior for summaries, rather than treating all raw nested samples
as equally probable.</p>
<p><a href="https://dynesty.readthedocs.io/en/latest/quickstart.html">dynesty user guide</a>
· <a href="https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.3132S/abstract">Speagle (2020)</a>.</p>
"""),
    'detrend': ('Transit detrending — help', """
<h1>Transit detrending</h1>
<p>Estimate a baseline trend and inspect the flattened light curve before using
it for transit fitting. The upper plot shows the data and estimated trend;
the lower plot lets you inspect the resulting data product.</p>
<h2>Try a method</h2>
<ol><li>Choose an offset, sliding-window, polynomial, spline, regression, or
GP method from the controls available in your installation.</li>
<li>Set the window length in the same units as the input times (normally days).
Break tolerance separates segments across observing gaps. GP controls instead
specify a kernel scale and, where applicable, a period.</li>
<li>Press <b>Try!</b>, wait for the calculation, and inspect both plots.
Choose a trend scale that preserves the transit shape. Compare different
settings rather than relying on the smallest residual scatter.</li>
<li>After rejecting or binning points, press <b>Try!</b> again before closing.
Use <b>Reset</b> to reload the window's input and <b>Save detr. data</b> to
export the resulting product.</li></ol>
<h2>Mask known transits</h2>
<p>Select a compatible method, check <b>Use transit mask</b>, and enter the
<b>Period [d]</b>, <b>T0 [BJD]</b>, and <b>Mask width [d]</b>. T0 is a mid-transit
epoch in the same time system as the plotted data. The width is the full
excluded interval, not a half-width: each mask spans T0 + n × Period ± width/2.
Use a width that covers ingress through egress, with a little margin; it must
be positive and shorter than the period. Leave enough unmasked data around
each transit to constrain the baseline.</p>
<p>The upper plot shows solid amber lines at transit centers and shaded bands
across the excluded intervals. These decorations update immediately as you
change the settings; press <b>Try!</b> to recalculate the trend and flattened
light curve using those settings.</p>
<p><b>Why does the trend continue through a shaded interval?</b> Masking excludes
those flux measurements from fitting the baseline; it does not remove them
from the light curve. The baseline fitted from unmasked data is still evaluated
inside the mask. The flattened flux is the original flux divided by that
baseline, including the transit points. A continuous black trend through a
shaded band is therefore expected and does not mean the transit was used in
the fit. Wide masks, data gaps, or windows containing too few unmasked points
can leave the baseline poorly constrained or produce undefined values.</p>
<p><b>Method compatibility:</b> the checkbox is disabled and unchecked for
cofiam, savgol, medfilt, pspline, supersmoother, ridge, lasso, and elasticnet,
which ignore the transit mask in the released Wotan 1.10 backend. The mask
fields are disabled and the shading is hidden when masking is unchecked. Mask-aware
choices include time-windowed estimators such as biweight, cosine, lowess,
rspline, hspline, and Gaussian processes. Mean/median normalization also
excludes masked points when calculating its baseline. Wotan documents cosine
and lowess support explicitly and treats most other methods as experimental;
compare results and inspect the transit shape. Wotan 1.10 already provides
masking for these compatible methods; a development version is not required.
The GUI uses this conservative compatibility list even with newer backends.</p>
<h2>Optional algorithms</h2>
<p>Some methods require extra packages in the same Python environment as
Exo-Striker: statsmodels for certain robust estimators; scikit-learn for GP,
Huber splines, ridge, and lasso; pygam for penalized splines; supersmoother
for its corresponding method. Only available methods appear in the selectors.</p>
<pre>python -m pip install statsmodels scikit-learn pygam supersmoother</pre>
<p><a href="https://wotan.readthedocs.io/en/latest/Interface.html">wotan interface guide</a>
· <a href="https://ui.adsabs.harvard.edu/abs/2019AJ....158..143H/abstract">Hippke et al. (2019)</a>.
Cite wotan when using these detrending algorithms in a publication.</p>
"""),
    'astrometry': ('Imaging astrometry — help', """
<h1>Imaging astrometry</h1>
<p>Use the imaging astrometry controls for relative measurements of a companion.
Assign each dataset to the appropriate included planet and select the
astrometry data type before fitting.</p>
<h2>Input file</h2>
<p>Use a whitespace-separated numeric file with five columns; comment lines
start with #. The imaging importer reads:</p>
<pre># BJD        x         sigma_x   y        sigma_y
2458000.5    120.0     1.0       45.0     1.0</pre>
<p>Use full Julian dates, finite coordinates, and positive uncertainties.
Match x and y to the selected astrometry convention and use the angular units
expected by that model; the importer does not infer units from the header.</p>
<h2>Check the input convention</h2>
<p>Verify the time, angular units, coordinate convention, and uncertainty
columns expected by the selected importer. Position-angle and separation data
must be interpreted consistently with the orbital orientation. Check the
ascending node and inclination against a plotted orbit before sampling.</p>
<h2>Joint fits</h2>
<p>For an RV + astrometry fit, keep the reference epoch, stellar mass, and orbital
parameters consistent across both datasets. Inspect residuals for each dataset
separately. Hipparcos and Gaia along-scan measurements use a different input
format; their dedicated READ ME button describes the required catalogue header
and scan geometry.</p>
<p>If the input convention is unclear, consult the bundled example datasets or
include a small example in an <a href="https://github.com/3fon3fonov/exostriker/issues">
Exo-Striker issue</a>.</p>
"""),
}


def data_options_help(kind):
    """Help shared by the RV and activity preprocessing windows."""
    label = 'Radial velocity' if kind == 'rv' else 'Activity'
    return (f'{label} data options — help', f"""
<h1>{label} data options</h1>
<p>Inspect and preprocess the selected dataset. The upper plot shows the
input data and trend; the lower plot shows the processed values.</p>
<h2>Offsets and trends</h2>
<p><b>No removal</b> keeps the values; <b>remove mean</b> and <b>remove median</b>
subtract a constant. Other enabled methods estimate a trend. The processed
series is constructed by subtracting the estimate returned by the selected
method; inspect it before using it in a fit.</p>
<p>Spline, polynomial, regression, and GP methods use their displayed window
or kernel settings. Time-window values follow the input time units.
The lowpass, highpass, bandpass, and bandstop controls use the displayed
order and reciprocal cutoff values; check their response carefully on
irregularly sampled observations.</p>
<h2>Preview and inspect</h2>
<ol><li>Choose a method and press <b>Try!</b>.</li>
<li>Inspect the trend, processed series, and statistics. Preprocessing can
alter signals you later want to fit.</li>
<li>Use point rejection, sigma clipping, or binning only after inspecting the
data. Recalculate with <b>Try!</b> after changing the points.</li>
<li>Use <b>Reset</b> to reload the window's input. <b>Save modif. data</b>
exports the processed product; inspect that file before reloading it.</li></ol>
<h2>Uncertainties and interpretation</h2>
<p>Trend subtraction retains the supplied measurement uncertainties; those
errors do not include uncertainty in the estimated trend. For RV analysis,
consider fitting offsets, trends, or a GP jointly with the planetary signal.
For activity data, keep the indicator's units and compare its variability
with RV residuals rather than assuming the same amplitude scale.</p>
<p><a href="https://wotan.readthedocs.io/en/latest/Interface.html">wotan method reference</a>
· <a href="https://github.com/3fon3fonov/exostriker/issues">Exo-Striker support</a>.</p>
""")


HELP_PAGES['rv'] = data_options_help('rv')
HELP_PAGES['activity'] = data_options_help('activity')


HELP_PAGES['rvbank'] = ('RVBank — help', """
<h1>Browsing RVBank data</h1>
<p>Choose <b>HARPS RVBank</b> or <b>HIRES NZP</b>, filter the target list by
name, and select a target and data product. Browsing retrieves data from the
remote archive and requires an internet connection.</p>
<h2>Choose the data product</h2>
<p>Check whether you selected velocities or an activity indicator before
loading the data. Keep the product's units and source corrections with the
exported file. When combining instruments, fit the appropriate dataset
offsets and jitter terms.</p>
<h2>Activity uncertainties and missing measurements</h2>
<p>Some activity columns do not have measured uncertainties. The importer
assigns placeholder errors in those cases; they are not formal measurement
errors. Inspect them before a weighted analysis. Some HIRES h-index values
are −1, indicating that the measurement could not be computed; exclude
those missing values rather than interpreting them as observations.</p>
<h2>Sources and citations</h2>
<ul><li>HARPS RVBank: <a href="https://ui.adsabs.harvard.edu/abs/2020arXiv200105942T/abstract">
Trifonov et al. (2020)</a> and the data references listed there.</li>
<li>Original HIRES velocities: <a href="https://ui.adsabs.harvard.edu/abs/2017AJ....153..208B/abstract">
Butler et al. (2017)</a>.</li>
<li>NZP-corrected HIRES velocities: <a href="https://ui.adsabs.harvard.edu/abs/2019MNRAS.484L...8T/abstract">
Tal-Or et al. (2019)</a>, together with Butler et al. (2017).</li></ul>
<p>If browsing fails, check the internet connection and try another target.
Include the target name, archive choice, and any error message when reporting
an <a href="https://github.com/3fon3fonov/exostriker/issues">Exo-Striker issue</a>.</p>
""")
