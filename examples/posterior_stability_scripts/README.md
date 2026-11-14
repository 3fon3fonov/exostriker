# Integrating posterior samples

Provided by **Trifon Trifonov with the aid of CODEX; use with caution**.

`integrate_posterior.py` replaces the system-specific workflow in `make_stab.py`
(the latter is retained as a historical reference). `HD128311.ses` contains three
active planets and nested-sampling posteriors. The notebook demonstrates the same
functions, joint θ/NV2016 posteriors, and the Ψ cos ψ versus δ structure diagram,
without a LaTeX installation. Its default 2:1 pair is GUI slots 2 and 3;
for TIC279401253 select its session and use slots 1 and 2.

## Run with your own `.ses` file

1. Activate the Python environment you use to run Exo-Striker. The commands below
   use `python`; if your working interpreter is `python3.10`, use that instead.
2. Save your session **after finishing MCMC or nested sampling**, with samples
   stored in the session. A best-fit-only session is insufficient. Use the original
   sampled session, before changing its planets or parameter configuration.
3. Install an Exo-Striker version containing `exostriker.lib.RV_mod.dynamical_posterior`
   in that environment. Copy `integrate_posterior.py` anywhere you like; it uses
   the installed package when copied outside the checkout, and needs no
   `example_tools.py`. Inside this repository it automatically uses the checkout. Open a
   terminal in the directory containing your copy.

4. Run the command below, replacing `/path/to/MySystem.ses` with your session's
   full path and `/path/to/MySystem_stability_test` with a **new directory** to
   create. Quote paths, especially if they contain spaces. Do not create the output
   directory beforehand.

   ```bash
   python integrate_posterior.py \
     "/path/to/MySystem.ses" \
     --output "/path/to/MySystem_stability_test" \
     --source ns \
     --count 3 \
     --years 10 \
     --timestep-days 1 \
     --output-years 1
   ```

   This short first run checks that your session can be loaded and integrated.
   `--source ns` selects nested sampling. **For MCMC use `--source mcmc` instead.**
   You can omit `--source` if the session contains only one sampler.

5. After completion, open this file through Exo-Striker's normal session-open
   dialog:

   ```text
   /path/to/MySystem_stability_test/MySystem_integrated.ses
   ```

   The output filename uses your input filename without `.ses`, followed by
   `_integrated.ses`. Inspect **Sample corr.** and its corner plot for the added
   integration columns. CSV tables and diagnostics are in the same directory.
   Your input `.ses` is unchanged.

## REBOUND integrators

Use `--integrator whfast` (or `mercurius`, `saba`, `ias15`, `trace`, `janus`)
to select REBOUND automatically. The default remains SWIFT `symba`; `mvs` and
`mvs_gr` are also available. Install REBOUND in the same Python environment with
`python -m pip install rebound` if needed.

```bash
python integrate_posterior.py "MySystem.ses" \
  --output "MySystem_rebound" --count 3 --years 10 --integrator ias15
```

The notebook/API accepts the same names, e.g. `IntegrationOptions(integrator='whfast')`.
Timesteps remain in days, durations/output intervals in years. IAS15 adapts its
step; `--timestep-days` supplies its initial step. Fixed-step integrators may
slightly overshoot output times; diagnostics use the actual times. Both backends
use Jacobi orbital histories and the same stability thresholds and exports.
The timeout is checked every REBOUND step. These REBOUND runs are Newtonian;
use SWIFT `mvs_gr` for the existing relativistic integration option.

## Older sessions and checkout versions

Older sampled sessions without `stellar_mass_str`, `mcmc_planet_slots`, or
`mcmc_dataset_slots` are supported when their saved free-parameter indices and
orbital labels still agree. Fixed stellar mass and fixed orbital elements are
retained. The script does not rebuild sampling metadata or modify the input file;
sessions with changed or ambiguous orbital mappings are rejected.

When run from this repository, the script automatically uses the matching
checkout, including in spawned integration workers. It prints the package path
before loading the session, so you can verify which version is running. No
`PYTHONPATH` setting is needed.

A script copied outside the checkout uses the package installed in its Python
environment; update that installation to obtain compatibility fixes.

## Choose how many posterior samples to randomly draw

Set **`--count N`** to the number of random posterior draws you want to integrate:

```bash
python integrate_posterior.py \
  "/path/to/MySystem.ses" \
  --output "/path/to/MySystem_random100" \
  --source ns --count 100 --seed 42 \
  --years 1000 --timestep-days 1 --output-years 1
```

Here `--count 100` randomly draws **100 samples**. Change it to `--count 1000`
for **1,000 draws**, or any positive integer. Use `--source mcmc` for MCMC.
The short test above uses only `--count 3`.

Draws follow the posterior probabilities: uniform for equal-weight MCMC samples,
and weighted for raw nested-sampling samples. Selection is **with replacement**,
so the same saved row can be drawn more than once; `N` is the number of draws,
not necessarily the number of distinct integrations. Repeated rows reuse their
integration result. Keep `--seed 42` to reproduce the selection; change the seed
for another random selection.

To integrate **all saved rows without random selection**, omit `--count` entirely.
Do not pass `--count 0`. Start with a small count to estimate runtime before
choosing the count for your full analysis.

## Longer run with your own system

After the short test works, choose settings appropriate to your system. For
example, this command tests 1,000 posterior draws for 100,000 years each:

```bash
python integrate_posterior.py \
  "/path/to/MySystem.ses" \
  --output "/path/to/MySystem_stability_long" \
  --source ns \
  --count 1000 \
  --years 100000 \
  --timestep-days 1 \
  --output-years 10 \
  --timeout-seconds 1500 \
  --integrator symba \
  --max-fractional-a-change 0.1 \
  --max-eccentricity 0.95 \
  --seed 42
```

These settings are illustrative and can take a long time. Change `--source ns` to
`--source mcmc` for MCMC. **To integrate every saved row, remove `--count 1000`
entirely**; do not use `--count 0`. A session can contain many thousands of rows.
Each repeated run needs a different output directory.

| Option | Meaning |
| --- | --- |
| `--count 1000` | Draw 1,000 samples from the posterior, with replacement; omit for every saved row. |
| `--years 100000` | Requested integration duration for each draw, in years. |
| `--timestep-days 1` | Numerical integration step, in days; choose well below the shortest period. |
| `--output-years 10` | Interval between recorded states and stability checks, in years. |
| `--timeout-seconds 1500` | Maximum allowed integration time per draw, in seconds. A timeout is recorded separately. |
| `--integrator symba` | SWIFT or REBOUND integrator; see above. |
| `--max-fractional-a-change 0.1` | Reject an orbit if its recorded semi-major axis changes by more than 10% from its initial value. |
| `--max-eccentricity 0.95` | Reject an orbit if its recorded eccentricity exceeds 0.95. |
| `--seed 42` | Reproducible posterior selection and GUI resampling. |

The script prints progress such as `1/1000: stable`. It saves
its output files after all selected draws finish. It does not currently checkpoint
or resume interrupted work. For all command-line options:

```bash
python integrate_posterior.py --help
```

## Parallel CPU integration

Use **`--workers N`** to run up to N independent integrations simultaneously.
The default is `--workers 1`. From **inside `examples/stability`**, for example:

```bash
python3.10 integrate_posterior.py HD128311.ses \
  --output HD128311_stab_parallel \
  --source ns --count 1000 --workers 8 \
  --years 1000 --timestep-days 0.5 --output-years 1
```

This uses up to **8 CPU processes**. Choose a number suitable for your CPU and
available memory; each worker loads Exo-Striker's Python dependencies. SWIFT uses
shared native state inside a process, so independent processes provide parallel
execution safely. Python threads would serialize the integrations.

The selected samples, their order and weights remain reproducible with the same
seed regardless of worker count. Progress reports completed draws in completion
order; repeated draws share one integration. Increasing workers can improve
throughput for long runs, but startup overhead can dominate short tests. The
per-draw timeout applies to its integration, not time waiting in the process queue.

For Python scripts, protect the call with a main guard for spawned processes:

```python
if __name__ == '__main__':
    result = integrate_posterior(fit, options, source='ns', count=1000, workers=8)
```

The provided command-line script already includes this guard. The notebook can
pass `workers=8` directly to the imported `integrate_posterior` function.

## Included HD128311 example

From the repository root, using your Exo-Striker Python environment:

```bash
python integrate_posterior.py examples/stability/HD128311.ses \
  --output /tmp/HD128311-smoke --source ns --count 3 --years 10 \
  --timestep-days 1 --output-years 1
```

Open `/tmp/HD128311-smoke/HD128311_integrated.ses` when it finishes. If that output
directory already exists, choose a different name.

This is a smoke test, **not a scientific stability determination**. For research,
choose an appropriate duration, integration step (well below the shortest orbital
period), output cadence, and thresholds. Check convergence by reducing the step
and output interval and increasing the duration. Defaults are examples, not a
universal prescription. All enabled planets are integrated, including those
marked Keplerian in a mixed RV fit. Initial elements are interpreted as Jacobi
osculating elements at the session epoch. Systems without RV amplitudes sufficient
to determine physical masses (for example transit-only fits) need an explicit
mass model and are not supported by this RV-to-mass adapter.

Select `--source mcmc` or `--source ns` if both are available. Without `--count`,
every saved posterior row is integrated with its posterior weight. With `--count`,
draws are sampled with replacement using their weights and receive equal weight;
repeated draws reuse their integration. `--seed` makes this selection reproducible.
Raw NS rows and log-weights are preferred when preserved; equal-weight samples are
never paired with raw weights. No likelihood values are exported because some old
sessions contain misaligned resampled likelihoods.

Free orbital parameters are read using saved packing metadata, including sparse
planet slots and h/k/lambda parameterization; fixed elements come from the session.
Do not change parameter labels, fixed elements, planet order, or configuration after
sampling: use the original sampled session. The adapter rejects detectable metadata
mismatches; historical sessions cannot detect every subsequent edit.

A draw is accepted only if every planet reaches the requested end time, every
recorded orbit remains finite and bound, eccentricity stays below the threshold,
and the fractional change in semi-major axis relative to its **initial** value
stays below the threshold. Events between output times may be missed. Discards,
mergers and truncated histories cannot count as stable. Timeout and integration
errors remain separate diagnostics, not evidence of physical instability.

Outputs in a **new** directory:

- `posterior.csv`: tested rows with original row index, weight and status.
- `stable_posterior.csv`: accepted rows with weights renormalized.
- `posterior.npz`: aligned arrays for plotting and downstream analysis.
- `diagnostics.json`: settings, seed, statuses, failure times and per-planet orbital ranges.
- `summary.csv`: weighted marginal quantiles before/after filtering. Angular
  intervals in this simple table are linear; inspect wraparound using circular statistics.

The result is a **stability-conditioned sampled posterior**, not a new N-body RV
fit, updated Bayesian evidence, or a proof of long-term stability. Inspect error
and timeout fractions before interpreting the accepted fraction.

The reusable API lives in `exostriker/lib/RV_mod/dynamical_posterior.py`:
`extract_posterior`, `orbital_draw`, `integrate_draw`, `integrate_posterior`, and
`export_posterior`. It has no Qt dependency and does not mutate the fit or sampler.
A future GUI worker can call it and relay `progress(done, total, record)` to a Qt
signal; the existing GUI tab is unchanged. Use `workers=N` in the API for parallel integrations in separate spawned
processes. Threads cannot parallelize the shared SWIFT native state. This example does not implement the deferred Stop feature.

## Output session

Both the command-line script and notebook also save `<session>_integrated.ses`
inside the output directory. Open it in Exo-Striker to inspect the augmented
sampler in **Sample corr.** and the corner plot. Like `make_stab.py`, it appends
per-planet integration columns and a final `max. time` column. The stability
colour overlay threshold is set just below the requested integration duration.
Timeout/error/incomplete rows have zero in this compatibility time column; their
actual status remains in `fit.dynamical_posterior['records']`.

The session retains all tested draws, with stable draws available through the
overlay and the separate `stable_posterior.csv`. Weighted tested rows are
systematically resampled for the GUI, which expects equal-weight samples; the
exact tested rows and weights remain in `fit.dynamical_posterior`. For NS the
original sampler/evidence is preserved in
`fit.dynamical_posterior_original_sampler`. The augmented sampler is a plotting
posterior, not a new nested-sampling run: its uniform log-weights are posterior
weights and must not be used for evidence calculations. Missing or untrustworthy
legacy likelihood alignment is represented by NaN and the median marker is used.
The input file is never overwritten, and an existing destination is rejected.

For use from code:

```python
from exostriker.lib.RV_mod.dynamical_posterior import save_integrated_session
save_integrated_session(fit, result, 'HD128311_integrated.ses')
```

## Default diagnostics: the compact make_stab.py list

By default, the output session and posterior tables append only:

- Mean outer/inner period ratio for the selected pair.
- Eccentricity mean and semi-amplitude for each enabled planet.
- Δω mean and semi-amplitude for the selected pair.
- Selected MMR angle means and semi-amplitudes, when `--mmr` is provided.
- `max. time` in the session for the existing stability overlay.

This follows `make_stab.py`, generalized to the enabled planets. Without `--pair`,
the first two enabled planets in increasing period order are used. Resonance
angles require an explicit pair and MMR; no resonance is guessed.

For the expanded set requested earlier, add **`--extended-diagnostics`**. Only
then are all orbital means/semi-amplitudes (a, e, i, P, omega, Omega, M, varpi,
lambda), all pair ratios/apsidal differences/mutual inclinations, angular
concentrations/excursions and history-coverage columns included. In Python use
`IntegrationOptions(..., extended_diagnostics=True)`. Default is `False`.
Semi-major axes are in AU, periods in days, and angles in degrees. Per-planet
min/max records remain in the JSON/session diagnostics rather than corner columns.

To inspect a particular resonance, supply **both** `--pair INNER OUTER` and
`--mmr P Q`. `P:Q` is the outer/inner period ratio, with P > Q. For the two
long-period planets in the supplied three-planet HD128311 session, from inside
`examples/stability`:

```bash
python3.10 integrate_posterior.py HD128311.ses \
  --output HD128311_resonance_out \
  --source ns --count 1000 --workers 8 \
  --years 1000 --timestep-days 0.5 --output-years 1 \
  --pair 2 3 --mmr 2 1
```

This requests the 2:1 MMR for GUI planets 2 and 3. For another pair or MMR, e.g.
planets 1 and 4 near 3:2, use `--pair 1 4 --mmr 3 2`. The inner planet must
come first. Use the simplest integer ratio (e.g. 2 1 rather than 4 2).
Without these arguments the compact default diagnostics are produced for the
first pair, without resonance-angle columns.

For order `r = P-Q`, the code evaluates `r+1` eccentricity-type angles:

```text
theta_k = Q*lambda_inner - P*lambda_outer
          + (r-k)*varpi_inner + k*varpi_outer,  k = 0,...,r
```

For 2:1 this gives `theta0 = lambda_inner - 2*lambda_outer + varpi_inner`
and `theta1 = lambda_inner - 2*lambda_outer + varpi_outer`. Columns include
`mean $\theta_{0,cd}$ (2:1) [deg]` and `ampl. $\theta_{0,cd}$ (2:1) [deg]`.
These are eccentricity-type resonant arguments; inclination-type and three-body
resonant arguments are not included.

**Angular semi-amplitudes** are half the shortest circular arc containing the
recorded values, so a 350°→0°→10° oscillation has a 10° semi-amplitude. With `--extended-diagnostics`, each
angular diagnostic also reports `resultant` (circular concentration from 0 to 1)
and `excursion` (the range of the unwrapped recorded angle, in degrees).
For a circulating angle the arc semi-amplitude typically approaches 180°; this
is not a measured libration amplitude or proof of resonance. Near-uniform angles
have no well-defined circular mean (reported as NaN when the resultant vanishes).

Summaries use the common recorded times of all planets, **strictly before the
first detected threshold failure**, or the full available common history when
there is no threshold failure. With `--extended-diagnostics`, `diagnostic_points` and `diagnostic_years` identify
the coverage. Means average the recorded values. Empty/failed histories produce
NaN diagnostics. Compare durations and statuses before combining these summaries.
Choose output cadence fine enough to resolve the angles: a coarse cadence can
alias circulation into apparent libration, and a short run may not cover a full
secular or resonant cycle. At e≈0, periastron angles are poorly defined; at i≈0,
nodes are poorly defined. Interpret those angle summaries accordingly.

In Python/notebooks use:

```python
options = IntegrationOptions(years=1000, timestep_days=0.5, output_years=1,
                             pair=(2, 3), mmr=(2, 1))
result = integrate_posterior(fit, options, source='ns', count=1000, workers=8)
```

## Planet designators in output labels

Session/corner labels and CSV/NPZ column labels follow the existing planet
parameter notation: GUI slots 1, 2, 3 correspond to b, c, d. For example,
`mean e$_c$`, `ampl. e$_c$`, `mean P$_d$/P$_c$`, and
`mean $\theta_{0,cd}$ (2:1) [deg]`. The `--pair` argument still takes numerical
GUI slots (`--pair 2 3` selects c and d). Internal JSON diagnostic keys retain
their machine-readable names. Regenerate the output session to obtain the new labels.

### Plotting angles across 0°/360°

Exo-Striker's corner and sample-correlation plots automatically put periodic
angle columns on a continuous branch centered on their circular mean. Thus a
posterior at 358°, 359°, 1°, 2° plots as −2°, −1°, 1°, 2°. Mean/median markers
and plotted intervals use the same branch. This is a display transformation;
saved posterior samples and semi-amplitudes are unchanged, and existing integrated
sessions do not need to be regenerated. Real separated modes remain separated.

## Sessions sampled in h/k/lambda

HKL sessions are integrated automatically using `h=e sin(omega)`,
`k=e cos(omega)`, and `lambda=M+omega`, whether or not you request additional
output columns. Free h/k/lambda values come from each posterior row; fixed
components come from the sampled session. `e` is `sqrt(h*h+k*k)` (not h*h+k*k).

To also include **initial/osculating derived e, omega, and mean anomaly** in
both the posterior tables and the output session, add these independent options:

```bash
python3.10 integrate_posterior.py MySystem.ses \
  --output MySystem_hkl_stability \
  --source ns --count 1000 --workers 8 \
  --years 1000 --timestep-days 0.5 --output-years 1 \
  --derive-e --derive-omega --derive-mean-anomaly
```

Run this command from `examples/stability`. For MCMC use `--source mcmc`.
`--derive-M0` is an alias for `--derive-mean-anomaly`. Omit any flag whose extra
columns you do not want. They default to off, preserving the compact output.
For a session sampled directly in e/omega/M, these flags add no duplicate columns.

Labels follow the original planet notation: `e$_b$`, `$\omega_b$ [deg]`, and
`MA$_b$ [deg]`. These are the derived **initial conditions**, distinct from the
integration-history summaries `mean e$_b$`, etc. Original h/k/lambda samples,
posterior weights and row alignment are retained. At exactly e=0, omega is
undefined (NaN in the derived output); the M conversion uses the omega=0
reference convention, so M=lambda. Integrations use that same circular reference.

In the notebook/API, use:

```python
options = IntegrationOptions(years=1000, timestep_days=0.5,
                             derive_e=True, derive_omega=True, derive_M0=True)
result = integrate_posterior(fit, options, source='ns', count=1000, workers=8)
```

These initial-condition output options are independent of `--extended-diagnostics`
and of `--pair`/`--mmr`.

## NV2016 resonance structure (Psi cos psi versus delta)

Append `--nv2016 --pair 1 2 --mmr 2 1` to also export the NV2016 resonance
variables, mapped trajectories, and a structure diagram. To make the epoch
posterior diagram from an original or integrated session **without rerunning
integrations**, use `--nv2016-only --pair 1 2 --mmr 2 1` and a new output directory.
See [the equations, output formats, examples and accuracy limits](../../docs/NV2016_RESONANCE.md).
The normal form assumes low eccentricities and near-coplanarity; its boundary
is approximate for eccentric giant-planet systems.
