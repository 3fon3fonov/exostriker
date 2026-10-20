# Exo-Striker examples

**Examples provided by Trifon Trifonov with the aid of CODEX; use with caution.**

These Python scripts and matching Jupyter notebooks follow the fitting → priors →
nested sampling → diagnostics → posterior export workflow in
`Alans_NB/NS_nu_Oph_eng_v3.ipynb`, using the current library API and bundled data.
They build fresh sessions, so no historical `.ses` file is required.

| Target | Script | Notebook | Analysis |
|---|---|---|---|
| HIP 5364 / Eta Ceti | [HIP5364_RV.py](HIP5364_RV.py) | [HIP5364_RV.ipynb](HIP5364_RV.ipynb) | Two RV planets; Keplerian then dynamical fit |
| 51 Peg | [51Peg_RV.py](51Peg_RV.py) | [51Peg_RV.ipynb](51Peg_RV.ipynb) | Circular one-planet RV fit; separate ELODIE/HARPS offsets and jitters |
| WASP-6 | [WASP6_transit.py](WASP6_transit.py) | [WASP6_transit.ipynb](WASP6_transit.ipynb) | Detrended light curve; period, transit geometry, offset and white noise |

Install the checkout into your Python environment from the repository root:

```bash
python -m pip install -e .
python examples/HIP5364_RV.py
python examples/51Peg_RV.py
python examples/WASP6_transit.py
```

The first import may compile the native backend and requires a working Fortran
compiler, Meson and Ninja. Plots use Matplotlib's own math text; external LaTeX is
not required. To display saved figures interactively, add `--show`. For headless
runs set `MPLBACKEND=Agg` and `QT_QPA_PLATFORM=offscreen`.

Each default run optimizes, saves model/residual and phase plots, and writes a
GUI-readable session under `examples/results/TARGET/`. To choose a destination,
use `--output /path/to/results`. Output files with the same name are replaced.

Nested sampling is explicit because it can take minutes to hours:

```bash
python examples/51Peg_RV.py --nested --threads 1
```

Open the notebooks from the repository root or `examples/` and select a kernel
with the same dependencies. Their default **Run All** performs the best-fit
analysis; change `RUN_NESTED=True` to enable the sampling cells. Start with one
process, especially in notebooks and on macOS/Windows.

The nested-sampling workflow checks finite prior bounds and starting values,
prints evidence and effective sample size, and writes:

- `TARGET_posterior.csv`: equal-weight parameter samples and aligned `logL`;
- `TARGET_posterior_raw.csv`: raw samples, `logL`, `logwt`, and normalized weights;
- `TARGET_posterior_summary.csv`: original labels, median, 68.3% intervals,
  maximum-likelihood parameters, evidence and ESS;
- `TARGET_corner.png` / `.pdf`: posterior corner plot;
- a `.ses` session retaining the sampler.

The current library retains aligned raw arrays in `samples_ns` and `logl_ns`.
The export helper resamples those arrays together with a fixed resampling seed.
This fixes the historical mismatch discussed in the ν Oph notebook. The seed
makes exports reproducible for a given sampler result; it does not seed dynesty.
If a short run leaves a posterior column with zero range, the helper warns and
skips the corner plot; it still saves the CSV diagnostics and sampler session.

Priors, fixed stellar masses, orbital constraints and limb darkening are explicit
**tutorial assumptions**. Inspect residuals, optimizer diagnostics and convergence
before scientific use. Do not treat an abbreviated sampling run as a resolved
posterior or evidence measurement. No full nested-sampling run is performed just
by opening or running all cells with the defaults.

Shared plotting and export code lives in [example_tools.py](example_tools.py).
Historical examples, session files and generated plots have been preserved in
[legacy/](legacy/) for reference; start with the tutorials above.
