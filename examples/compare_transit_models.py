"""Run inside ExoStriker's Jupyter shell with its existing `fit` object:

%run -i /home/trifonov/exostriker_codex/examples/compare_transit_models.py

Choose a zero-based dataset slot below, or keep None for the first loaded one.
The comparison preserves fitted parameters, the old GUI results, and dense curves;
new predictions are stored separately at the actual transit data timestamps.
"""
import numpy as np
import matplotlib.pyplot as plt

# Change to 0, 1, ..., 59 to select a specific transit dataset.
dataset = None
if dataset is None:
    dataset = next((i for i in range(60) if len(fit.tra_data_sets[i])), None)
if dataset is None:
    raise ValueError('Load a transit dataset first.')
if not len(fit.transit_results):
    raise ValueError('Initialize the old GUI transit model first.')

# The legacy GUI structure preserves the original dataset slots.
old_times = np.asarray(fit.transit_results[1][0][dataset], dtype=float)
old_flux = np.asarray(fit.transit_results[1][3][dataset], dtype=float)
if old_times.size == 0 or old_times.shape != old_flux.shape:
    raise ValueError('The selected dataset has no matching old model; initialize it first.')

# Match batman's GUI t0 parameterization for Keplerian models. N-body initial
# conditions remain M0 at the integration epoch.
phase = 'M0' if fit.mod_dynamical else 't0'
fit.update_transit_data_models(dataset=dataset, phase=phase)
new_times = fit.fit_results.transit_data_jd[dataset]
if not np.array_equal(new_times, old_times):
    raise ValueError('Old model timestamps differ from the current dataset; initialize the old GUI model first.')
# Include the same baseline/dilution as the legacy GUI model, without the
# observation-frame offset. No interpolation or dense epoch grid is used.
new_flux = fit.fit_results.transit_data_models_with_baseline[dataset]
finite = np.isfinite(old_times) & np.isfinite(old_flux) & np.isfinite(new_flux)
if not np.any(finite):
    raise ValueError('No finite points to compare.')
order = np.argsort(old_times[finite], kind='stable')
comparison_times = old_times[finite][order]
comparison_old = old_flux[finite][order]
comparison_new = new_flux[finite][order]
difference_ppm = (comparison_new - comparison_old) * 1e6

print(f'Transit dataset slot {dataset} (GUI dataset {dataset + 1})')
print(f'Compared points: {len(difference_ppm)}')
print(f'Mean (new - old): {np.mean(difference_ppm):.6g} ppm')
print(f'RMS difference:   {np.sqrt(np.mean(difference_ppm**2)):.6g} ppm')
print(f'Max |difference|: {np.max(np.abs(difference_ppm)):.6g} ppm')

origin = comparison_times[0]
fig, axes = plt.subplots(2, 1, sharex=True, figsize=(10, 6))

axes[0].plot(comparison_times-origin, comparison_old, '-', label='Old GUI model')
axes[0].plot(comparison_times-origin, comparison_new, '--', label=f'New Fortran ({phase}): same timestamps')
axes[0].set_ylabel('Normalized flux')
axes[0].legend()
axes[1].plot(comparison_times-origin, difference_ppm, '.-', markersize=3)
axes[1].axhline(0., color='black', linewidth=.7)
axes[1].set_ylabel('New - old [ppm]')
axes[1].set_xlabel(f'Time - {origin:.8f} [days]')
fig.suptitle(f'Transit model comparison: dataset {dataset + 1}')
fig.tight_layout()
plt.show()

print('Residuals use identical timestamps, not interpolated stored curves.')
print(f'Phase: {phase}; baseline and dilution included to match the old GUI frame.')
print('Derived M0 at RV epoch [deg]:', fit.fit_results.transit_data_mean_anomaly)
