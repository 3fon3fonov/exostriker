#!/usr/bin/env python3
"""Integrate sessions using this checkout, or the installed package when copied.

Provided by Trifon Trifonov with the aid of CODEX; use with caution.
This script can be copied and run independently of the examples directory.
"""
from pathlib import Path
import argparse
import sys

# Prefer this checkout when the example is run from its source tree.
# Standalone copies use the installed package; spawned workers inherit sys.path.
_checkout = Path(__file__).resolve().parent.parent.parent
if (_checkout / 'pyproject.toml').is_file() and (
        _checkout / 'exostriker/lib/RV_mod/dynamical_posterior.py').is_file():
    sys.path.insert(0, str(_checkout))

import dill
import exostriker
from exostriker.lib.RV_mod.dynamical_posterior import (
    INTEGRATORS, IntegrationOptions, integrate_posterior, export_posterior, save_integrated_session)


def print_nv2016_summary(result):
    records = [r.get('resonance_nv2016', {}) for r in result['records']]
    mapped = sum('summary' in r for r in records)
    limited = sum(any(not note.startswith('Epoch-only') for note in r.get('warnings', [])) for r in records)
    print(f'NV2016: mapped {mapped}/{len(records)} draws; {limited} have model-applicability notes. See resonance.csv.', flush=True)


def main():
    parser = argparse.ArgumentParser(description='Integrate sampled Exo-Striker orbits with SWIFT or REBOUND; export stability-conditioned posteriors.')
    parser.add_argument('session', type=Path)
    parser.add_argument('--output', type=Path, required=True, help='New output directory')
    parser.add_argument('--source', choices=['auto', 'ns', 'mcmc'], default='auto')
    parser.add_argument('--count', type=int, help='Weighted posterior draws; default: every row')
    parser.add_argument('--derive-e', action='store_true', help='Append e=sqrt(h^2+k^2) for HKL samples')
    parser.add_argument('--derive-omega', action='store_true', help='Append omega=atan2(h,k) for HKL samples')
    parser.add_argument('--derive-M0', '--derive-mean-anomaly', dest='derive_M0', action='store_true', help='Append M=lambda-omega for HKL samples')
    parser.add_argument('--nv2016', action='store_true', help='Export NV2016 Psi, psi, delta, mapped histories and resonance-structure plot; requires first-order --mmr and --pair')
    parser.add_argument('--nv2016-only', action='store_true', help='Export only epoch NV2016 coordinates/structure from original or integrated session; no new integrations')
    parser.add_argument('--nv2016-points', type=int, default=256, help='Maximum mapped trajectory points saved per unique posterior draw (default: 256; summaries use all output samples)')
    parser.add_argument('--extended-diagnostics', action='store_true', help='Include all orbital/pair means, amplitudes and angular diagnostics; default: compact make_stab.py columns')
    parser.add_argument('--pair', nargs=2, type=int, default=(), metavar=('INNER', 'OUTER'), help='1-based GUI planet slots for MMR')
    parser.add_argument('--mmr', nargs=2, type=int, default=(), metavar=('P', 'Q'), help='Outer/inner period ratio P:Q, e.g. --mmr 2 1')
    parser.add_argument('--workers', type=int, default=1, help='Parallel CPU processes (default: 1)')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--years', type=float, default=1000)
    parser.add_argument('--timestep-days', type=float, default=1)
    parser.add_argument('--output-years', type=float, default=1)
    parser.add_argument('--timeout-seconds', type=float, default=60)
    parser.add_argument('--integrator', choices=INTEGRATORS, default='symba',
                        help='Integrator name selects SWIFT or REBOUND automatically (default: symba)')
    parser.add_argument('--max-fractional-a-change', type=float, default=.1)
    parser.add_argument('--max-eccentricity', type=float, default=.95)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be a positive integer')
    if args.output.exists():
        parser.error('Output directory already exists; choose a new directory')
    print(f'Using Exo-Striker: {exostriker.__file__}', flush=True)
    with args.session.open('rb') as stream:
        fit = dill.load(stream)
    args.nv2016 = args.nv2016 or args.nv2016_only
    if args.nv2016 and 'nv2016' not in IntegrationOptions.__dataclass_fields__:
        parser.error('The installed Exo-Striker lacks NV2016 support; install an updated version first.')
    options = IntegrationOptions(**{k: getattr(args, k) for k in IntegrationOptions.__dataclass_fields__})
    if args.nv2016_only:
        from exostriker.lib.RV_mod.dynamical_posterior import map_nv2016_posterior
        from exostriker.lib.RV_mod.resonance import export_resonance
        result = map_nv2016_posterior(fit, options, args.source, args.count, args.seed)
        args.output.mkdir(parents=True, exist_ok=False)
        export_resonance(result, args.output)
        print_nv2016_summary(result)
        print('Epoch resonance mapping only; no new integrations performed.')
        print(args.output.resolve())
        return
    result = integrate_posterior(fit, options, args.source, args.count, args.seed,
                                progress=lambda n, total, r: print(f'{n}/{total}: {r["status"]}', flush=True), workers=args.workers)
    export_posterior(result, args.output)
    if args.nv2016:
        print_nv2016_summary(result)
    save_integrated_session(fit, result, args.output / (args.session.stem + '_integrated.ses'))
    print(f'Stable weight fraction: {result["weights"][result["stable"]].sum():.6g}')
    print(args.output.resolve())


if __name__ == '__main__':
    main()
