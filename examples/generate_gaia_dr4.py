#!/usr/bin/env python3
"""Generate maximum-planet CSV models, optionally using sym_small_error scans."""
import argparse
import csv
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'exostriker/lib'))
from gaia_synthetic import GetGaiaDR4, HEADER_FIELDS, REFERENCE_EPOCH


# Explicit aliases between the supplied CSV and reference dataset filenames.
REFERENCE_NAMES = {
    'GJ179': 'GJ 179', 'GJ2126': 'GJ 2126', 'GJ317': 'GJ 317',
    'GJ433': 'GJ 433', 'GJ649': 'GJ 649', 'GJ676': 'GJ 676 A',
    'GJ83.1': 'GJ 9066', 'GJ849': 'GJ 849', 'GL687': 'GJ 687',
    'HD1326': 'GJ 15 A', 'HD204961': 'GJ 832', 'HIP109638': 'GJ 9773',
    'HIP57050': 'GJ 1148', 'HIP79431': 'HIP 79431',
    'J06024+498': 'G 192-15', 'J08413+594': 'GJ 3512',
    'J12230+640': 'GJ 463', 'J13209+342': 'GJ 508.2',
    'J21221+229': 'TYC 2187-512-1',
}


def read_models(path):
    """Select the largest complete model per target; reject ambiguous ties."""
    groups = {}
    with Path(path).open(newline='') as stream:
        for row in csv.DictReader(stream):
            models = groups.setdefault(row['target'], {})
            models.setdefault(row['model'], []).append(row)
    selected = {}
    for target, models in groups.items():
        maximum = max(int(rows[0]['n_planets']) for rows in models.values())
        candidates = [rows for rows in models.values()
                      if int(rows[0]['n_planets']) == maximum]
        if len(candidates) != 1:
            raise ValueError(f'{target}: multiple models have {maximum} planets.')
        rows = sorted(candidates[0], key=lambda row: int(row['planet_slot']))
        if [int(row['planet_slot']) for row in rows] != list(range(1, maximum + 1)):
            raise ValueError(f'{target}: model has missing or duplicate planet slots.')
        for field in ('n_planets', 'epoch_BJD', 'gaia_dr3_id', 'stellar_mass_Msun'):
            if len({row[field] for row in rows}) != 1:
                raise ValueError(f'{target}: inconsistent {field} within the model.')
        selected[target] = rows
    return selected


def orbital_parameters(rows):
    """Map the CSV's physical columns to GetModel's named orbital inputs."""
    columns = dict(K='K_m_s', P='P_days', e='eccentricity', omega='omega_deg',
                   M0='M0_deg', inclination='inclination_deg', node='Omega_deg')
    return {f"Pl_{row['planet_slot']}": {
        key: float(row[column]) for key, column in columns.items()} for row in rows}


def read_reference(path):
    """Reuse the reference header, observing geometry and per-row errors."""
    with Path(path).open() as stream:
        for line in stream:
            names = line.lstrip('# ').split()
            if names[:3] == ['RAdeg', 'DEdeg', 'Plx']:
                values = next(stream).lstrip('# ').split()
                catalogue = {key: float(values[names.index(key)]) for key in HEADER_FIELDS}
                break
        else:
            raise ValueError(f'{path}: no stellar solution header.')
    # These reference centroids were constructed relative to J2017.5.
    catalogue['ref_epoch'] = 2017.5
    data = np.loadtxt(path, ndmin=2)
    orbits = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith('#'):
            continue
        fields = line[1:].split()
        if len(fields) == 7:
            orbits.append([float(value) for value in fields])
    return catalogue, data, orbits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', type=Path, default=ROOT / 'examples/orbital_results_from_json_2.csv')
    parser.add_argument('--output', type=Path, default=ROOT / 'examples/results/gaia_dr4_csv')
    parser.add_argument('--reference-directory', type=Path,
                        help='offline: reuse sym_small_error catalogue, scans and uncertainties')
    parser.add_argument('--target', action='append', help='CSV target name; repeat to select several')
    parser.add_argument('--sigma-al', type=float, default=0.05,
                        help='mas per nominal transit in archive mode (default 0.05)')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--no-noise', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    models = read_models(args.csv)
    if args.target:
        unknown = set(args.target) - models.keys()
        if unknown:
            parser.error('Unknown targets: ' + ', '.join(sorted(unknown)))
    args.output.mkdir(parents=True, exist_ok=True)
    report = []
    # Spawn by the complete sorted target list so subsets retain the same seeds.
    children = np.random.SeedSequence(args.seed).spawn(len(models))
    for (target, rows), child in zip(sorted(models.items()), children):
        if args.target and target not in args.target:
            continue
        name = REFERENCE_NAMES.get(target, target)
        status = dict(target=target, model=rows[0]['model'], n_planets=len(rows),
                      source_id=rows[0]['gaia_dr3_id'], status='', output='',
                      reference='', rms_difference_mas='',
                      stored_orbit_residual_rms_mas='')
        options = {}
        reference = None
        if args.reference_directory:
            path = args.reference_directory / f'{name}_dr4sym.d'
            if not path.exists():
                status['status'] = 'no reference file; use archive mode'
                report.append(status)
                print(f'{target}: {status["status"]}')
                continue
            catalogue, reference, stored_orbits = read_reference(path)
            options = dict(catalogue=catalogue,
                           scan_data=(reference[:, 2], np.deg2rad(reference[:, 6]), reference[:, 5]),
                           sigma_al=reference[:, 4])
            status['reference'] = str(path)
        else:
            options['sigma_al'] = args.sigma_al
        output = args.output / f'{name}_dr4sym.d'
        if output.resolve().parent == (ROOT / 'examples/sym_small_error').resolve():
            parser.error('Use a separate output directory to preserve the original reference files.')
        result = GetGaiaDR4(
            orbital_parameters(rows), source_id=rows[0]['gaia_dr3_id'],
            epoch=float(rows[0]['epoch_BJD']), stellar_mass=float(rows[0]['stellar_mass_Msun']),
            seed=int(child.generate_state(1)[0]), add_noise=not args.no_noise,
            output=output, overwrite=args.overwrite, **options)
        status['status'] = 'generated'
        status['output'] = str(output)
        if reference is not None:
            scan = np.deg2rad(reference[:, 6])
            epoch_shift = (catalogue['ref_epoch'] - REFERENCE_EPOCH) * (
                catalogue['pm_RA'] * np.sin(scan) + catalogue['pm_DE'] * np.cos(scan))
            aligned_reference = reference[:, 3] + epoch_shift
            status['rms_difference_mas'] = float(
                np.sqrt(np.mean((result.model - aligned_reference) ** 2)))
            # Diagnose the old generator using its own explicitly saved stellar
            # orbits, independent of the updated CSV fits and omega conversion.
            from ExTRA.astrometry import orbit_total

            east, north = orbit_total(stored_orbits, reference[:, 2])
            stored_model = result.baseline_al - epoch_shift + east * np.sin(scan) + north * np.cos(scan)
            status['stored_orbit_residual_rms_mas'] = float(
                np.sqrt(np.mean((stored_model - reference[:, 3]) ** 2)))
        report.append(status)
        print(f'{target}: {len(rows)} planets -> {output}')
    if report:
        with (args.output / 'generation_report.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(report[0]))
            writer.writeheader()
            writer.writerows(report)


if __name__ == '__main__':
    main()
