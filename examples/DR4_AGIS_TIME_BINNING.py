#!/usr/bin/env python3
"""Extract Gaia CCD astrometry, optionally average by visit or TCB time window.

The eight output columns match DR4_AGIS_REDUCTION.py. XML input needs
numpy, astropy, pandas and gaiasupdate; reduced .d input needs only Python.
"""
import argparse
import math
from pathlib import Path


def read_reduced(path):
    header, rows = [], []
    with open(path) as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith('#'):
                header.append(line)
            elif line.strip():
                fields = line.split()
                if len(fields) != 8:
                    raise ValueError(f'Line {number}: expected eight columns')
                rows.append((int(fields[0]), int(fields[1]),
                             *(float(x) for x in fields[2:7]), int(fields[7])))
    return ''.join(header), rows


def read_xml(path, source_id, remove_model):
    import numpy as np
    import astropy.units as u
    from astropy.table import Table
    from astropy.time import Time
    from gaiasupdate.epoch_astrometry import GaiaEpochAstrometryArchive

    frame = Table.read(path).to_pandas()
    selected = frame[frame['source_id'] == source_id]
    if selected.empty:
        raise ValueError(f'Source {source_id} is absent from {path}')
    solution = GaiaEpochAstrometryArchive.supdate(frame, source_id)
    parameters = np.asarray(solution['parameters'])
    errors = np.asarray(solution['parameters_formal_uncertainty'])
    ra, dec = float(selected.iloc[0]['ra0']), float(selected.iloc[0]['dec0'])
    # Keep the twelve-line header layout expected by the ExoStriker reader.
    header = '\n'.join([
        '# Gaia DR4 pre-release epoch astrometry (TCB JD; angles in degrees)',
        f'# Source {source_id}', '# Exported for ExoStriker', '#', '#',
        '# HIP    MCE    NRES NC isol_n SCE  F2     F1', '# dummy',
        '# Hp      B-V    VarAnn NOB NR', '# dummy',
        '# RAdeg DEdeg Plx pm_RA pm_DE e_RA e_DE e_Plx e_pmRA e_pmDE',
        '# ' + ' '.join(str(x) for x in [ra, dec, *parameters[2:5], *errors[:5]])
        + ' ' + ' '.join(['---'] * 13), '#']) + '\n'
    rows = []
    epoch = Time(2010.0, format='jyear', scale='tcb')
    for _, visit in selected.iterrows():
        times = (epoch + np.asarray(visit['obs_time_tcb']) * u.ns).jd
        for ccd, used in enumerate(visit['used_by_agis_al']):
            if used:
                rows.append((int(visit['transit_id']), ccd, float(times[ccd]),
                             float(visit['centroid_pos_al'][ccd]),
                             float(visit['centroid_pos_error_al'][ccd]),
                             float(visit['parallax_factor_al']),
                             float(visit['scan_pos_angle'][ccd]), 0))
    if remove_model:
        residuals = np.asarray(solution['residuals']).reshape(-1)
        if len(residuals) != len(rows):
            raise ValueError('AGIS residual count does not match selected CCDs')
        rows = [(*row[:3], float(value), *row[4:])
                for row, value in zip(rows, residuals)]
    return header, rows


def angle_distance(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def bin_rows(rows, window, origin, max_scan_spread, mode='time'):
    """Use [origin+k*window, origin+(k+1)*window) bins, then split by angle.

    In time mode every pair of angles must be within max_scan_spread.
    Visit mode averages all valid CCDs with the same transit ID.
    None mode preserves individual valid CCD measurements and their IDs.
    Formal errors assume independent CCD uncertainties, without a noise floor.
    """
    if mode not in ('none', 'visit', 'time'):
        raise ValueError(f'Unknown binning mode: {mode}')
    clean = [row for row in rows if row[7] == 0 and row[4] > 0
             and all(math.isfinite(value) for value in row[2:7])]
    clean.sort(key=lambda row: row[2])
    if mode == 'none':
        return clean, len(clean)
    buckets = {}
    for row in clean:
        key = row[0] if mode == 'visit' else math.floor((row[2] - origin) / window)
        groups = buckets.setdefault(key, [])
        if mode == 'visit':
            if not groups:
                groups.append([])
            groups[0].append(row)
            continue
        for group in groups:
            if all(angle_distance(row[6], other[6]) <= max_scan_spread
                   for other in group):
                group.append(row)
                break
        else:
            groups.append([row])
    result = []
    for groups in buckets.values():
        for group in groups:
            # Scale weights to avoid overflow with very small uncertainties.
            minimum_error = min(row[4] for row in group)
            weights = [(minimum_error / row[4]) ** 2 for row in group]
            total = math.fsum(weights)
            def mean(column):
                base = group[0][column]
                return base + math.fsum(w * (row[column] - base)
                                        for w, row in zip(weights, group)) / total
            sine = math.fsum(w * math.sin(math.radians(row[6]))
                             for w, row in zip(weights, group))
            cosine = math.fsum(w * math.cos(math.radians(row[6]))
                               for w, row in zip(weights, group))
            scan = math.degrees(math.atan2(sine, cosine))
            # ID is the first contributing transit; CCD 0 denotes a bin.
            result.append((group[0][0], 0, mean(2), mean(3),
                           minimum_error / math.sqrt(total), mean(5), scan, 0))
    return sorted(result, key=lambda row: row[2]), len(clean)


def positive(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError('must be finite and positive')
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-p', '--path', type=Path, required=True,
                        help='raw XML or eight-column reduced .d file')
    parser.add_argument('-s', '--source_id', type=int, help='required for XML')
    parser.add_argument('-m', '--model_flag', type=str.lower, choices=['y', 'n'],
                        default='n', help='y: remove single-star model for XML (default n)')
    parser.add_argument('--bin-mode', choices=['none', 'visit', 'time'], default='time',
                        help='none: individual CCDs; visit: one row per transit ID; '
                             'time: fixed windows (default)')
    parser.add_argument('-b', '--bin-days', type=positive, default=0.01,
                        help='fixed window width in days (default 1; 0.25 = 6 hours)')
    parser.add_argument('--origin-jd', type=float, default=2451544.5,
                        help='TCB JD window origin (default J2000 midnight, 2451544.5)')
    parser.add_argument('--max-scan-spread', type=positive, default=1.0,
                        help='maximum pairwise scan-angle separation in degrees (default 1)')
    parser.add_argument('-o', '--output', type=Path,
                        help='output file or existing directory (default beside input)')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if not math.isfinite(args.origin_jd):
        parser.error('--origin-jd must be finite')
    if args.max_scan_spread >= 180:
        parser.error('--max-scan-spread must be less than 180 degrees')
    try:
        if args.path.suffix.lower() == '.xml':
            if args.source_id is None:
                parser.error('--source_id is required for XML input')
            header, rows = read_xml(args.path, args.source_id, args.model_flag == 'y')
            stem = f'dr4_reduction_id_{args.source_id}_' + (
                'no_sss' if args.model_flag == 'y' else 'with_sss')
        else:
            if args.model_flag == 'y':
                parser.error('-m y applies only to XML; supply an existing no_sss export instead')
            header, rows = read_reduced(args.path)
            stem = args.path.stem
        binned, accepted = bin_rows(rows, args.bin_days, args.origin_jd,
                                    args.max_scan_spread, args.bin_mode)
        if not binned:
            raise ValueError('No valid AGIS measurements with positive errors')
        suffix = {'none': 'unbinned', 'visit': 'visitbin',
                  'time': f'timebin_{args.bin_days:g}d'}[args.bin_mode]
        name = f'{stem}_{suffix}.d'
        output = args.output or args.path.with_name(name)
        if output.is_dir():
            output = output / name
        if output.resolve() == args.path.resolve():
            raise ValueError('Output must differ from input')
        with open(output, 'w' if args.overwrite else 'x') as stream:
            stream.write(header)
            for row in binned:
                stream.write(' '.join(str(value) for value in row) + '\n')
        label = 'CCD measurements' if args.bin_mode == 'none' else 'bins'
        print(f'{accepted} valid CCD measurements -> {len(binned)} {label}: {output}')
        if accepted != len(rows):
            print(f'Excluded {len(rows) - accepted} flagged or invalid measurements')
        if args.bin_mode != 'none':
            print('Errors are formal independent-measurement errors; no correlated-noise floor applied.')
    except (OSError, ValueError, ImportError) as error:
        parser.exit(1, f'Error: {error}\n')


if __name__ == '__main__':
    main()
