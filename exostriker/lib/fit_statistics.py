"""Build the statistics dialog without importing GUI or numerical dependencies."""
import copy
import html
import math


COMPONENTS = ('RV', 'Transit', 'TTV', 'Astrometry')


def information_criteria(loglik, n, k):
    """Standard criteria using the full likelihood and unique free parameters."""
    if (loglik is None or n is None or k is None or n <= 0 or k < 0
            or not all(math.isfinite(value) for value in (loglik, n, k))):
        return {}
    return {'AIC': -2 * loglik + 2 * k,
            'BIC': -2 * loglik + k * math.log(n)}


def datasets_values(datasets):
    return datasets.values() if isinstance(datasets, dict) else datasets


def component_parameter_counts(prepared):
    """Nominal counts of free packed parameters entering each likelihood.

    This follows prepare_for_mcmc's packing order, rather than parameter labels
    (which users can rename). Shared parameters belong to multiple components.
    These counts are diagnostics, not leverage-based effective degrees of freedom.
    """
    free = set(int(index) for index in prepared.f_for_mcmc)
    roles = {name: set() for name in COMPONENTS}
    cursor = 0

    def block(size, names=()):
        nonlocal cursor
        indices = list(range(cursor, cursor + size))
        cursor += size
        for name in names:
            roles[name].update(indices)
        return indices

    rv_slots = prepared.mcmc_dataset_slots
    planet_slots = prepared.mcmc_planet_slots
    transit_slots = [i for i in range(len(prepared.tra_data_sets))
                     if len(prepared.tra_data_sets[i])]
    dynamical = bool(prepared.mod_dynamical)
    transits = bool(prepared.type_fit.get('Transit', False))
    free_transit_times = bool(prepared.get_TTVs[0])
    block(2 * len(rv_slots), ('RV',))
    phases = []
    for _ in planet_slots:
        orbit = block(7, ('TTV', 'Astrometry'))
        roles['RV'].update(orbit if dynamical else orbit[:5])
        if dynamical:
            roles['Transit'].update(orbit)
        else:
            roles['Transit'].update(orbit[index] for index in (1, 2, 3, 5))
        phases.append(orbit[4])
    block(2, ('RV',))  # RV linear and quadratic trends
    gp_sizes = {'RotKernel': 4, 'SHOKernel': 3, 'Matern32': 3,
                'RealTerm': 2, 'dSHOKernel': 6}
    rv_gp = block(gp_sizes[prepared.gp_kernel] if prepared.rtg[1] else 0, ('RV',))
    for planet, _ in enumerate(planet_slots):
        t0, scale, radius = block(3)
        roles['Transit'].add(radius)
        if not prepared.rho or planet == 0:
            roles['Transit'].add(scale)
        if not dynamical and not free_transit_times:
            roles['Transit'].add(t0)
        if transits and not dynamical:
            # model_loglik replaces M0/lambda by the phase derived from t0.
            for name in ('RV', 'TTV', 'Astrometry'):
                roles[name].discard(phases[planet])
                roles[name].add(t0)
    block(2 * len(transit_slots), ('Transit',))  # flux offsets and jitters
    transit_gp = block(gp_sizes[prepared.tra_gp_kernel] if prepared.rtg[3] else 0,
                       ('Transit',))
    for j, linked in enumerate(prepared.link_RV_GP):
        if linked and prepared.rtg[1] and prepared.rtg[3] and j < len(rv_gp) and j < len(transit_gp):
            roles['RV'].discard(rv_gp[j])
            roles['RV'].add(transit_gp[j])
    block(2 * len(transit_slots), ('Transit',))  # flux linear/quadratic trends
    block(len(planet_slots), ('RV',))  # apsidal precession in the RV model
    for i in transit_slots:
        if prepared.ld_gr[i] == i:
            block({'linear': 1, 'quadratic': 2, 'nonlinear': 4}.get(prepared.ld_m[i], 0),
                  ('Transit',))
    if free_transit_times:
        for i in planet_slots:
            block(len(prepared.tra_ttv[i]), ('Transit',))
    block(5, ('Astrometry',))  # catalogue position, parallax and proper-motion corrections
    stellar = block(1)
    if dynamical:
        roles['RV'].update(stellar)
        roles['Transit'].update(stellar)
    roles['TTV'].update(stellar)
    if not prepared.use_ast_hipp_gaia:
        roles['Astrometry'].update(stellar)
    if cursor != len(prepared.parameters):
        raise ValueError('Parameter layout does not match statistics accounting')
    return {name: len(indices & free) for name, indices in roles.items()}


def joint_information_criteria(fit):
    """Headline criteria without re-evaluating models or running integrations."""
    try:
        prepared = copy.deepcopy(fit)
        prepared.prepare_for_mcmc(rtg=prepared.rtg)
        k = len(prepared.par_for_mcmc)
        n = 0
        if fit.type_fit.get('RV', False):
            count = len(fit.fit_results.o_c)
            if not count:
                return {}
            n += count
        if fit.type_fit.get('Transit', False):
            count = int(fit.transit_results[4][4])
            if count <= 0:
                return {}
            n += count
        if fit.type_fit.get('TTV', False):
            count = sum(len(dataset[2]) * (1 + bool(dataset[-2]) + bool(dataset[-1]))
                        for dataset in datasets_values(fit.ttv_data_sets)
                        if len(dataset) and dataset[4])
            if not count:
                return {}
            n += count
        if fit.type_fit.get('AST', False):
            count = 0
            for datasets in (fit.ast_data_sets_hipp, fit.ast_data_sets_gaia):
                count += sum(len(dataset[0][-1])
                             for dataset in datasets_values(datasets) if len(dataset))
            if not count and not getattr(fit, 'use_ast_hipp_gaia', True):
                count = sum(len(dataset[2]) + len(dataset[4])
                            for dataset in datasets_values(fit.ast_data_sets)
                            if len(dataset) and dataset[6] and 1 <= dataset[5] <= fit.npl)
            if not count:
                return {}
            n += count
        return information_criteria(float(fit.loglik), n, k)
    except Exception:
        return {}


def gaussian_statistics(loglik, errors):
    """Recover chi-square from the normalized independent Gaussian likelihood."""
    errors = [float(error) for error in errors]
    if not errors or any(not math.isfinite(e) or e <= 0 for e in errors):
        raise ValueError('Measurement uncertainties must be finite and positive')
    loglik = float(loglik)
    if not math.isfinite(loglik):
        raise ValueError('Likelihood is not finite')
    normalization = math.fsum(math.log(2 * math.pi) + 2 * math.log(e) for e in errors)
    chi2 = -2 * loglik - normalization
    if chi2 < -1e-7 * max(1, abs(normalization), abs(loglik)):
        raise ValueError('Likelihood and uncertainties are inconsistent')
    return dict(lnL=loglik, N=len(errors), chi2=max(0.0, chi2))


def collect_fit_statistics(fit, backend):
    """Read RV/transit results and evaluate TTV/astrometry at current parameters.

    Parameter preparation runs on a copy, keeping the live session unchanged.
    Component DOF uses N_component - K_component as a nominal diagnostic;
    the total counts joint parameters once and is not the sum of component DOF.
    """
    rows = {name: {} for name in COMPONENTS}
    notes = []
    active = [name for name in COMPONENTS
              if fit.type_fit.get('AST' if name == 'Astrometry' else name, False)]
    prepared = None
    k = None
    try:
        prepared = copy.deepcopy(fit)
        prepared.prepare_for_mcmc(rtg=prepared.rtg)
        k = len(prepared.par_for_mcmc)
    except Exception as error:
        prepared = None
        notes.append(f'Parameter count unavailable: {error}')

    if 'RV' in active:
        try:
            result = fit.fit_results
            residuals = list(result.o_c)
            errors = [math.hypot(float(e), float(fit.jitt[int(i)]))
                      for e, i in zip(result.rv_err, result.idset)]
            if not residuals or len(residuals) != len(errors):
                raise ValueError('No matching RV residuals')
            chi2 = math.fsum((float(r) / e) ** 2 for r, e in zip(residuals, errors))
            normalization = math.fsum(math.log(2 * math.pi) + 2 * math.log(e) for e in errors)
            loglik = -0.5 * (chi2 + normalization)
            if fit.rtg[1]:
                # Use the full covariance likelihood, not the post-GP residual likelihood.
                gp_fit = copy.deepcopy(fit)
                backend.get_RV_gps_model(gp_fit, get_lnl=True)
                loglik = gp_fit.loglik
                chi2 = gp_fit.fit_results.chi2
                notes.append('RV GP chi² uses residuals after subtracting the GP prediction.')
            rows['RV'] = dict(lnL=loglik, N=len(residuals), chi2=float(chi2))
        except Exception as error:
            notes.append(f'RV statistics unavailable: {error}')

    if 'Transit' in active:
        try:
            result = fit.transit_results
            if int(result[4][4]) <= 0:
                raise ValueError('No fitted transit observations')
            rows['Transit'] = dict(lnL=float(result[0]), N=int(result[4][4]),
                                   chi2=float(result[4][0]))
            if fit.rtg[3]:
                notes.append('Transit GP chi² uses residuals after subtracting the GP prediction.')
        except Exception as error:
            notes.append(f'Transit statistics unavailable: {error}')

    if prepared is not None:
        rv_sets = prepared.rv_data_sets
        vel_files = [dataset for dataset in
                     (rv_sets.values() if isinstance(rv_sets, dict) else rv_sets)
                     if len(dataset)]
        common = (prepared.parameters, vel_files)
        model_fit = prepared.fit_results if fit.type_fit.get('RV', False) else False
        if 'TTV' in active:
            try:
                datasets = prepared.ttv_data_sets
                errors = []
                for dataset in datasets.values() if isinstance(datasets, dict) else datasets:
                    if not len(dataset) or not dataset[4]:
                        continue
                    errors.extend(dataset[2])
                    if dataset[-2]:
                        errors.extend(dataset[7])
                    if dataset[-1]:
                        errors.extend(dataset[9])
                loglik = backend.ttvs_loglik(*common, datasets, prepared.npl,
                    prepared.params.stellar_mass, prepared.ttv_times, prepared.hkl,
                    fit_results=model_fit, return_model=False)
                rows['TTV'] = gaussian_statistics(loglik, errors)
                notes.append('TTV N includes enabled depth and duration measurements, if present.')
            except Exception as error:
                notes.append(f'TTV statistics unavailable: {error}')
        if 'Astrometry' in active:
            try:
                errors, loglik = [], 0.0
                # These are the along-scan components used by model_loglik.
                for attribute, function in (
                        ('ast_data_sets_hipp', backend.ast_loglik_hipp),
                        ('ast_data_sets_gaia', backend.ast_loglik_gaia)):
                    datasets = getattr(prepared, attribute)
                    if not any(len(dataset) for dataset in datasets_values(datasets)):
                        continue
                    for dataset in datasets.values() if isinstance(datasets, dict) else datasets:
                        if len(dataset):
                            errors.extend(dataset[0][-1])
                    loglik += function(*common, datasets, prepared.npl,
                        prepared.params.stellar_mass, prepared.ast_times, prepared.hkl,
                        fit_results=model_fit, return_model=False)
                if not errors and not getattr(prepared, 'use_ast_hipp_gaia', True):
                    datasets = prepared.ast_data_sets
                    for dataset in datasets.values() if isinstance(datasets, dict) else datasets:
                        if len(dataset) and dataset[6] and 1 <= dataset[5] <= prepared.npl:
                            # The relative-astrometry likelihood has independent x/y terms.
                            errors.extend(dataset[2])
                            errors.extend(dataset[4])
                    loglik = backend.ast_loglik(*common, datasets, prepared.npl,
                        prepared.params.stellar_mass, prepared.ast_times, prepared.hkl,
                        fit_results=model_fit, return_model=False)
                rows['Astrometry'] = gaussian_statistics(loglik, errors)
            except Exception as error:
                notes.append(f'Astrometry statistics unavailable: {error}')

    counts = {}
    if prepared is not None:
        try:
            counts = component_parameter_counts(prepared)
        except Exception as error:
            notes.append(f'Component parameter counts unavailable: {error}')
    for name in active:
        row = rows[name]
        if name in counts:
            row['K'] = counts[name]
        if row.get('N') is not None and name in counts:
            row['DOF'] = row['N'] - counts[name]
            if row['DOF'] > 0 and row.get('chi2') is not None:
                row['red_chi2'] = row['chi2'] / row['DOF']
    total = {}
    if active:
        for field in ('lnL', 'N', 'chi2'):
            values = [rows[name].get(field) for name in active]
            if all(value is not None and math.isfinite(value) for value in values):
                total[field] = sum(values)
        if total.get('N') is not None and k is not None:
            total['DOF'] = total['N'] - k
            if total['DOF'] > 0 and total.get('chi2') is not None:
                total['red_chi2'] = total['chi2'] / total['DOF']
    if k is not None:
        total['K'] = k
    rows['Total'] = total
    total.update(information_criteria(total.get('lnL'), total.get('N'), k))
    return rows, k, notes


def statistics_html(rows, k, notes):
    def value(row, field):
        number = row.get(field)
        if number is None or not math.isfinite(number):
            return 'n.a.'
        return str(int(number)) if field in ('N', 'K', 'DOF') else f'{number:.5f}'
    headers = ('Component', 'lnL', 'N data', 'K', 'χ²', 'Reduced χ²', 'DOF', 'AIC', 'BIC')
    cells = [(name, *(value(row, field) for field in
              ('lnL', 'N', 'K', 'chi2', 'red_chi2', 'DOF', 'AIC', 'BIC')))
             for name, row in rows.items()]
    text = '<h3>Fit quality by component</h3>' + framed_table(headers, cells, bold_last=True)
    text += '<p>Joint free parameters K = ' + ('n.a.' if k is None else str(k))
    text += ('. Component DOF = N component − K component (nominal diagnostic); '
             'K component counts only free parameters entering that likelihood. Shared parameters '
             'count in each affected component but once in the total, so component DOF do not add up. '
             'Total DOF = N total − K joint. These are nominal counts, not effective DOF from '
             'the leverage or joint parameter covariance. '
             'Reduced χ² is n.a. when DOF ≤ 0. Unfitted or unavailable components are n.a.</p>')
    text += ('<p>Joint AIC = −2 lnL total + 2K; BIC = −2 lnL total + K ln(N total). '
             'Shown only for the total; no per-component penalty is assigned to shared parameters. '
             'K includes free orbital, stellar, offset, jitter, trend and GP parameters. '
             'Use maximum-likelihood fits and compare models on the same data, units and binning; '
             'lower values are preferred. Values at posterior means/medians are diagnostics, '
             'not maximum-likelihood criteria. BIC uses the number of scalar observations; '
             'with correlated data or GP models it is an approximation.</p>')
    for note in notes:
        text += '<p>' + html.escape(note) + '</p>'
    return text


def framed_table(headers, rows, bold_last=False):
    """Use table attributes supported by Qt's rich-text renderer."""
    text = '<table border="1" cellspacing="0" cellpadding="6" width="100%">'
    text += '<tr>' + ''.join(
        '<th align="left" bgcolor="#e8eef5" style="color: #182335;">'
        + html.escape(str(header)) + '</th>' for header in headers) + '</tr>'
    for i, row in enumerate(rows):
        text += '<tr>'
        for column, value in enumerate(row):
            cell = html.escape(str(value))
            if bold_last and i == len(rows) - 1:
                cell = '<b>' + cell + '</b>'
            alignment = 'left' if column == 0 else 'right'
            text += f'<td align="{alignment}">{cell}</td>'
        text += '</tr>'
    return text + '</table>'


def residual_statistics(residuals, errors):
    residuals = [float(value) for value in residuals]
    errors = [float(value) for value in errors]
    if not residuals or len(residuals) != len(errors):
        raise ValueError('No matching residuals and uncertainties')
    if any(not math.isfinite(r) for r in residuals) or any(
            not math.isfinite(e) or e <= 0 for e in errors):
        raise ValueError('Invalid residuals or uncertainties')
    rms = math.sqrt(math.fsum(r * r for r in residuals) / len(residuals))
    scale = min(errors)
    weights = [(scale / e) ** 2 for e in errors]
    wrms = math.sqrt(math.fsum(w * r * r for w, r in zip(weights, residuals))
                     / math.fsum(weights))
    return len(residuals), rms, wrms


def dataset_statistics_html(fit):
    sections = []
    for component, attribute, units in (
            ('RV', 'rv_data_sets', 'm/s'),
            ('Transit', 'tra_data_sets', 'relative flux')):
        if not fit.type_fit.get(component, False):
            continue
        datasets = getattr(fit, attribute, {})
        items = datasets.items() if isinstance(datasets, dict) else enumerate(datasets)
        rows = []
        for slot, dataset in items:
            if not len(dataset):
                continue
            label = str(dataset[4] if component == 'RV' else dataset[-1])
            try:
                if component == 'RV':
                    result = fit.fit_results
                    selected = [index for index, identity in enumerate(result.idset)
                                if int(identity) == slot]
                    residuals = [result.o_c[index] for index in selected]
                    errors = [result.rv_err[index] for index in selected]
                else:
                    # Dataset[4] can hold detrending residuals rather than fit residuals.
                    residuals = fit.transit_results[1][4][slot]
                    errors = dataset[2]
                n, rms, wrms = residual_statistics(residuals, errors)
                rows.append((label, n, f'{rms:.5f}', f'{wrms:.5f}', units))
            except (AttributeError, IndexError, KeyError, TypeError, ValueError):
                rows.append((label, 'n.a.', 'n.a.', 'n.a.', units))
        if rows:
            sections.append(f'<h3>{component} statistics by dataset</h3>'
                            + framed_table(('Dataset', 'N data', 'RMS', 'WRMS', 'Units'), rows))
    if sections:
        sections.append('<p>RMS and WRMS use fitted-model residuals before GP subtraction. '
                        'WRMS weights are 1/σ² using measurement uncertainties, '
                        'without added jitter. Unavailable fitted residuals are n.a.</p>')
    return ''.join(sections)
