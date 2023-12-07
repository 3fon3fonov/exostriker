from __future__ import division, print_function
from os import path
#import multiprocessing
from pathos import multiprocessing
import numpy
import sys
import warnings
from functools import partial
from tqdm import tqdm

# TLS parts
from transitleastsquares_ES.results import transitleastsquaresresults
import transitleastsquares_ES.tls_constants as tls_constants
from transitleastsquares_ES.stats import (
    FAP,
    rp_rs_from_depth,
    period_uncertainty,
    spectra,
    final_T0_fit,
    model_lightcurve,
    all_transit_times,
    calculate_transit_duration_in_days,
    calculate_stretch,
    calculate_fill_factor,
    intransit_stats,
    snr_stats,
    count_stats,
)
from transitleastsquares_ES.catalog import catalog_info
from transitleastsquares_ES.helpers import resample, transit_mask
from transitleastsquares_ES.helpers import impact_to_inclination
from transitleastsquares_ES.grid import duration_grid, period_grid
from transitleastsquares_ES.core import (
    edge_effect_correction,
    lowest_residuals_in_this_duration,
    out_of_transit_residuals,
    fold,
    foldfast,
    search_period,
)
from transitleastsquares_ES.transit import reference_transit, fractional_transit, get_cache
from transitleastsquares_ES.validate import validate_inputs, validate_args


class transitleastsquares_ES(object):
    """Compute the transit least squares of limb-darkened transit models"""

    def __init__(self, t, y, dy=None, verbose=True):
        self.t, self.y, self.dy = validate_inputs(t, y, dy)
        self.verbose = verbose

    def power(self, **kwargs):
        """Compute the periodogram for a set of user-defined parameters"""
        self, kwargs = validate_args(self, kwargs)

        if self.verbose:
            print(tls_constants.TLS_VERSION)
        

        periods = period_grid(
            R_star=self.R_star,
            M_star=self.M_star,
            time_span=numpy.max(self.t) - numpy.min(self.t),
            period_min=self.period_min,
            period_max=self.period_max,
            oversampling_factor=self.oversampling_factor,
            n_transits_min=self.n_transits_min,
        )

        durations = duration_grid(
            periods, shortest=1 / len(self.t), log_step=self.duration_grid_step
        )

        maxwidth_in_samples = int(numpy.max(durations) * numpy.size(self.y))
        if maxwidth_in_samples % 2 != 0:
            maxwidth_in_samples = maxwidth_in_samples + 1
        lc_cache_overview, lc_arr = get_cache(
            durations=durations,
            maxwidth_in_samples=maxwidth_in_samples,
            per=self.per,
            rp=self.rp,
            a=self.a,
            inc=self.inc,
            ecc=self.ecc,
            w=self.w,
            u=self.u,
            limb_dark=self.limb_dark,
            verbose=self.verbose
        )

        if self.verbose:
            print(
                "Searching "
                + str(len(self.y))
                + " data points, "
                + str(len(periods))
                + " periods from "
                + str(round(min(periods), 3))
                + " to "
                + str(round(max(periods), 3))
                + " days"
            )

        # Python 2 multiprocessing with "partial" doesn't work
        # For now, only single-threading in Python 2 is supported
        if sys.version_info[0] < 3:
            self.use_threads = 1
            warnings.warn("This TLS version supports no multithreading on Python 2")

        if self.verbose:
            if self.use_threads == multiprocessing.cpu_count():
                print("Using all " + str(self.use_threads) + " CPU threads")
            else:
                print(
                    "Using "
                    + str(self.use_threads)
                    + " of "
                    + str(multiprocessing.cpu_count())
                    + " CPU threads"
                )

        if self.show_progress_bar:
            bar_format = "{desc}{percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} periods | {elapsed}<{remaining}"
            pbar = tqdm(total=numpy.size(periods), smoothing=0.3, bar_format=bar_format)

        if tls_constants.PERIODS_SEARCH_ORDER == "ascending":
            periods = reversed(periods)
        elif tls_constants.PERIODS_SEARCH_ORDER == "descending":
            pass  # it already is
        elif tls_constants.PERIODS_SEARCH_ORDER == "shuffled":
            periods = numpy.random.permutation(periods)
        else:
            raise ValueError("Unknown PERIODS_SEARCH_ORDER")

        # Result lists now (faster), convert to numpy array later
        test_statistic_periods = []
        test_statistic_residuals = []
        test_statistic_rows = []
        test_statistic_depths = []

        if self.use_threads > 1:  # Run multi-core search
            pool = multiprocessing.Pool(processes=self.use_threads)
            params = partial(
                search_period,
                t=self.t,
                y=self.y,
                dy=self.dy,
                transit_depth_min=self.transit_depth_min,
                R_star_min=self.R_star_min,
                R_star_max=self.R_star_max,
                M_star_min=self.M_star_min,
                M_star_max=self.M_star_max,
                lc_arr=lc_arr,
                lc_cache_overview=lc_cache_overview,
                T0_fit_margin=self.T0_fit_margin,
            )
            for data in pool.imap_unordered(params, periods):
                test_statistic_periods.append(data[0])
                test_statistic_residuals.append(data[1])
                test_statistic_rows.append(data[2])
                test_statistic_depths.append(data[3])
                if self.show_progress_bar:
                    pbar.update(1)
            pool.close()
        else:
            for period in periods:
                data = search_period(
                    period=period,
                    t=self.t,
                    y=self.y,
                    dy=self.dy,
                    transit_depth_min=self.transit_depth_min,
                    R_star_min=self.R_star_min,
                    R_star_max=self.R_star_max,
                    M_star_min=self.M_star_min,
                    M_star_max=self.M_star_max,
                    lc_arr=lc_arr,
                    lc_cache_overview=lc_cache_overview,
                    T0_fit_margin=self.T0_fit_margin,
                )
                test_statistic_periods.append(data[0])
                test_statistic_residuals.append(data[1])
                test_statistic_rows.append(data[2])
                test_statistic_depths.append(data[3])
                if self.show_progress_bar:
                    pbar.update(1)

        if self.show_progress_bar:
            pbar.close()

        # imap_unordered delivers results in unsorted order ==> sort
        test_statistic_periods = numpy.array(test_statistic_periods)
        sort_index = numpy.argsort(test_statistic_periods)
        test_statistic_periods = test_statistic_periods[sort_index]
        test_statistic_residuals = numpy.array(test_statistic_residuals)[sort_index]
        test_statistic_rows = numpy.array(test_statistic_rows)[sort_index]
        test_statistic_depths = numpy.array(test_statistic_depths)[sort_index]

        idx_best = numpy.argmin(test_statistic_residuals)
        best_row = test_statistic_rows[idx_best]
        duration = lc_cache_overview["duration"][best_row]
        maxwidth_in_samples = int(numpy.max(durations) * numpy.size(self.t))

        if max(test_statistic_residuals) == min(test_statistic_residuals):
            no_transits_were_fit = True
            warnings.warn('No transit were fit. Try smaller "transit_depth_min"')
        else:
            no_transits_were_fit = False

        # Power spectra variants
        chi2 = test_statistic_residuals
        degrees_of_freedom = 4
        chi2red = test_statistic_residuals / (len(self.t) - degrees_of_freedom)
        chi2_min = numpy.min(chi2)
        chi2red_min = numpy.min(chi2red)

        if no_transits_were_fit:
            power_raw = numpy.zeros(len(chi2))
            power = numpy.zeros(len(chi2))
            period = numpy.nan
            depth = 1
            SR = 0
            SDE = 0
            SDE_raw = 0
            T0 = 0
            transit_times = numpy.nan
            transit_duration_in_days = numpy.nan
            internal_samples = (
                int(len(self.y)) * tls_constants.OVERSAMPLE_MODEL_LIGHT_CURVE
            )
            folded_phase = numpy.nan
            folded_y = numpy.nan
            folded_dy = numpy.nan
            model_folded_phase = numpy.nan
            model_folded_model = numpy.nan
            model_transit_single = numpy.nan
            model_lightcurve_model = numpy.nan
            model_lightcurve_time = numpy.nan
            depth_mean_odd = numpy.nan
            depth_mean_even = numpy.nan
            depth_mean_odd_std = numpy.nan
            depth_mean_even_std = numpy.nan
            all_flux_intransit_odd = numpy.nan
            all_flux_intransit_even = numpy.nan
            per_transit_count = numpy.nan
            transit_depths = numpy.nan
            transit_depths_uncertainties = numpy.nan
            all_flux_intransit = numpy.nan
            snr_per_transit = numpy.nan
            snr_pink_per_transit = numpy.nan
            depth_mean = numpy.nan
            depth_mean_std = numpy.nan
            snr = numpy.nan
            rp_rs = numpy.nan
            depth_mean_odd = numpy.nan
            depth_mean_even = numpy.nan
            depth_mean_odd_std = numpy.nan
            depth_mean_even_std = numpy.nan
            odd_even_difference = numpy.nan
            odd_even_std_sum = numpy.nan
            odd_even_mismatch = numpy.nan
            transit_count = numpy.nan
            empty_transit_count = numpy.nan
            distinct_transit_count = numpy.nan
            duration = numpy.nan
            in_transit_count = numpy.nan
            after_transit_count = numpy.nan
            before_transit_count = numpy.nan
        else:
            SR, power_raw, power, SDE_raw, SDE = spectra(chi2, self.oversampling_factor)
            index_highest_power = numpy.argmax(power)
            period = test_statistic_periods[index_highest_power]
            depth = test_statistic_depths[index_highest_power]
            T0 = final_T0_fit(
                signal=lc_arr[best_row],
                depth=depth,
                t=self.t,
                y=self.y,
                dy=self.dy,
                period=period,
                T0_fit_margin=self.T0_fit_margin,
                show_progress_bar=self.show_progress_bar,
                verbose=self.verbose
            )
            transit_times = all_transit_times(T0, self.t, period)

            transit_duration_in_days = calculate_transit_duration_in_days(
                self.t, period, transit_times, duration
            )
            phases = fold(self.t, period, T0=T0 + period / 2)
            sort_index = numpy.argsort(phases)
            folded_phase = phases[sort_index]
            folded_y = self.y[sort_index]
            folded_dy = self.dy[sort_index]
            # Model phase, shifted by half a cadence so that mid-transit is at phase=0.5
            model_folded_phase = numpy.linspace(
                0 + 1 / numpy.size(self.t) / 2,
                1 + 1 / numpy.size(self.t) / 2,
                numpy.size(self.t),
            )
            # Folded model / model curve
            # Data phase 0.5 is not always at the midpoint (not at cadence: len(y)/2),
            # so we need to roll the model to match the model so that its mid-transit
            # is at phase=0.5
            fill_factor = calculate_fill_factor(self.t)
            fill_half = 1 - ((1 - fill_factor) * 0.5)
            stretch = calculate_stretch(self.t, period, transit_times)
            internal_samples = (
                int(len(self.y) / len(transit_times))
            ) * tls_constants.OVERSAMPLE_MODEL_LIGHT_CURVE

            # Folded model flux
            model_folded_model = fractional_transit(
                duration=duration * maxwidth_in_samples * fill_half,
                maxwidth=maxwidth_in_samples / stretch,
                depth=1 - depth,
                samples=int(len(self.t / len(transit_times))),
                per=self.per,
                rp=self.rp,
                a=self.a,
                inc=self.inc,
                ecc=self.ecc,
                w=self.w,
                u=self.u,
                limb_dark=self.limb_dark,
            )
            # Full unfolded light curve model
            model_transit_single = fractional_transit(
                duration=(duration * maxwidth_in_samples),
                maxwidth=maxwidth_in_samples / stretch,
                depth=1 - depth,
                samples=internal_samples,
                per=self.per,
                rp=self.rp,
                a=self.a,
                inc=self.inc,
                ecc=self.ecc,
                w=self.w,
                u=self.u,
                limb_dark=self.limb_dark,
            )
            model_lightcurve_model, model_lightcurve_time = model_lightcurve(
                transit_times, period, self.t, model_transit_single
            )
            depth_mean_odd, depth_mean_even, depth_mean_odd_std, depth_mean_even_std, all_flux_intransit_odd, all_flux_intransit_even, per_transit_count, transit_depths, transit_depths_uncertainties = intransit_stats(
                self.t, self.y, transit_times, transit_duration_in_days
            )
            all_flux_intransit = numpy.concatenate(
                [all_flux_intransit_odd, all_flux_intransit_even]
            )
            snr_per_transit, snr_pink_per_transit = snr_stats(
                t=self.t,
                y=self.y,
                period=period,
                duration=duration,
                T0=T0,
                transit_times=transit_times,
                transit_duration_in_days=transit_duration_in_days,
                per_transit_count=per_transit_count,
            )
            intransit = transit_mask(self.t, period, 2 * duration, T0)
            flux_ootr = self.y[~intransit]
            depth_mean = numpy.mean(all_flux_intransit)
            depth_mean_std = numpy.std(all_flux_intransit) / numpy.sum(
                per_transit_count
            ) ** (0.5)
            snr = ((1 - depth_mean) / numpy.std(flux_ootr)) * len(
                all_flux_intransit
            ) ** (0.5)
            rp_rs = rp_rs_from_depth(depth=1 - depth, law=self.limb_dark, params=self.u)

            if len(all_flux_intransit_odd) > 0:
                depth_mean_odd = numpy.mean(all_flux_intransit_odd)
                depth_mean_odd_std = numpy.std(all_flux_intransit_odd) / numpy.sum(
                    len(all_flux_intransit_odd)
                ) ** (0.5)
            else:
                depth_mean_odd = numpy.nan
                depth_mean_odd_std = numpy.nan

            if len(all_flux_intransit_even) > 0:
                depth_mean_even = numpy.mean(all_flux_intransit_even)
                depth_mean_even_std = numpy.std(all_flux_intransit_even) / numpy.sum(
                    len(all_flux_intransit_even)
                ) ** (0.5)
            else:
                depth_mean_even = numpy.nan
                depth_mean_even_std = numpy.nan

            in_transit_count, after_transit_count, before_transit_count = count_stats(
                self.t, self.y, transit_times, transit_duration_in_days
            )

            # Odd even mismatch in standard deviations
            odd_even_difference = abs(depth_mean_odd - depth_mean_even)
            odd_even_std_sum = depth_mean_odd_std + depth_mean_even_std
            odd_even_mismatch = odd_even_difference / odd_even_std_sum

            transit_count = len(transit_times)
            empty_transit_count = numpy.count_nonzero(per_transit_count == 0)
            distinct_transit_count = transit_count - empty_transit_count

            duration = transit_duration_in_days

            if empty_transit_count / transit_count >= 0.33:
                text = (
                    str(empty_transit_count)
                    + " of "
                    + str(transit_count)
                    + " transits without data. The true period may be twice the given period."
                )
                warnings.warn(text)

        return transitleastsquaresresults(
            SDE,
            SDE_raw,
            chi2_min,
            chi2red_min,
            period,
            period_uncertainty(test_statistic_periods, power),
            T0,
            duration,
            depth,
            (depth_mean, depth_mean_std),
            (depth_mean_even, depth_mean_even_std),
            (depth_mean_odd, depth_mean_odd_std),
            transit_depths,
            transit_depths_uncertainties,
            rp_rs,
            snr,
            snr_per_transit,
            snr_pink_per_transit,
            odd_even_mismatch,
            transit_times,
            per_transit_count,
            transit_count,
            distinct_transit_count,
            empty_transit_count,
            FAP(SDE),
            in_transit_count,
            after_transit_count,
            before_transit_count,
            test_statistic_periods,
            power,
            power_raw,
            SR,
            chi2,
            chi2red,
            model_lightcurve_time,
            model_lightcurve_model,
            model_folded_phase,
            folded_y,
            folded_dy,
            folded_phase,
            model_folded_model,
        )
