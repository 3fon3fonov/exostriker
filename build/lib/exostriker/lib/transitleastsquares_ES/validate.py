from __future__ import division, print_function
import numpy
import warnings
import multiprocessing
from transitleastsquares_ES.helpers import cleaned_array, impact_to_inclination
import transitleastsquares_ES.tls_constants as tls_constants


def validate_inputs(t, y, dy):
    """Check the consistency of the inputs"""

    # Clean array
    if dy is None:
        t, y = cleaned_array(t, y)
    else:
        t, y, dy = cleaned_array(t, y, dy)
        # Normalize dy to act as weights in least squares calculatio
        dy = dy / numpy.mean(dy)

    duration = max(t) - min(t)
    if duration <= 0:
        raise ValueError("Time duration must positive")
    if numpy.size(y) < 3 or numpy.size(t) < 3:
        raise ValueError("Too few values in data set")
    if numpy.mean(y) > 1.01 or numpy.mean(y) < 0.99:
        text = (
            "Warning: The mean flux should be normalized to 1"
            + ", but it was found to be "
            + str(numpy.mean(y))
        )
        warnings.warn(text)

    if min(y) < 0:
        raise ValueError("Flux values must be positive")
    if max(y) >= float("inf"):
        raise ValueError("Flux values must be finite")

    # If no dy is given, create it with the standard deviation of the flux
    if dy is None:
        dy = numpy.full(len(y), numpy.std(y))
    if numpy.size(t) != numpy.size(y) or numpy.size(t) != numpy.size(dy):
        raise ValueError("Arrays (t, y, dy) must be of the same dimensions")
    if t.ndim != 1:  # Size identity ensures dimensional identity
        raise ValueError("Inputs (t, y, dy) must be 1-dimensional")

    return t, y, dy


def validate_args(self, kwargs):

    self.verbose = kwargs.get("verbose", tls_constants.VERBOSE)

    # Warn user if unknown parameters
    for key, value in kwargs.items():
        if key not in tls_constants.VALID_PARAMETERS:
            text = "Ignoring unknown parameter: " + str(key)
            warnings.warn(text)

    """Validate **kwargs and set to defaults where missing"""
    self.show_progress_bar = kwargs.get("show_progress_bar", True)
    self.transit_depth_min = kwargs.get(
        "transit_depth_min", tls_constants.TRANSIT_DEPTH_MIN
    )
    self.R_star = kwargs.get("R_star", tls_constants.R_STAR)
    self.M_star = kwargs.get("M_star", tls_constants.M_STAR)
    self.oversampling_factor = kwargs.get(
        "oversampling_factor", tls_constants.OVERSAMPLING_FACTOR
    )
    self.period_max = kwargs.get("period_max", float("inf"))
    self.period_min = kwargs.get("period_min", 0)
    self.n_transits_min = kwargs.get("n_transits_min", tls_constants.N_TRANSITS_MIN)

    self.R_star_min = kwargs.get("R_star_min", tls_constants.R_STAR_MIN)
    self.R_star_max = kwargs.get("R_star_max", tls_constants.R_STAR_MAX)
    self.M_star_min = kwargs.get("M_star_min", tls_constants.M_STAR_MIN)
    self.M_star_max = kwargs.get("M_star_max", tls_constants.M_STAR_MAX)
    self.duration_grid_step = kwargs.get(
        "duration_grid_step", tls_constants.DURATION_GRID_STEP
    )

    self.use_threads = kwargs.get("use_threads", multiprocessing.cpu_count())

    self.per = kwargs.get("per", tls_constants.DEFAULT_PERIOD)
    self.rp = kwargs.get("rp", tls_constants.DEFAULT_RP)
    self.a = kwargs.get("a", tls_constants.DEFAULT_A)

    self.T0_fit_margin = kwargs.get("T0_fit_margin", tls_constants.T0_FIT_MARGIN)

    # If an impact parameter is given, it overrules the supplied inclination
    if "b" in kwargs:
        self.b = kwargs.get("b")
        self.inc = impact_to_inclination(b=self.b, semimajor_axis=self.a)
    else:
        self.inc = kwargs.get("inc", tls_constants.DEFAULT_INC)

    self.ecc = kwargs.get("ecc", tls_constants.DEFAULT_ECC)
    self.w = kwargs.get("w", tls_constants.DEFAULT_W)
    self.u = kwargs.get("u", tls_constants.DEFAULT_U)
    self.limb_dark = kwargs.get("limb_dark", tls_constants.DEFAULT_LIMB_DARK)

    self.transit_template = kwargs.get("transit_template", "default")
    if self.transit_template == "default":
        self.per = tls_constants.DEFAULT_PERIOD
        self.rp = tls_constants.DEFAULT_RP
        self.a = tls_constants.DEFAULT_A
        self.inc = tls_constants.DEFAULT_INC

    elif self.transit_template == "grazing":
        self.b = tls_constants.GRAZING_B
        self.inc = impact_to_inclination(b=self.b, semimajor_axis=self.a)

    elif self.transit_template == "box":
        self.per = tls_constants.BOX_PERIOD
        self.rp = tls_constants.BOX_RP
        self.a = tls_constants.BOX_A
        self.b = tls_constants.BOX_B
        self.inc = tls_constants.BOX_INC
        self.u = tls_constants.BOX_U
        self.limb_dark = tls_constants.BOX_LIMB_DARK

    else:
        raise ValueError(
            'Unknown transit_template. Known values: \
            "default", "grazing", "box"'
        )

    """Validations to avoid (garbage in ==> garbage out)"""

    # Stellar radius
    # 0 < R_star < inf
    if self.R_star <= 0 or self.R_star >= float("inf"):
        raise ValueError("R_star must be positive")

    # Assert (0 < R_star_min <= R_star)
    if self.R_star_min > self.R_star:
        raise ValueError("R_star_min <= R_star is required")
    if self.R_star_min <= 0 or self.R_star_min >= float("inf"):
        raise ValueError("R_star_min must be positive")

    # Assert (R_star <= R_star_max < inf)
    if self.R_star_max < self.R_star:
        raise ValueError("R_star_max >= R_star is required")
    if self.R_star_max <= 0 or self.R_star_max >= float("inf"):
        raise ValueError("R_star_max must be positive")

    # Stellar mass
    # Assert (0 < M_star < inf)
    if self.M_star <= 0 or self.M_star >= float("inf"):
        raise ValueError("M_star must be positive")

    # Assert (0 < M_star_min <= M_star)
    if self.M_star_min > self.M_star:
        raise ValueError("M_star_min <= M_star is required")
    if self.M_star_min <= 0 or self.M_star_min >= float("inf"):
        raise ValueError("M_star_min must be positive")

    # Assert (M_star <= M_star_max < inf)
    if self.M_star_max < self.M_star:
        raise ValueError("M_star_max >= M_star required")
    if self.M_star_max <= 0 or self.M_star_max >= float("inf"):
        raise ValueError("M_star_max must be positive")

    # Period grid
    if self.period_min < 0:
        raise ValueError("period_min >= 0 required")
    if self.period_min >= self.period_max:
        raise ValueError("period_min < period_max required")
    if not isinstance(self.n_transits_min, int):
        raise ValueError("n_transits_min must be an integer value")
    if self.n_transits_min < 1:
        raise ValueError("n_transits_min must be an integer value >= 1")

    if not isinstance(self.use_threads, int) or self.use_threads < 1:
        raise ValueError("use_threads must be an integer value >= 1")

    # Assert 0 < T0_fit_margin < 0.1
    if self.T0_fit_margin < 0:
        self.T0_fit_margin = 0
    elif self.T0_fit_margin > 0.1:  # Sensible limit 10% of transit duration
        self.T0_fit_margin = 0.1
    return self, kwargs
