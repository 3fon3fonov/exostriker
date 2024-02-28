# -*- coding: utf-8 -*-

__all__ = ["corner", "hist2d", "quantile", "overplot_lines", "overplot_points"]

from .core import hist2d, overplot_lines, overplot_points, quantile
from .corner import corner
#from .corner_version import version as __version__
from .corner_version import __version__  # noqa
