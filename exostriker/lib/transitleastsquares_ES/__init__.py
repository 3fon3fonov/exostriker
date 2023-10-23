#  Optimized algorithm to search for transits of small extrasolar planets
#                                                                            /
#       ,        AUTHORS                                                   O/
#    \  :  /     Michael Hippke (1) [michael@hippke.org]                /\/|
# `. __/ \__ .'  Rene' Heller (2) [heller@mps.mpg.de]                      |
# _ _\     /_ _  _________________________________________________________/ \_
#    /_   _\
#  .'  \ /  `.   (1) Sonneberg Observatory, Sternwartestr. 32, Sonneberg
#    /  :  \     (2) Max Planck Institute for Solar System Research,
#       '            Justus-von-Liebig-Weg 3, 37077 G\"ottingen, Germany

from __future__ import division, print_function
from transitleastsquares_ES.main import transitleastsquares_ES
from transitleastsquares_ES.helpers import cleaned_array, resample, transit_mask
from transitleastsquares_ES.grid import duration_grid, period_grid
from transitleastsquares_ES.stats import FAP
from transitleastsquares_ES.catalog import catalog_info
from transitleastsquares_ES.core import fold
