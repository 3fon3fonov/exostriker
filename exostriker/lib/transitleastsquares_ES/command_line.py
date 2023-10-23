from __future__ import division, print_function
from transitleastsquares_ES import transitleastsquares_ES
from transitleastsquares_ES.helpers import cleaned_array
import transitleastsquares_ES.tls_constants as tls_constants
import numpy
import os
import sys
from configparser import ConfigParser

try:
    import argparse
except:
    raise ImportError("Could not import package argparse")


def main():
    pass


print(tls_constants.TLS_VERSION)
parser = argparse.ArgumentParser()
parser.add_argument("lightcurve", help="path to lightcurve file")
parser.add_argument("-o", "--output", help="path to output directory")
parser.add_argument("-c", "--config", help="path to configuration file")
args = parser.parse_args()

# Read config file if possible
use_config_file = False
if args.config is not None:
    try:
        config = ConfigParser()
        config.read(args.config)
        R_star = float(config["Grid"]["R_star"])
        R_star_min = float(config["Grid"]["R_star_min"])
        R_star_max = float(config["Grid"]["R_star_max"])
        M_star = float(config["Grid"]["M_star"])
        M_star_min = float(config["Grid"]["M_star_min"])
        M_star_max = float(config["Grid"]["M_star_max"])
        period_min = float(config["Grid"]["period_min"])
        period_max = float(config["Grid"]["period_max"])
        n_transits_min = int(config["Grid"]["n_transits_min"])
        transit_template = config["Template"]["transit_template"]
        duration_grid_step = float(config["Speed"]["duration_grid_step"])
        transit_depth_min = float(config["Speed"]["transit_depth_min"])
        oversampling_factor = int(config["Speed"]["oversampling_factor"])
        T0_fit_margin = float(config["Speed"]["T0_fit_margin"])
        use_threads = int(config["Speed"]["use_threads"])
        delimiter = int(config["File"]["delimiter"])
        use_config_file = True
        print("Using TLS configuration from config file", args.config)
    except:
        print(
            "Using default values because of broken or missing configuration file",
            args.config,
        )
else:
    print("No config file given. Using default values")

# Load data
if use_config_file:
    data = numpy.genfromtxt(args.lightcurve, delimiter=args.delimiter)
else:
    data = numpy.genfromtxt(args.lightcurve, delimiter=",")

t = data[:, 0]
y = data[:, 1]


# Initiate transitleastsquares_ES model
try:
    dy = data[:, 2]
    model = transitleastsquares_ES(t, y, dy)
except:
    model = transitleastsquares_ES(t, y)

if use_config_file:
    results = model.power(
        R_star=R_star,
        R_star_min=R_star_min,
        R_star_max=R_star_max,
        M_star=M_star,
        M_star_min=M_star_min,
        M_star_max=M_star_max,
        period_min=period_min,
        period_max=period_max,
        n_transits_min=n_transits_min,
        transit_template=transit_template,
        duration_grid_step=duration_grid_step,
        transit_depth_min=transit_depth_min,
        oversampling_factor=oversampling_factor,
        T0_fit_margin=T0_fit_margin,
        use_threads=use_threads,
    )
else:
    results = model.power()


# Save results to CSV files

# Determine path and file names of output files
if args.output is None:
    file_stats = args.lightcurve + "_statistics.csv"
    file_power = args.lightcurve + "_power.csv"
else:
    file_stats = os.path.join(args.output, args.lightcurve + "_statistics.csv")
    file_power = os.path.join(args.output, args.lightcurve + "_power.csv")

# Save
try:
    numpy.savetxt(
        file_power,
        numpy.column_stack(
            [
                list(dict(list(results.items())[25:26]).values())[0],
                list(dict(list(results.items())[26:27]).values())[0],
            ]
        ),
        delimiter=",",
        fmt="%1.6f",
    )
    print("SDE-ogram saved to", file_stats)

    statistics = dict(list(results.items())[0:25])
    numpy.set_printoptions(precision=8, threshold=10e10)
    with open(file_stats, "w") as f:
        for key in statistics.keys():
            f.write("%s %s\n" % (key, statistics[key]))
    print("Statistics saved to", file_stats)
except IOError:
    print("Error saving result file")
