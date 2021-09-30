#! /usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (c) 2012 Sebastian Schröter, Stefan Czesla, and Mathias Zechmeister

# The MIT License (MIT)

# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# Note : The software is also available as part of the PyAstronomy package.
#        See: http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/index.html
#             https://github.com/sczesla/PyAstronomy/blob/master/src/pyTiming/pyPeriod/gls.py

from __future__ import print_function, division
import numpy as np
from numpy import sum, pi, cos, sin, arctan2, exp, log, sqrt,\
                  dot, arange

__version__ = '2018-10-17'
__author__ = 'Mathias Zechmeister, Stefan Czesla'

class Gls:
    """
    Compute the Generalized Lomb-Scargle (GLS) periodogram.
    The *Gls* class computes the error-weighted Lomb-Scargle periodogram as
    developed by [ZK09]_ using various possible normalizations.
    The constructor of *Gls* takes a *TimeSeries* instance (i.e., a light curve)
    as first argument. The constructor allows to pass keywords to adjust the
    `freq` array, which will be used to calculate the periodogram.
    The main result of the calculation, i.e., the power, are stored in the
    class property `power`.
    Parameters
    ----------
    lc : tuple or list or str or TimeSeries object
        The light curve data either in the form of a TimeSeries object (or any
        object providing the attributes time, flux, and error) or a tuple or list
        or a filename as str providing time as first element, flux as second
        element, and optionally, the error as third element. 
    fbeg, fend : float, optional
        The beginning and end frequencies for the periodogram
        (inverse units of time axis).
    Pbeg, Pend : float, optional
        The beginning and end periods for the periodogram
        (same units as for time axis).
    ofac : int
        Oversampling factor of frequency grid (default=10).
    hifac : float
        Maximum frequency `freq` = `hifac` * (average Nyquist frequency)
        (default=1).
    freq : array, optional
        Contains the frequencies at which to calculate the periodogram.
        If given, fast and verbose option are not available.
        If not given, a frequency array will be automatically generated.
    norm : string, optional
        The normalization; either of "ZK", "Scargle", "HorneBaliunas", "Cumming", "wrms", "chisq".
        The default is unity ("ZK").
    ls : boolean, optional
        If True, the conventional Lomb-Scargle periodogram will be computed
        (default is False).
    fast : boolean, optional
        If True, recursive relations for trigonometric functions will be used
        leading to faster evaluation (default is False).
    verbose : boolean, optional
        Set True to obtain some statistical output (default is False).
    Attributes
    ----------
    power : array
        The normalized power of the GLS.
    freq : array
        The frequency array.
    ofac : int
        The oversampling factor of frequency grid.
    hifac : float
        The maximum frequency.
    t : array
        The abscissa data values.
    y : array
        The ordinate data values.
    e_y : array
        The errors of the data values.
    norm : string, {'ZK', 'Scargle', 'HorneBaliunas', 'Cumming', 'wrms', 'chisq'}
        The used normalization.
    Examples
    --------
    Create 1000 unevenly sampled data points with frequency=0.1,
    measurement error and Gaussian noise
    >>> time = np.random.uniform(54000., 56000., 1000)
    >>> flux = 0.15 * np.sin(2. * np.pi * time / 10.)
    Add some noise
    >>> error = 0.5 * np.ones(time.size)
    >>> flux += np.random.normal(0, error)
    Compute the full error-weighted Lomb-Periodogram
    in 'ZK' normalization and calculate the significance
    of the maximum peak.
    >>> gls = Gls((time, flux, error), verbose=True)
    >>> maxPower = gls.pmax
    >>> print("GLS maximum power: ", maxPower)
    >>> print("GLS statistics of maximum power peak: ", gls.stats(maxPower))
    >>> gls.plot(block=True)
    """
    # Available normalizations
    norms = ['ZK', 'Scargle', 'HorneBaliunas', 'Cumming', 'wrms', 'chisq', 'lnL', 'dlnL']

    def __init__(self, lc, fbeg=None, fend=None, Pbeg=None, Pend=None, ofac=10, hifac=1, freq=None, norm="ZK", ls=False, fast=False, verbose=False, **kwargs):

        self._freq = freq
        self.fbeg = fbeg
        self.fend = fend
        self.Pbeg = Pbeg
        self.Pend = Pend
        self.ofac = ofac
        self.hifac = hifac
        self.ls = ls
        self.norm = norm
        self.fast = fast
        self.label = {'title': 'Generalized Lomb Periodogram',
                      'xlabel': 'Frequency'}

        self._normcheck(norm)

        self._assignTimeSeries(lc)
        self._buildFreq()
        self._calcPeriodogram()
        self.pnorm(norm)
        self._peakPeriodogram()

        self._calcWF()

        # Output statistics
        if verbose:
            self.info()

    def _assignTimeSeries(self, lc):
        """
        A container class that holds the observed light curve.
        Parameters
        ----------
        time : array
            The time array.
        flux : array
            The observed flux/data.
        error : array, optional
            The error of the data values.
        """
        if isinstance(lc, str):
            # A data file has been given.
            try:
               lc = np.loadtxt(lc, unpack=True)[0:3]
            except Exception as e:
               print("An error occurred while trying to read data file:")
               print("  " + str(e))

        if isinstance(lc, (tuple, list, np.ndarray)):
            # t, y[, e_y] were given as list or tuple.
            if len(lc) in (2, 3):
                self.t = np.ravel(lc[0])
                self.y = np.ravel(lc[1])
                self.e_y = None
                if len(lc) == 3 and lc[2] is not None:
                    # Error has been specified.
                    self.e_y = np.ravel(lc[2])
            else:
                raise(ValueError("lc is a list or tuple with " + str(len(lc)) + " elements. Needs to have 2 or 3 elements." + \
                                   " solution=Use 2 or 3 elements (t, y[, e_y]) or an instance of TimeSeries"))
        else:
            # Assume lc is an instance of TimeSeries.
            self.t, self.y, self.e_y = lc.time, lc.flux, lc.error

        self.th = self.t - self.t.min()
        self.tbase = self.th.max()
        self.N = len(self.y)

        # Re-check array length compatibility.
        if (len(self.th) != self.N) or ((self.e_y is not None) and (len(self.e_y) != self.N)):
            raise(ValueError("Incompatible dimensions of input data arrays (time and flux [and error]). Current shapes are: " + \
                             ', '.join(str(np.shape(x)) for x in (self.t, self.y, self.e_y))))

    def _buildFreq(self):
        """
        Build frequency array (`freq` attribute).
        Attributes
        ----------
        fnyq : float
            Half of the average sampling frequency of the time series.
        """
        self.fstep = 1 / self.tbase / self.ofac   # frequency sampling depends on the time span, default for start frequency
        self.fnyq = 0.5 / self.tbase * self.N     # Nyquist frequency
        self.f = self._freq

        if self.f is None:
            # Build frequency array if not present.
            if self.fbeg is None:
                self.fbeg = self.fstep if self.Pend is None else 1 / self.Pend
            if self.fend is None:
                self.fend = self.fnyq * self.hifac if self.Pbeg is None else 1 / self.Pbeg

            if self.fend <= self.fbeg:
                raise(ValueError("fend is smaller than (or equal to) fbeg but it must be larger." + \
                               "Choose fbeg and fend so that fend > fbeg."))

            self.f = arange(self.fbeg, self.fend, self.fstep)
        elif self.fast:
            raise(ValueError("freq and fast cannot be used together."))

        self.freq = self.f   # alias name
        self.nf = len(self.f)

        # An ad-hoc estimate of the number of independent frequencies (Eq. (24) in ZK_09).
        self.M = (self.fend-self.fbeg) * self.tbase

    def _calcPeriodogram(self):

        if self.e_y is None:
            w = np.ones(self.N)
        else:
            w = 1 / (self.e_y * self.e_y)
        self.wsum = w.sum()
        w /= self.wsum

        self._Y = dot(w, self.y)       # Eq. (7)
        wy = self.y - self._Y          # Subtract weighted mean
        self._YY = dot(w, wy**2)       # Eq. (10), weighted variance with offset only
        wy *= w                        # attach errors

        C, S, YC, YS, CC, CS = np.zeros((6, self.nf))

        if self.fast:
            # Prepare trigonometric recurrences.
            eid = exp(2j * pi * self.fstep * self.th)  # cos(dx)+i sin(dx)

        for k, omega in enumerate(2.*pi*self.f):
            # Circular frequencies.
            if self.fast:
                if k % 1000 == 0:
                    # init/refresh recurrences to stop error propagation
                    eix = exp(1j * omega * self.th)  # exp(ix) = cos(x) + i*sin(x)
                cosx = eix.real
                sinx = eix.imag
            else:
                x = omega * self.th
                cosx = cos(x)
                sinx = sin(x)

            C[k] = dot(w, cosx)         # Eq. (8)
            S[k] = dot(w, sinx)         # Eq. (9)

            YC[k] = dot(wy, cosx)       # Eq. (11)
            YS[k] = dot(wy, sinx)       # Eq. (12)
            wcosx = w * cosx
            CC[k] = dot(wcosx, cosx)    # Eq. (13)
            CS[k] = dot(wcosx, sinx)    # Eq. (15)

            if self.fast:
               eix *= eid              # increase freq for next loop

        SS = 1. - CC
        if not self.ls:
            CC -= C * C            # Eq. (13)
            SS -= S * S            # Eq. (14)
            CS -= C * S            # Eq. (15)
        D = CC*SS - CS*CS          # Eq. (6)

        self._a = (YC*SS-YS*CS) / D
        self._b = (YS*CC-YC*CS) / D
        self._off = -self._a*C - self._b*S

        # power
        self.p = (SS*YC*YC + CC*YS*YS - 2.*CS*YC*YS) / (self._YY*D)   # Eq. (5) in ZK09


    def _calcWF(self):
        

        ######################## DFT (Window) ##############################

        WF_power = []
        for omi in 2.*pi*self.f: 
            phase = self.th * omi
            WC = np.sum(np.cos(phase))
            WS = np.sum(np.sin(phase))
            WF_power.append((WC**2 + WS**2)/self.N**2) 

        self.WF_power = np.array(WF_power)
        self.WF_omega = self.f


    def _normcheck(self, norm):
        """
        Check normalization
        Parameters
        ----------
        norm : string
            Normalization string
        """
        if norm not in self.norms:
            raise(ValueError("Unknown norm: " + str(norm) + ". " + \
                "Use either of " + ', '.join(self.norms)))

    def pnorm(self, norm="ZK"):
        """
        Assign or modify normalization (can be done afterwards).
        Parameters
        ----------
        norm : string, optional
            The normalization to be used (default is 'ZK').
        Examples
        --------
        >>> gls.pnorm('wrms')
        """
        self._normcheck(norm)
        self.norm = norm
        p = self.p
        power = p   # default ZK
        self.label["ylabel"] = "Power ("+norm+")"

        if norm == "Scargle":
            popvar = input('pyTiming::gls - Input a priori known population variance:')
            power = p / float(popvar)
        elif norm == "HorneBaliunas":
            power = (self.N-1)/2. * p
        elif norm == "Cumming":
            power = (self.N-3)/2. * p / (1.-self.p.max())
        elif norm == "chisq":
            power = self._YY *self.wsum * (1.-p)
            self.label["ylabel"] = "chisq"
        elif norm == "wrms":
            power = sqrt(self._YY*(1.-p))
            self.label["ylabel"] = "wrms"
        elif norm == "lnL":
            chi2 = self._YY *self.wsum * (1.-p)
            power = -0.5*chi2 - 0.5*np.sum(np.log(2*np.pi * self.e_y**2))
            self.label["ylabel"] = "lnL"
        elif norm == "dlnL":
            # dlnL = lnL - lnL0 = -0.5 chi^2 + 0.5 chi0^2 = 0.5 (chi0^2 - chi^2) = 0.5 chi0^2 p
            power = 0.5 * self._YY * self.wsum * p
            self.label["ylabel"] = "dlnL"

        self.power = power

    def _peakPeriodogram(self):
        """
        Analyze the highest periodogram peak.
        """
        # Index of maximum power
        k = self.p.argmax()
        # Maximum power
        self.pmax = pmax = self.p[k]
        self.rms = rms = sqrt(self._YY*(1.-pmax))

        # Statistics of highest peak
        self.hpstat = self.best = p = {}   # alias name for best and hpstat

        # Best parameters
        p["f"] = fbest = self.f[k]
        p["P"] = 1. / fbest
        p["amp"] = amp = sqrt(self._a[k]**2 + self._b[k]**2)
        p["ph"] = ph = arctan2(self._a[k], self._b[k]) / (2.*pi)
        p["T0"]  = self.t.min() - ph/fbest
        p["offset"] = self._off[k] + self._Y            # Re-add the mean.

        # Error estimates
        p["e_amp"] = sqrt(2./self.N) * rms
        p["e_ph"] = e_ph = sqrt(2./self.N) * rms/amp/(2.*pi)
        p["e_T0"] = e_ph / fbest
        p["e_offset"] = sqrt(1./self.N) * rms

        # Get the curvature in the power peak by fitting a parabola y=aa*x^2
        if 1 < k < self.nf-2:
            # Shift the parabola origin to power peak
            xh = (self.f[k-1:k+2] - self.f[k])**2
            yh = self.p[k-1:k+2] - pmax
            # Calculate the curvature (final equation from least square)
            aa = dot(yh, xh) / dot(xh, xh)
            p["e_f"] = e_f = sqrt(-2./self.N / aa * (1.-self.pmax))
            p["e_P"] = e_f / fbest**2
        else:
            p["e_f"] = np.nan
            p["e_P"] = np.nan
            print("WARNING: Highest peak is at the edge of the frequency range.\n"\
                  "No output of frequency error.\n"\
                  "Increase frequency range to sample the peak maximum.")

    def sinmod(self, t=None):
        """
        Calculate best-fit sine curve.
        The parameters of the best-fit sine curve can be accessed via
        the dictionary attribute `best`. Specifically, "amp" holds the
        amplitude, "fbest" the best-fit frequency, "T0" the reference time
        (i.e., time of zero phase), and "offset" holds the additive offset
        of the sine wave. 
        Parameters
        ----------
        t : array
            Time array at which to calculate the sine.
            If None, the time of the data are used.
        Returns
        -------
        Sine curve : array
            The best-fit sine curve (i.e., that for which the
            power is maximal).
        """
        if t is None:
            t = self.t

        try:
            p = self.best
            return p["amp"] * sin(2*np.pi*p["f"]*(t-p["T0"])) + p["offset"]
        except Exception as e:
            print("Failed to calcuate best-fit sine curve.")
            raise(e)

    def info(self, stdout=True):
        """
        Prints some basic statistical output screen.
        """
        lines = ("Generalized LS - statistical output",
           "-----------------------------------",
           "Number of input points:     %6d" % self.N,
           "Weighted mean of dataset:   %f"  % self._Y,
           "Weighted rms of dataset:    %f"  % sqrt(self._YY),
           "Time base:                  %f"  % self.tbase,
           "Number of frequency points: %6d" % self.nf,
           "",
           "Maximum power p [%s]: %f" % (self.norm, self.power.max()),
           "RMS of residuals:     %f" % self.rms)
        if self.e_y is not None:
           lines += "  Mean weighted internal error:  %f" % (sqrt(self.N/sum(1./self.e_y**2))),
        lines += (
           "Best sine frequency:  {f:f} +/- {e_f:f}",
           "Best sine period:     {P:f} +/- {e_P:f}",
           "Amplitude:            {amp:f} +/- {e_amp:f}",
           #"Phase (ph):          {ph:f} +/- {e_ph:f}",
           #"Phase (T0):          {T0:f} +/- {e_T0:f}",
           #"Offset:              {offset:f} +/- {e_offset:f}",
           "-----------------------------------")
        text = "\n".join(lines).format(**self.best)
        if stdout:
           print(text)
        else:
           return text

    def plot(self, block=False, period=False):
        """
        Create a plot.
        """
        try:
            import matplotlib
            import matplotlib.pylab as mpl
        except ImportError:
            raise(ImportError("Could not import matplotlib.pylab."))

        fbest, T0 = self.best["f"], self.best["T0"]

        fig = mpl.figure()
        fig.canvas.set_window_title('GLS periodogram')
        fig.subplots_adjust(hspace=0.05, wspace=0.04, right=0.97, bottom=0.09, top=0.84)
        fs = 10   # fontsize

        # Periodogram
        plt = fig.add_subplot(3, 1, 1)
        plt.tick_params(direction='in')
        if period:
           plt.set_xscale("log")
           plt.set_xlabel("Period P")
        else:
           plt.set_xlabel("Frequency f")

        plt.set_ylabel(self.label["ylabel"])
        plt.plot(1/self.f if period else self.f, self.power, 'b-', linewidth=.5)
        # mark the highest peak
        plt.plot(1/fbest if period else fbest, self.power[self.p.argmax()], 'r.', label="1/f = %f" % (1/fbest))
        plt.legend(numpoints=1, fontsize=fs, frameon=False)

        x2tics = 1/np.array([0.5, 1, 2, 3, 5, 10, 20., 100])
        mx2tics = 1/np.array([0.75, 1.5, 2.5, 4, 15, 40, 60., 80, 100])
        def tick_function(X):
           return ["%g" % (1/z) for z in X]

        plt.tick_params(direction='in', which='both', top=True, right=True)
        plt.minorticks_on()
        plt.autoscale(enable=True, axis='x', tight=True)
        if not period:
           ax2 = plt.twiny()
           ax2.tick_params(direction='in', which='both')
           ax2.format_coord = lambda x,y: "x=%g, x2=%g, y=%g"% (x, 1/x, y)
           ax2.set_xticks(x2tics)
           ax2.set_xticks(mx2tics, minor=True)
           ax2.set_xticklabels(tick_function(x2tics))
           ax2.set_xlim(plt.get_xlim())
           ax2.set_xlabel("Period")
           plt.tick_params(top=False)

        # Data and model
        col = mpl.cm.rainbow(mpl.Normalize()(self.t))
        def plot_ecol(plt, x, y):
           # script for scatter plot with errorbars and time color-coded
           datstyle = dict(color=col, marker='.', edgecolor='k', linewidth=0.5, zorder=2)
           if self.e_y is not None:
              errstyle = dict(yerr=self.e_y, marker='', ls='', elinewidth=0.5)
              if matplotlib.__version__ < '2.' :
                 errstyle['capsize'] = 0.
                 datstyle['s'] = 8**2   # requires square size !?
              else:
                 errstyle['ecolor'] = col
              _, _, (c,) = plt.errorbar(x, y, **errstyle)
              if matplotlib.__version__ < '2.':
                 c.set_color(col)
           plt.scatter(x, y, **datstyle)

        def phase(t):
           #return (t-T0)*fbest % 1
           return (t-T0) % (1/fbest)

        # Time series
        tt = arange(self.t.min(), self.t.max(), 0.01/fbest)
        ymod = self.sinmod(tt)
        yfit = self.sinmod()
        plt1 = fig.add_subplot(3, 2, 3)
        plt1.set_ylabel("Data")
        mpl.setp(plt1.get_xticklabels(), visible=False)
        plot_ecol(plt1, self.t, self.y)
        plt1.plot(tt, ymod, 'k-', zorder=0)

        # Phase folded data
        tt = arange(T0, T0+1/fbest, 0.01/fbest)
        yy = self.sinmod(tt)
        plt2 = fig.add_subplot(3, 2, 4, sharey=plt1)
        mpl.setp(plt2.get_xticklabels(), visible=False)
        mpl.setp(plt2.get_yticklabels(), visible=False)
        plot_ecol(plt2, phase(self.t), self.y)
        xx = phase(tt)
        ii = np.argsort(xx)
        plt2.plot(xx[ii], yy[ii], 'k-')
        plt2.format_coord = lambda x,y: "x=%g, x2=%g, y=%g"% (x, x*fbest, y)

        # Time serie of residuals
        yres = self.y - yfit
        plt3 = fig.add_subplot(3, 2, 5, sharex=plt1)
        plt3.set_xlabel("Time")
        plt3.set_ylabel("Residuals")
        plot_ecol(plt3, self.t, yres)
        plt3.plot([self.t.min(), self.t.max()], [0,0], 'k-')

        # Phase folded residuals
        plt4 = fig.add_subplot(3, 2, 6, sharex=plt2, sharey=plt3)
        plt4.set_xlabel("Phase")
        mpl.setp(plt4.get_yticklabels(), visible=False)
        plot_ecol(plt4, phase(self.t), yres)
        plt4.plot([0,1/fbest], [0,0], 'k-')
        plt4.format_coord = lambda x,y: "x=%g, x2=%g, y=%g"% (x, x*fbest, y)

        for x in [plt1, plt2, plt3, plt4]:
           x.tick_params(direction='in', which='both', top=True, right=True)
           x.minorticks_on()
           x.autoscale(enable=True, tight=True)

        if hasattr(mpl.get_current_fig_manager(), 'toolbar'):
            # check seems not needed when "TkAgg" is set
            mpl.get_current_fig_manager().toolbar.pan()
        #t = fig.canvas.toolbar
        #mpl.ToggleTool(mpl.wx_ids['Pan'], False)
        if block:
           print("Close the plot to continue.")
        else:
           mpl.ion()

        fig.tight_layout()   # to get the left margin
        marleft = fig.subplotpars.left * fig.get_figwidth() * fig.dpi / fs
        def tighter():
           # keep margin tight when resizing
           xdpi = fs / (fig.get_figwidth() * fig.dpi)
           ydpi = fs / (fig.get_figheight() * fig.dpi)
           fig.subplots_adjust(bottom=4.*ydpi, top=1-8.*ydpi, right=1-1*xdpi, wspace=4*xdpi, hspace=4*ydpi, left=marleft*xdpi)
           if matplotlib.__version__ < '2.':
              ax2.set_position(plt.get_position().translated(0,4*ydpi))
           plt.set_position(plt.get_position().translated(0,4*ydpi))

        #fig.canvas.mpl_connect("resize_event", lambda _: (fig.tight_layout()))
        fig.canvas.mpl_connect("resize_event", lambda _: (tighter()))
        mpl.show()
        # mpl.show(block=block) # unexpected keyword argument 'block' in older matplotlib
        return mpl

    def prob(self, Pn):
        """
        Probability of obtaining the given power.
        Calculate the probability to obtain a power higher than
        `Pn` from the noise, which is assumed to be Gaussian.
        .. note:: Normalization
          (see [ZK09]_ for further details).
          - `Scargle`:
          .. math::
            exp(-Pn)
          - `HorneBaliunas`:
          .. math::
            \\left(1 - 2 \\times \\frac{Pn}{N-1} \\right)^{(N-3)/2}
          - `Cumming`:
          .. math::
            \\left(1+2\\times \\frac{Pn}{N-3}\\right)^{-(N-3)/2}
        Parameters
        ----------
        Pn : float
            Power threshold.
        Returns
        -------
        Probability : float
            The probability to obtain a power equal or
            higher than the threshold from the noise.
        """
        self._normcheck(self.norm)
        if self.norm == "ZK": return (1.-Pn)**((self.N-3.)/2.)
        if self.norm == "Scargle": return exp(-Pn)
        if self.norm == "HorneBaliunas": return (1.-2.*Pn/(self.N-1.))**((self.N-3.)/2.)
        if self.norm == "Cumming": return (1.+2.*Pn/(self.N-3.))**(-(self.N-3.)/2.)
        if self.norm == "wrms": return (Pn**2/self._YY)**((self.N-3.)/2.)
        if self.norm == "chisq": return (Pn/self._YY/self.wsum)**((self.N-3.)/2.)

    def probInv(self, Prob):
        """
        Calculate minimum power for given probability.
        This function is the inverse of `Prob(Pn)`.
        Returns the minimum power for a given probability threshold `Prob`.
        Parameters
        ----------
        Prob : float
            Probability threshold.
        Returns
        -------
        Power threshold : float
            The minimum power for the given false-alarm probability threshold.
        """
        self._normcheck(self.norm)
        if self.norm == "ZK": return 1.-Prob**(2./(self.N-3.))
        if self.norm == "Scargle": return -log(Prob)
        if self.norm == "HorneBaliunas": return (self.N-1) / 2. * (1.-Prob**(2./(self.N-3)))
        if self.norm == "Cumming": return (self.N-3) / 2. * (Prob**(-2./(self.N-3.))-1.)
        if self.norm == "wrms": return sqrt(self._YY * Prob**(2./(self.N-3.)))
        if self.norm == "chisq": return self._YY * self.wsum * Prob**(2./(self.N-3.))

    def FAP(self, Pn=None):
        """
        Obtain the false-alarm probability (FAP).
        The FAP denotes the probability that at least one out of M independent
        power values in a prescribed search band of a power spectrum computed
        from a white-noise time series is as large as or larger than the
        threshold, `Pn`. It is assessed through
        .. math:: FAP(Pn) = 1 - (1-Prob(P>Pn))^M \\; ,
        where "Prob(P>Pn)" depends on the type of periodogram and normalization
        and is calculated by using the *prob* method; *M* is the number of
        independent power values and is computed internally.
        Parameters
        ----------
        Pn : float
            Power threshold. If None, the highest periodogram peak is used.
        Returns
        -------
        FAP : float
            False alarm probability.
        """
        if Pn is None:
           Pn = self.pmax
        prob = self.M * self.prob(Pn)
        if prob > 0.01:
           return 1. - (1.-self.prob(Pn))**self.M
        return prob

    def powerLevel(self, FAPlevel):
        """
        Power threshold for FAP level.
        Parameters
        ----------
        FAPlevel : float or array_like
              "False Alarm Probability" threshold
        Returns
        -------
        Threshold : float or array
            The power threshold pertaining to a specified false-alarm
            probability (FAP). Powers exceeding this threshold have FAPs
            smaller than FAPlevel.
        """
        Prob = 1. - (1.-FAPlevel)**(1./self.M)
        return self.probInv(Prob)

    def stats(self, Pn):
        """
        Obtain basic statistics for power threshold.
        Parameters
        ----------
        Pn : float
            Power threshold.
        Returns
        -------
        Statistics : dictionary
            A dictionary containing {'Pn': *Pn*, 'Prob': *Prob(Pn)* ,
            'FAP': *FAP(Pn)*} for the specified power threshold, *Pn*.
        """
        return {'Pn': Pn, 'Prob': self.prob(Pn), 'FAP': self.FAP(Pn)}

    def toFile(self, ofile, header=True):
        """
        Write periodogram to file.
        The output file is a standard text file with two columns,
        viz., frequency and power.
        Parameters
        ----------
        ofile : string
            Name of the output file.
        """
        with open(ofile, 'w') as f:
            if header:
               f.write("# Generalized Lomb-Scargle periodogram\n")
               f.write("# Parameters:\n")
               f.write("#    Data file: %s\n" % self.df)
               f.write("#    ofac     : %s\n" % self.ofac)
               f.write("#    norm     : %s\n" % self.norm)
               f.write("# 1) Frequency, 2) Normalized power\n")
            for line in zip(self.f, self.power):
               f.write("%f  %f\n" % line)

        print("Results have been written to file: ", ofile)


def example():
    # Run the example in the Gls class.
    print("--- EXAMPLE CALCULATION ---")
    import doctest
    exec(doctest.script_from_examples(Gls.__doc__))
    print("----------------------------------------------------")


if __name__ == "__main__":

  import argparse

  parser = argparse.ArgumentParser(description='Generalized Lomb-Scargle periodogram.', add_help=False)
  argadd = parser.add_argument   # function short cut
  argadd('df', nargs='?',
                 help='Data file (three columns: time, data, error). If not specified example will be shown.')
  argadd('-fbeg', '--fbeg', type=float, help="Starting frequency for periodogram.")
  argadd('-fend', '--fend', type=float, help="Stopping frequency for periodogram.")
  argadd('-Pbeg', '--Pbeg', type=float, help="Starting period for periodogram.")
  argadd('-Pend', '--Pend', type=float, help="Stopping period for periodogram.")
  argadd('-ofac', '--ofac', type=float, help="Oversampling factor (default=10).", default=10)
  argadd('-hifac', '--hifac', type=float, help="Maximum frequency (default=1).", default=1)
  argadd('-fast', '--fast', help="Use trigonometric recurrences.", action='store_true')
  argadd('-norm', '--norm', help="The normalization (default=ZK).", choices=Gls.norms, default='ZK')
  argadd('-ofile', '--ofile', type=str, help="Output file for results.")
  argadd('-noplot', '--noplot', help="Suppress plots.", dest='plot', action='store_false')
  argadd('-nostat', '--nostat', help="Switch off statistical output on screen.", dest='verbose',
                 action='store_false')
  argadd('-?', '-h', '-help', '--help', help='Show this help message and exit.', action='help')

  args = vars(parser.parse_args())
  df = args.pop('df')
  ofile = args.pop('ofile')
  plot = args.pop('plot')

  if df is None:
    # No data file given. Show example:
    example()
    print("Available options:")
    parser.print_help()
    exit(0)

  gls = Gls(df, **args)

  if plot:
     gls.plot(block=True)

  if ofile:
     gls.df = df
     gls.toFile(ofile)
