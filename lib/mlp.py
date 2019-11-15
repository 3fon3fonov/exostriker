#! /usr/bin/python
# differentail corrrection for
#
# mlpd /home/astro115/carmenes/data/svn/zero/CARM_VIS/J16167+672S.avc.dat  /home/astro115/carmenes/data/svn/other-rvs/J16167+672S/J16167+672S_hires.dat --fbeg 0.00022


# -*- coding: utf-8 -*-
# Copyright (c) 2012 Sebastian Schrter, Stefan Czesla, and Mathias Zechmeister

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

from __future__ import print_function, division
import numpy as np
from numpy import sum, pi, cos, sin, arctan2, exp, log, sqrt,\
                  dot, argmax, arange
from pause import *
from gplot import *
from scipy import optimize as op

__version__ = '2017-11-08'
__author__ = 'Mathias Zechmeister, Stefan Czesla'

def mod_abc(x, a):
   ''' sine model with multiple data set and offsets
   x - list of tuples with cosx and sinx terms
   '''
   y = [a[0]*cosx + a[1]*sinx + a[2+i] for i,(cosx,sinx) in enumerate(x)]
   return y

def mod_c(x, a):
   ''' sine model with multiple data set and offsets
   x - list of tuples with cosx and sinx terms
   '''
   y = [0*cosx + a[i] for i,cosx in enumerate(x)]
   return y

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
    lc : TimeSeries object or tuple or list
        The light curve data either in the form of a TimeSeries object (or any
        object providing the attributes time, flux, and error) or a tuple or list
        providing time as first element, flux as second element, and optionally,
        the error as third element.
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
        The normalization; either of "lnL", "Scargle", "HorneBaliunas", "Cumming", "wrms", "chisq".
        The default is unity ("lnL").
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
    norm : string, {'lnL', 'Scargle', 'HorneBaliunas', 'Cumming', 'wrms', 'chisq'}
        The used normalization.
    Examples
    --------
    Create 1000 unevenly sampled data points with frequency=0.1,
    measurement error and Gaussian noise
    >>> time = np.random.uniform(54000., 56000., 1000)
    >>> flux = 0.15 * np.sin(2. * np.pi * time / 10.)
    Add some noise
    >>> error = 0.3 * np.ones(time.size)
    >>> flux += np.random.normal(0, error+0.2)
    Compute the full error-weighted Lomb-Periodogram
    in 'lnL' normalization and calculate the significance
    of the maximum peak.
    >>> gls = Gls((time, flux, error), verbose=True)
    >>> maxPower = gls.pmax
    >>> print("GLS maximum power: ", maxPower)
    >>> print("GLS statistics of maximum power peak: ", gls.stats(maxPower))
    >>> gls.plot(block=True)
    """
    # Available normalizations
    norms = ['dlnL', 'lnL', 'Scargle', 'HorneBaliunas', 'Cumming', 'wrms', 'chisq']

    def __init__(self, lc, fbeg=None, fend=None, Pbeg=None, Pend=None, ofac=10, hifac=1, freq=None, norm="dlnL", ls=False, fast=False, verbose=False, **kwargs):

        self.freq = freq
        self.fbeg = fbeg
        self.fend = fend
        self.Pbeg = Pbeg
        self.Pend = Pend
        self.ofac = ofac
        self.hifac = hifac
        self.ls = ls
        self.norm = norm
        self.fast = fast
        self.label = {'title': 'Maximum Likelihood Periodogram',
                      'xlabel': 'Frequency'}
        if "stats" in kwargs:
          print("Warning: 'stats' option is outdated. Please use 'verbose' instead.")
          verbose = kwargs["stats"]

        self._normcheck(norm)

        self._assignTimeSeries(lc)
        self._buildFreq()
        self._calcPeriodogram()
        self.pnorm(norm)
        self._peakPeriodogram()

        # Output statistics
        if verbose:
            self.info()

    def _assignTimeSeries(self, lcs):
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
      self.t = []
      self.y = []
      self.e_y = []
      self.N = 0
      self.Nj = len(lcs)
      for lc in lcs:
         if isinstance(lc, (tuple, list)):
             # t, y[, e_y] were given as list or tuple.
             if len(lc) in (2, 3):
                 t = np.ravel(lc[0])
                 y = np.ravel(lc[1])
                 e_y = None
                 if len(lc) == 3 and lc[2] is not None:
                     # Error has been specified.
                     e_y = np.ravel(lc[2])
             else:
                 raise(ValueError("lc is a list or tuple with " + str(len(lc)) + " elements. Needs to have 2 or 3 elements." + \
                                    " solution=Use 2 or 3 elements (t, y[, e_y]) or an instance of TimeSeries"))
         else:
             # Assume lc is an instance of TimeSeries.
             t, y, e_y = lc.time, lc.flux, lc.error
         self.t += [t]
         self.y += [y]
         self.e_y += [e_y]

         N = len(y)
         self.N += N

         # Re-check array length compatibility
         if (len(t) != N) or ((e_y is not None) and (len(e_y) != N)):
             raise(ValueError("Incompatible dimensions of input data arrays (time and flux [and error]). Current shapes are: " + \
                              ', '.join(str(np.shape(x)) for x in (t, y, e_y))))

      self.data = zip(self.t, self.y, self.e_y)
      self.tmin = min(map(min,self.t))
      self.th = [t - self.tmin for t in self.t]
      self.tbase = max(map(max,self.th))


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
        self.f = self.freq

        if self.freq is None:
            # Build frequency array if not present.
            if self.fbeg is None:
                self.fbeg = self.fstep if self.Pend is None else 1 / self.Pend
            if self.fend is None:
                self.fend = self.fnyq * self.hifac if self.Pbeg is None else 1 / self.Pbeg

            if self.fend <= self.fbeg:
                raise(ValueError("fend is smaller than (or equal to) fbeg but it must be larger." + \
                               "Choose fbeg and fend so that fend > fbeg."))

            self.freq = arange(self.fbeg, self.fend, self.fstep)
        elif self.fast:
            raise(ValueError("freq and fast cannot be used together."))

        self.nf = len(self.freq)

        # An ad-hoc estimate of the number of independent frequencies (Eq. (24) in ZK_09).
        self.M = (self.fend-self.fbeg) * self.tbase

    def lnL(self, theta, X, Y, e_Y, func):
        # the log-likelihood
       global L, chisqr, wtrms # a blob
       L = []
       Ymod = func(X, theta)
       chisqr = []
       wtrms = [] #weighted rms
       N = len(X)
       for y,e_y,ymod,lnf in zip(Y,e_Y,Ymod,theta[-N:]):
          # loop over data sets
          sigma2 = e_y**2 + lnf**2
          #weight = 1/sigma2
          chisqr += [np.sum((y-ymod)**2/sigma2)]
          L += [-0.5 * np.sum((y-ymod)**2/sigma2 + np.log(2*pi*sigma2))]
          wtrms += [np.sqrt(np.sum((y-ymod)**2/sigma2)/np.sum(1/sigma2))]
       return sum(L)

    def _calcPeriodogram(self):

        self._a, self._b, self.p, self.lnML = np.zeros((4, self.nf))
        self.par = []
        self._off = np.zeros((self.nf, self.Nj))
        self._lnMLj = np.zeros((self.nf, self.Nj))
        self._chij = np.zeros((self.nf, self.Nj))
        self._chisqr = np.zeros((self.nf, self.Nj))
        self._chi = np.zeros(self.nf)
        self._wrmsj = np.zeros((self.nf, self.Nj))
        self._wrms = np.zeros(self.nf)
        self._wtrms = np.zeros((self.nf, self.Nj))

        global L, chisqr, wtrms 
        nll = lambda *args: -self.lnL(*args)
        #a0 = [0.]*self.Nj + [3.]*self.Nj # start guess for first frequency
        a0 = map(np.mean, self.y) + map(np.std, self.y) 
        print(a0)

        # The model with only offset c
        a0 = op.fmin_powell(nll, a0, args=(self.th, self.y, self.e_y, mod_c), disp=False)
        self.lnML0 = sum(L)
        self.lnML0j = L
        self.a0 = a0

        self.chisqr = chisqr
        self.wtrms = wtrms
        print(sum(L), L, a0)
        a = [0., 0.] + list(a0) # start guess for first frequency
        a = [np.median(a0[len(a0)//2:])]*2 + list(a0) # start guess for first frequency
        for k, omega in enumerate(2.*pi*self.freq):
          # Circular frequencies
          X = [omega*th for th in self.th]
          cosX = map(cos, X)
          sinX = map(sin, X)
          # self.lnL(theta, zip(cosX,sinX), self.y, self.e_y, mod_abc)
          #a = op.fmin(nll, a, args=(zip(cosX,sinX), self.y, self.e_y, mod_abc), disp=False)
          a = op.fmin_powell(nll, a, args=(zip(cosX,sinX), self.y, self.e_y, mod_abc), disp=False)
          self.par.append(a)
          self._a[k]= a[0]
          self._b[k]= a[1]
          self._off[k]= a[2:2+self.Nj]
          self._lnMLj[k] = L
          self._chisqr[k] = chisqr
          self.lnML[k] = sum(L)
          self._wtrms[k] = wtrms
          print(k, self.nf, L)
          #pause()
        self.p = self.lnML
        #pause()

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

    def pnorm(self, norm="dlnL"):
        """
        Assign or modify normalization (can be done afterwards).
        Parameters
        ----------
        norm : string, optional
            The normalization to be used (default is 'lnL').
        Examples
        --------
        >>> gls.pnorm('wrms')
        """
        self._normcheck(norm)
        self.norm = norm
        p = self.p
        power = p   # default lnL
        self.label["ylabel"] = norm

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
        elif norm == "dlnL":
            self.powerj = self._lnMLj - self.lnML0j
            power = self.lnML - self.lnML0
            self.label["ylabel"] = "$\Delta$lnL"
        else:
            self.powerj = self._lnMLj.T
            power = self.lnML

        self.power = power

    def _peakPeriodogram(self):
        """
        Analyze the highest periodogram peak.
        """
        # Index with maximum power
        k = argmax(self.p)
        # Maximum power
        self.pmax = pmax = self.p[k]
        self.rms = self._wrmsj[k]
        rms = self.rms[0]
        # Statistics of highest peak
        self.hpstat = p = {}

        # Best parameters
        p["fbest"] = fbest = self.freq[k]
        p["amp"] = amp = sqrt(self._a[k]**2 + self._b[k]**2)
        p["ph"] = ph = arctan2(self._a[k], self._b[k]) / (2.*pi)
        p["T0"]  = self.tmin - ph/fbest
        p["offset"] = self._off[k]
        p["jitter"] = self.par[k][-self.Nj:]

        # Error estimates
        p["amp_err"] = sqrt(2./self.N) * rms
        p["ph_err"] = ph_err = sqrt(2./self.N) * rms/amp/(2.*pi)
        p["T0_err"] = ph_err / fbest
        p["offset_err"] = sqrt(1./self.N) * rms

        # Get the curvature in the power peak by fitting a parabola y=aa*x^2
        if 1 < k < self.nf-2:
            # Shift the parabola origin to power peak
            xh = (self.freq[k-1:k+2] - self.freq[k])**2
            yh = self.p[k-1:k+2] - pmax
            # Calculate the curvature (final equation from least square)
            aa = dot(yh, xh) / dot(xh, xh)
            p["f_err"] = e_f = sqrt(-2./self.N / aa * (1.-self.pmax))
            p["Psin_err"] = e_f / fbest**2
        else:
            self.hpstat["f_err"] = np.nan
            self.hpstat["Psin_err"] = np.nan
            print("WARNING: Highest peak is at the edge of the frequency range.\nNo output of frequency error.\nIncrease frequency range to sample the peak maximum.")

    def sinmod(self, t):
        """
        Calcuate best-fit sine curve.
        Parameters
        ----------
        t : array
            Time array at which to calculate the sine.
        Returns
        -------
        Sine curve : array
            The best-fit sine curve (i.e., that for which the
            power is maximal).
        """
        try:
            p = self.hpstat
            if isinstance(t, (list, tuple)):
               return [p["amp"] * sin(2*np.pi*p["fbest"]*(tj-p["T0"])) + offj for tj,offj in zip(t,p["offset"])]
            else:
               return p["amp"] * sin(2*np.pi*p["fbest"]*(t-p["T0"]))
        except Exception as e:
            print("Failed to calcuate best-fit sine curve.")
            raise(e)

    def info(self):
        """
        Prints some basic statistical output screen.
        """
        print("Generalized LS - statistical output")
        print("-----------------------------------")
        print("Number of input points:     %6d" % self.N)
        print("Weighted mean of dataset:   %f"  % self.rms[0])
        print("Weighted rms of dataset:    %f"  % self.rms[0])
        print("Time base:                  %f"  % self.tbase)
        print("Number of frequency points: %6d" % self.nf)
        print()
        k = argmax(self.p)
        Yfit = self.sinmod([t for t,y,e in self.data])
        W = [1/(e_y**2+jit**2) for e_y,jit in zip(self.e_y, self.hpstat["jitter"])]
        wrmsj = [np.sqrt(np.sum(w*(yfit-y)**2)/np.sum(w)) for y,w,yfit in zip(self.y, W, Yfit)]

        header = "%s:   %10s %10s  %10s %10s" % ("j", "lnL0", "dlnL", "wrms","jit")
        fmt = "%d:   %10.3f %10.3f  %10.3f %10.3f"
        if self.e_y is not None:
           header += " %10s" % "internal error"
           fmt += " %10.3f"
        print(header)

        col = list(self.lnML0j), list(self._lnMLj[k]-self.lnML0j)
        col += wrmsj, self.hpstat["jitter"]
        if self.e_y is not None:
           col += [np.mean(1./e_y**2)**-0.5 for e_y in self.e_y],

        for j,line in enumerate(zip(*col)):
           print(fmt % ((j,)+line))

        wrmsall = np.sqrt(np.sum([np.sum(w*(yfit-y)**2) for y,w,yfit in zip(self.y, W, Yfit)])/sum(map(sum,W)))
        print("-"*60)
        print("all: %10.3f %10.3f  %10.3f\n" % (self.lnML0, self.lnML.max()-self.lnML0, wrmsall))

        print("Best sine frequency:  %f +/- %f" % (self.hpstat["fbest"], self.hpstat["f_err"]))
        print("Best sine period:     %f +/- %f" % (1./self.hpstat["fbest"], self.hpstat["Psin_err"]))
        print("Amplitude:            %f +/- %f" % (self.hpstat["amp"], self.hpstat["amp_err"]))
        print("Phase (ph):           %f +/- %f" % (self.hpstat["ph"], self.hpstat["ph_err"]))
        print("Phase (T0):           %f +/- %f" % (self.hpstat["T0"], self.hpstat["T0_err"]))
        for j,off in enumerate(self.hpstat["offset"]):
           print("Offset %d:             %f +/- %f" % (j, off, self.hpstat["offset_err"]))
        #for j,jit in enumerate(self.hpstat["jitter"]:
        #   print("Jitter %d:             %f +/- %f" % (j, jit, self.hpstat["offset_err"]))
        pause()
        print("-----------------------------------")

    def plot(self, block=False, period=False):
        """
        Create a plot.
        """
        try:
            import matplotlib
            if (matplotlib.get_backend() != "TkAgg"):
                matplotlib.use("TkAgg")
            import matplotlib.pyplot as plt
            from matplotlib.ticker import FormatStrFormatter
        except ImportError:
            raise(ImportError("Could not import matplotlib.pylab."))

        fig = plt.figure()
        fig.subplots_adjust(hspace=0.15, wspace=0.08, right=0.97, top=0.95)
        ax = fig.add_subplot(3, 1, 1)
        ax.set_title("Maximum likelihood periodogram")
        if period:
           ax.set_xscale("log")
           ax.set_xlabel("Period")
        else:
           ax.set_xlabel("Frequency")
        ax.set_ylabel(self.label["ylabel"])
        #pause()
        for pj in self.powerj.T:
           ax.plot(1/self.freq if period else self.freq, pj, '-')
        ax.plot(1/self.freq if period else self.freq, self.power, 'k-')

        fbest, T0 = self.hpstat["fbest"], self.hpstat["T0"]
        # Data and model
        datstyle = {'fmt':'.', 'capsize':0}
        tt = arange(self.tmin, self.tmin+self.tbase, 0.01/fbest)
        ymod = self.sinmod(tt)
        yfit = self.sinmod([t for t,y,e in self.data])
        ax1 = fig.add_subplot(3, 2, 3)
        # ax1.set_xlabel("Time")
        ax1.set_ylabel("Data")
        plt.setp(ax1.get_xticklabels(), visible=False)
        #ax1.errorbar(self.t, self.y, **datstyle)
        for (tj,yj,e_yj),off in zip(self.data, self.hpstat['offset']):
           ax1.errorbar(tj, yj-off, yerr=e_yj, **datstyle)
        ax1.plot(tt, ymod, 'k-')

        tt = arange(T0, T0+1/fbest, 0.01/fbest)
        yy = self.sinmod(tt)
        ax2 = fig.add_subplot(3, 2, 4, sharey=ax1)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax2.get_yticklabels(), visible=False)
        # ax2.set_xlabel("Time")
        # ax2.set_ylabel("Data")
        for (tj,yj,e_yj),off in zip(self.data, self.hpstat['offset']):
           ax2.errorbar(tj*fbest % 1, yj-off, yerr=e_yj, **datstyle)
        xx = tt*fbest % 1
        ii = np.argsort(xx)
        ax2.plot(xx[ii], yy[ii], 'k-')

        # Residuals
        #yres = self.y - yfit
        ax3 = fig.add_subplot(3, 2, 5, sharex=ax1)
        ax3.set_xlabel("Time")
        ax3.set_ylabel("Residuals")
        for (tj,yj,e_yj),yjfit in zip(self.data, yfit):
           ax3.errorbar(tj, yj-yjfit, yerr=e_yj, **datstyle)
        ax3.plot([self.tmin, self.tmin+self.tbase], [0,0], 'k-')

        ax4 = fig.add_subplot(3, 2, 6, sharex=ax2, sharey=ax3)
        # ax4.set_title("Data")
        ax4.set_xlabel("Phase")
        # ax4.set_ylabel("Data")
        plt.setp(ax4.get_yticklabels(), visible=False)
        for (tj,yj,e_yj),yjfit in zip(self.data, yfit):
           ax4.errorbar(tj*fbest % 1,  yj-yjfit, yerr=e_yj, **datstyle)
        ax4.plot([0,1], [0,0], 'k-')

        if hasattr(plt.get_current_fig_manager(), 'toolbar'):
            # check seems not needed when "TkAgg" is set
            plt.get_current_fig_manager().toolbar.pan()
        #t = fig.canvas.toolbar
        #plt.ToggleTool(plt.wx_ids['Pan'], False)
        if block:
           print("Close the plot to continue.")
        else:
           plt.ion()
        plt.show()
        # plt.show(block=block) # unexpected keyword argument 'block' in older matplotlib
        pause()
        return plt

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
        if self.norm == "lnL": return (1.-Pn)**((self.N-3.)/2.)
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
        if self.norm == "lnL": return 1.-Prob**(2./(self.N-3.))
        if self.norm == "Scargle": return -log(Prob)
        if self.norm == "HorneBaliunas": return (self.N-1) / 2. * (1.-Prob**(2./(self.N-3)))
        if self.norm == "Cumming": return (self.N-3) / 2. * (Prob**(-2./(self.N-3.))-1.)
        if self.norm == "wrms": return sqrt(self._YY * Prob**(2./(self.N-3.)))
        if self.norm == "chisq": return self._YY * self.wsum * Prob**(2./(self.N-3.))

    def FAP(self, Pn):
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
            Power threshold.
        Returns
        -------
        FAP : float
            False alarm probability.
        """
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
        Parameters
        ----------
        ofile : string
            Name of the output file.
        """
        with open(ofile, 'w') as f:
            if header:
               f.write("# Generalized Lomb-Scargle periodogram\n")
               f.write("# Parameters:\n")
               if hasattr(self, 'df'):
                  f.write("#    Data file: %s\n" % self.df)
               f.write("#    ofac     : %s\n" % self.ofac)
               f.write("#    norm     : %s\n" % self.norm)
               f.write("# 1) Frequency, 2) Normalized power\n")
            for line in zip(self.freq, self.power, *self.powerj.T.tolist()):
               f.write(("%f "*len(line)) % line + "\n")

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
  argadd('df', nargs='*',
                 help='Data file (three columns: time, data, error). If not specified example will be shown.')
  argadd('-fbeg', '--fbeg', type=float, help="Starting frequency for periodogram.")
  argadd('-fend', '--fend', type=float, help="Stopping frequency for periodogram.")
  argadd('-Pbeg', '--Pbeg', type=float, help="Starting period for periodogram.")
  argadd('-Pend', '--Pend', type=float, help="Stopping period for periodogram.")
  argadd('-ofac', '--ofac', type=float, help="Oversampling factor (default=10).", default=10)
  argadd('-hifac', '--hifac', type=float, help="Maximum frequency (default=1).", default=1)
  argadd('-fast', '--fast', help="Use trigonometric recurrences.", action='store_true')
  argadd('-nojit', '--nojit', help="Optimise jitter.", dest='jit', action='store_false')
  argadd('-norm', '--norm', help="The normalization (default=dlnL).", choices=Gls.norms, default='dlnL')
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

  # A data file has been given.
  try:
   DAT = []
   for dfi in df:
     dat = np.loadtxt(dfi, unpack=True)
     tye = None
     if len(dat) > 1:
        tye = dat[0], dat[1]
     if len(dat) > 2:
        tye += dat[2],
     DAT += [tye]
   tye = DAT
  except Exception as e:
     print("An error occurred while trying to read data file: ")
     print("  " + str(e))
     exit(9)

  gls = Gls(tye, **args)

  if ofile:
     gls.df = df
     gls.toFile(ofile)

  if plot:
     gls.plot(block=True)

#import numpy as np
#x = np.genfromtxt('/home/astro115/carmenes/data/svn/serval/CARM_VIS/J11421+267/J11421+267.rvc.dat', usecols=(0,1,2)).T
#print(x)

##g = Gls(tuple(x),fbeg=0.37,fend=0.384, plot=True)
#ml = Gls(tuple(x),fend=1., plot=True)
##g.plot()
##g=Gls(tuple(x), plot=True)


##import numpy as np; import mlp; g=mlp.Gls(tuple(x),fbeg=0.37,fend=0.384); g.plot()
##import numpy as np; import gls; g=gls.Gls(tuple(x),fbeg=0.37,fend=0.384);
#import numpy as np; import gls; g=gls.Gls(tuple(x),fend=1.);
#from gplot import*
#gplot(ml.freq,ml.p,',',g.freq,g.p)
#g.plot()