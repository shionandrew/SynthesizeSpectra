import numpy as np
from smooth_gauss import smooth_gauss
import pandas
import numpy as np
import os
import sys
import glob
import math

def smooth_gauss_wrapper(lambda1, spec1, lambda2, dlam_in):
    """
    A wrapper around the Fortran routine smooth_gauss.f, which
    interpolates the synthetic spectrum onto the wavelength array of the
    observed spectrum, while smoothing it to the specified resolution of the
    observed spectrum.
    Adapted into Python from IDL (E. Kirby)

    Parameters
    ----------
    lambda1: array-like: synthetic spectrum wavelength array
    spec1: array-like: synthetic spectrum normalized flux values
    lambda2: array-like: observed wavelength array
    dlam_in: float, or array-like: full-width half max resolution in Angstroms
             to smooth the synthetic spectrum to, or the FWHM as a function of wavelength

    Returns
    -------
    spec2: array-like: smoothed and interpolated synthetic spectrum, matching observations
    """

    if not isinstance(lambda1, np.ndarray): lambda1 = np.array(lambda1)
    if not isinstance(lambda2, np.ndarray): lambda2 = np.array(lambda2)
    if not isinstance(spec1, np.ndarray): spec1 = np.array(spec1)

    #Make sure the synthetic spectrum is within the range specified by the
    #observed wavelength array
    n2 = lambda2.size; n1 = lambda1.size

    def findex(u, v):
        """
        Return the index, for each point in the synthetic wavelength array, that corresponds
        to the bin it belongs to in the observed spectrum
        e.g., lambda1[i-1] <= lambda2 < lambda1[i] if lambda1 is monotonically increasing
        The minus one changes it such that lambda[i] <= lambda2 < lambda[i+1] for i = 0,n2-2
        in accordance with IDL
        """
        result = np.digitize(u, v)-1
        w = [long((v[i] - u[result[i]])/(u[result[i]+1] - u[result[i]]) + result[i]) for i in range(n2)]
        return np.array(w)

    f = findex(lambda1, lambda2)

    #Make it such that smooth_gauss.f takes an array corresponding to the resolution
    #each point of the synthetic spectrum will be smoothed to
    if isinstance(dlam_in, list) or isinstance(dlam_in, np.ndarray): dlam = dlam_in
    else: dlam = np.full(n2, dlam_in)
    dlam = np.array(dlam)

    dlambda1 = np.diff(lambda1)
    dlambda1 = dlambda1[dlambda1 > 0.]
    halfwindow = long(np.ceil(1.1*5.*dlam.max()/dlambda1.min()))

    #pure-Python implementation of smooth_gauss.f
    """
    temp = np.zeros(500); gauss = np.zeros(500)
    spec2 = np.zeros(n2)

    for i in range(n2):

        low = f[i] - halfwindow
        if low < 0: low = 0
        high = f[i] + halfwindow
        if (high < 0): high = 0
        if (high > n1 - 1): high = n1 - 1

        if (low < n1) and (low < high):

            temp2 = 0.
            temp3 = 0.
            temp4 = 0.

            for j in range(low,high+1):

                temp5 = lambda1[j] - lambda2[i]

                if (np.abs(temp5) < dlam[i]*40.):
                    gauss = np.exp(-1.*temp5**2./(2.*dlam[i]**2.))
                    temp2 += gauss
                    temp3 += spec1[j]*gauss
                    temp4 += gauss**2.

            if temp2 > 0.:
                spec2[i] = temp3 / temp2
    """

    #Python wrapped fortran implementation of smooth gauss
    spec2 = smooth_gauss(lambda1, spec1, lambda2, dlam, f, halfwindow)

    return spec2


def main():
    data = pandas.read_csv('Test.csv', delimiter=',')
    flux = []
    wvl = []
    for index, row in data.iterrows():
        wavelength = float(row['wavelength'])
        wvl.append(wavelength)
        flux_ = float(row['flux'])
        flux.append(flux_)

    dlam = 0.5
    smooth_gauss_wrapper.smooth_gauss_wrapper(wvl, flux, wvl2, dlam)

if __name__== "__main__":
    main()
