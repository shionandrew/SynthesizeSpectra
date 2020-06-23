from os.path import isfile, join
from functools import reduce
import scipy.ndimage
import numpy as np
import itertools
import gzip
import glob
import sys
import os

def enforce_grid_check(teff, logg, feh, alphafe):

    """
    Enforce the grid constraints from Kirby et al. 2009 on a given combination
    of atmospheric model parameters.

    Parameters
    ----------
    teff: float: effective temperature
    logg: float: surface gravity
    feh: float: [Fe/H]
    alphafe: float: [alpha/Fe]

    Returns
    -------
    in_grid: boolean: if True, the specified parameters are within the K08 grid range
    """

    #Check that the effective temperature is within limits
    teff_lims = [3500., 8000.]
    if (teff < teff_lims[0]) or (teff > teff_lims[1]):
        in_grid = False; key = 0
        return in_grid, key

    logg_hi = 5.0
    if teff < 7000.: logg_lo = 0.0
    else: logg_lo = 0.5
    logg_lims = [logg_lo, logg_hi]

    #Check if the specified surface gravity is within limits
    if (logg < logg_lims[0]) or (logg > logg_lims[1]):
        in_grid = False; key = 1
        return in_grid, key

    #Check that the specified metallicity is within limits
    feh_lims = [-5., 0.]

    #Put checks in place based on the limits of the grid imposed by
    #difficulty in model atmosphere convergence
    teff_vals = np.array([3600., 3700., 3800., 3900., 4000., 4100.])
    feh_thresh = np.array([-4.9, -4.8, -4.8, -4.7, -4.4, -4.6])
    logg_thresh = np.array([[1.5], [2.5, 3.5, 4.0, 4.5], [4.0, 4.5],
                           [2.5, 3.0, 3.5, 4.0, 4.5, 5.0], [4.5, 5.0], [4.5, 5.0]])

    if teff in teff_vals:
        where_teff = np.where(teff_vals == teff)[0][0]
        if logg in logg_thresh[where_teff]:
            if (feh < feh_thresh[where_teff]) or (feh > feh_lims[1]):
                in_grid = False; key = 2
                return in_grid, key
    else:
        if (feh < feh_lims[0]) or (feh > feh_lims[1]):
            in_grid = False; key = 2
            return in_grid, key

    #Check that the alpha enhancement is within limits
    alpha_lims = [-0.8, 1.2]
    if (alphafe < alpha_lims[0]) or (alphafe > alpha_lims[1]):
        in_grid = False; key = 3
        return in_grid, key

    in_grid = True; key = 4
    return in_grid, key

def read_synth(Teff = np.nan, Logg = np.nan, Feh = np.nan, Alphafe = np.nan, fullres=True,\
    filename= '', gauss_sigma = 1, lambda_sigma = -1, data_path='', verbose=False,\
    start=4100., sstop=6300., file_ext = '.bin', title='', numiso=2):

    """
    Read the ".bin.gz" file containing the synthetic spectrum information
    that it output by running MOOGIE

    Parameters:
    ------------
    Teff: float: effective temperature (K) of synthetic spectrum
    Logg: float: surface gravity (log cm s^(-2)) of synthetic spectrum
    Feh: float: iron abundance ([Fe/H]) (dex) of synthetic spectrum
    Alphafe: float : alpha-to-iron ratio [alpah/Fe] (dex) of synthetic spectrum
    fullres: boolean: if True, then use the unbinned, full-resolution version of the
                      synthetic spectrum
    filename: string: the filename of the desired synthetic spectrum. If the parameters
                       describing the synthetic spectrum are given, then a filename is
                       not necessary
    data_path: string: the path leading to the parent directory containing the synthetic
                       spectrum data
    verbose: boolean: if True, then print out statements during run time
    start: float: start wavelength of synthetic spectrum
    sstop: float: stop wavelength of synthetic spectrum
    file_ext: string: file extension of the filename to be read, default '.bin'

    Returns:
    ---------
    wvl: array: the wavelength range covered by the synthetic spectrum, depending on
                whether it is full resolution or binned, and on stop and start wavelengths
    relflux: array: the normalized flux of the synthetic spectrum

    """

    if file_ext != '.bin': linestart = 3 + numiso
    else: linestart = 3

    #Determine which directory to point to (binned/full resolution)
    if file_ext == '.bin':
        if fullres == True:
            directory = 'synths/'
            step = 0.02
        else:
            directory = 'bin/'
            step = 0.14
    else:
        directory = ''
        step = 0.01

    #If given, determine the degree of gaussian smoothing of the spectrum
    if lambda_sigma > 0.:
        gauss_sigma = round(lambda_sigma/step)

    if filename == '':

        #Check if the parameters are specified, if the filename is not
        if np.all(np.isnan([Teff,Logg,Feh,Alphafe])) == True:
            print("Error: must define teff, logg, feh, and alphafe")
            return np.nan, np.nan, None, None

        path = data_path+directory #full path name

        #Redefine the parameters according to the file naming convention
        title, filename = construct_title_filename(Teff, Logg, Feh, Alphafe)

        if file_ext == '.bin': bin_gz_file = filename
        else: out_file = filename

        filename = path+filename

    #Otherwise, if a filename is specified
    else:
        if file_ext == '.bin': bin_gz_file = filename
        else: out_file = filename

    if verbose == True: print(filename)

    if file_ext == '.bin':

        #Open and read the contents of the gzipped binary file without directly
        #unzipping, for enhanced performance
        with gzip.open(filename, 'rb') as f:
            bstring = f.read()
            flux = np.fromstring(bstring, dtype=np.float32)

        #Alternative method of opening and reading in the file that requires
        #explicit unzipping and zipping
        """
        filename_decompress = filename[:-3] #get rid of the filename extension
        success = os.system('gunzip %s'%filename) #attempt to unzip the file

        if success != 0:
            if verbose==True:
                print "Error unzipping %s"%filename
            return np.nan, np.nan, None, None

        flux = np.fromfile(filename_decompress,'f4')
        success = os.system('gzip %s'%filename_decompress) #now zip the file back up

        if success != 0:
            if verbose==True:
                print "Error zipping %s"%filename
            return np.nan, np.nan
        """

    else:
        with open(filename, 'r') as f:
            lines = f.readlines()
            f.close()
        flux = []
        for line in lines[linestart:]:
            i = 0
            while i < len(line)-2:
                el = float(line[i:i+7].strip())
                flux.append(el)
                i += 7

    wvl_range = np.arange(start, sstop+step, step)
    wvl = 0.5*(wvl_range[1:] + wvl_range[:-1])

    if gauss_sigma != 0.0:
        relflux = scipy.ndimage.filters.gaussian_filter(1.0 - flux, gauss_sigma)
    else:
        relflux = 1. - flux

    if file_ext == '.bin':
        return wvl, relflux, title, bin_gz_file[11:-7]
    else:
        return wvl, relflux, title, out_file[:-5]

def construct_title_filename(Teff=np.nan, Logg=np.nan, Feh=np.nan, Alphafe=np.nan,
file_ext='.bin', interp=False, Dlam=np.nan):

    #Redefine the parameters according to the file naming convention
    if not interp:
        teff = round(Teff/100.)*100
        logg = round(Logg*10.)
        feh = round(Feh*10.)
        alphafe = round(Alphafe*10.)
    else:
        teff = np.round(Teff, decimals=0)
        logg = np.round(Logg*10., decimals=2)
        feh = np.round(Feh*10., decimals=2)
        alphafe = np.round(Alphafe*10., decimals=2)
        dlam = np.round(Dlam*10.,decimals=2)

    if logg >= 0.:
        gsign = '_'
    else: gsign = '-'

    if feh >= 0.:
        fsign = '_'
    else: fsign = '-'

    if alphafe >= 0.:
        asign = '_'
    else: asign = '-'

    title = r'T$_{eff}$=%g, log(g)=%g, [Fe/H]=%g, [$\alpha$/Fe]=%g,\
    $\Delta\lambda$=%g'%(np.round(Teff, decimals=0), np.round(Logg, decimals=2),\
    np.round(Feh, decimals=2), np.round(Alphafe, decimals=2), np.round(Dlam, decimals=2))

    if file_ext == '.bin':

        bin_gz_file = "t%2i/g%s%2i/t%2ig%s%2if%s%2ia%s%2i.bin.gz"%(teff,gsign,abs(logg),teff,gsign,abs(logg),fsign,abs(feh),asign,abs(alphafe))
        bin_gz_file = bin_gz_file.replace(" ","0")
        filename = bin_gz_file
    else:
        out_file = "t%2ig%s%2if%s%2ia%s%2i.out2"%(teff,gsign,abs(logg),fsign,abs(feh),asign,abs(alphafe))
        out_file = out_file.replace(" ","0")
        filename = out_file

    return title, filename

def read_interp_synth(teff=np.nan, logg=np.nan, feh=np.nan, alphafe=np.nan,
fullres=False, data_path='', start=4100., sstop=6300., npar=4, gauss_sigma=0., hash=None):

    """
    Construct a synthetic spectrum in between grid points based on linear interpolation
    of synthetic spectra in the MOOGIE grid

    Parameters:
    -----------
    Teff: float: effective temperature (K) of synthetic spectrum
    Logg: float: surface gravity (log cm s^(-2)) of synthetic spectrum
    Feh: float: iron abundance ([Fe/H]) (dex) of synthetic spectrum
    Alphafe: float : alpha-to-iron ratio [alpah/Fe] (dex) of synthetic spectrum
    fullres: boolean: if True, then use the unbinned, full-resolution version of the
                      synthetic spectrum
    filename: string: the filename of the desired synthetic spectrum. If the parameters
                       describing the synthetic spectrum are given, then a filename is
                       not necessary
    data_path: string: the path leading to the parent directory containing the synthetic
                       spectrum data
    verbose: boolean: if True, then print out statements during run time
    start: float: start wavelength of synthetic spectrum
    sstop: float: stop wavelength of synthetic spectrum
    file_ext: string: file extension of the filename to be read, default '.bin'
    npar: integer: number of parameters used to describe a synthetic spectrum
    hash: dict, optional: a dictionary to use to store memory concerning which synthetic
          spectra have been read in. Should be initliazed externally as an empty dict.

    Returns:
    --------
    wvl: array: the wavelength range covered by the synthetic spectrum, depending on
                whether it is full resolution or binned, and on stop and start wavelengths
    relflux: array: the normalized flux of the synthetic spectrum
    """

    #Define the points of the 4D grid
    teff_arr = np.arange(3500., 5600., 100.).tolist() + np.arange(5600., 8200., 200.).tolist()
    teff_arr = np.round(np.array(teff_arr), decimals=0)

    logg_arr = np.round(np.arange(0., 5.5, 0.5), decimals=1)
    feh_arr = np.round(np.arange(-5., 0.1, 0.1), decimals=2)
    alphafe_arr = np.round(np.arange(-0.8, 1.3, 0.1), decimals=2)
    alphafe_arr[8] = 0.

    #First check that given synthetic spectrum parameters are in range
    #If not in range, throw an error
    #Consider implementing a random selection of nearby grid points if the parameters
    #do go out of range, for values near the edge of the grid

    #in_grid, key = enforce_grid_check(teff, logg, feh, alphafe)
    in_grid,_ = enforce_grid_check(teff, logg, feh, alphafe)
    if not in_grid: return

    """
    if not in_grid:

        sys.stderr.write('\nInput parameters {} = {} out of grid range\n'.format(['Teff',\
        'logg', '[Fe/H]', '[alpha/Fe]'], [teff, logg, feh, alphafe]))

        while not in_grid:

            if key == 0: teff = np.random.uniform(teff_arr[0], teff_arr[-1])
            if key == 1: logg = np.random.uniform(logg_arr[0], logg_arr[-1])
            if key == 2: feh = np.random.uniform(feh_arr[0], feh_arr[-1])
            if key == 3: alphafe = np.random.uniform(alphafe_arr[0], alphafe_arr[-1])

            sys.stderr.write('Selecting new parameters {} = {}\n'.format(['Teff', 'Logg',\
            '[Fe/H]', '[alpha/Fe]'], [teff, logg, feh, alphafe]))

            in_grid, key = enforce_grid_check(teff, logg, feh, alphafe)
    """

    params = np.array([teff, logg, feh, alphafe])
    params_grid = np.array([teff_arr, logg_arr, feh_arr, alphafe_arr])

    #Now identify the nearest grid points to the specified parameter values
    ds = []; nspecs = []; iparams = []
    for i in range(npar):

        #The specified parameter value is a grid point
        w = np.digitize(params[i], params_grid[i])
        if params[i] in params_grid[i]:
            iparam = np.array([w-1, w-1])
            d = [1.]
            nspec = 1
            ds.append(d)

        #The specified parameter value is in between grid points
        else:
            if w == (len(params_grid[i])): w -= 1
            iparam = np.array([w-1, w])
            d = params_grid[i][iparam] - params[i]
            d_rev = np.abs(d[::-1])
            nspec = 2
            ds.append(d_rev)

        nspecs.append(nspec)
        iparams.append(iparam)

    #Now, based on the nearest grid points, construct the linearly interpolated
    #synthetic spectrum

    #Determine the number of pixels in a synthetic spectrum based on whether
    #the spectrum is binned or unbinned
    if fullres: step = 0.02
    else: step = 0.14

    #Calculate the number of pixels in the synthetic spectrum, and initialize the
    #interpolated synthetic spectrum array
    npixels = len(np.arange(start, sstop, step))
    synth_interp = np.zeros(npixels)

    #Function for loading a specified synthetic spectrum
    def load_synth(p):

        teffi, loggi, fehi, alphafei = p

        if hash is not None:

            #First construct the filename corresponding to the parameters, to use for
            #testing whether we should read in the specified synthetic spectrum

            _,filename = construct_title_filename(Teff=teff_arr[teffi],
            Logg=logg_arr[loggi], Feh=feh_arr[fehi], Alphafe=alphafe_arr[alphafei])
            key = filename[11:-7]

            if key not in hash.keys():

                _,synthi,_,_ = read_synth(Teff=teff_arr[teffi], Logg=logg_arr[loggi],
                Feh=feh_arr[fehi], Alphafe=alphafe_arr[alphafei], fullres=fullres,
                verbose=False, gauss_sigma=gauss_sigma, data_path=data_path)

                hash[key] = synthi

            #If the key is already present in the hash table, then find it and load the data
            else: synthi = hash[key]

        else:
            _,synthi,_,_ = read_synth(Teff=teff_arr[teffi], Logg=logg_arr[loggi],
            Feh=feh_arr[fehi], Alphafe=alphafe_arr[alphafei], fullres=fullres,
            verbose=False, gauss_sigma=gauss_sigma, data_path=data_path)

        return synthi

    """
    #Function for calculating the product of the prefactors
    def dprod(q):
        ds0, ds1, ds2, ds3 = q
        return ds0*ds1*ds2*ds3

    #Load each nearby synthetic spectrum on the grid to linearly interpolate
    #to calculate the interpolated synthetic spectrum

    params_tup = list(itertools.product(*iparams))
    ds_tup = list(itertools.product(*ds))

    synthis = map(load_synth, params_tup)
    dprods = map(dprod, ds_tup)

    for i in range(len(dprods)):
        for m in range(npixels):
            synth_interp[m] += dprods[i]*synthis[i][m]
    """

    #Load each nearby synthetic spectrum on the grid to linearly interpolate
    #to calculate the interpolated synthetic spectrum
    for i in range(nspecs[0]):
        for j in range(nspecs[1]):
            for k in range(nspecs[2]):
                for l in range(nspecs[3]):

                    p = [iparams[0][i], iparams[1][j], iparams[2][k], iparams[3][l]]
                    synthi = load_synth(p)

                    """
                    _,synthi,_,_ = read_synth(Teff=teff_arr[iparams[0][i]],
                    Logg=logg_arr[iparams[1][j]], Feh = feh_arr[iparams[2][k]],
                    Alphafe = alphafe_arr[iparams[3][l]], fullres=fullres, verbose=False,
                    gauss_sigma=0., data_path=data_path)
                    """

                    for m in range(npixels):
                        synth_interp[m] += ds[0][i]*ds[1][j]*ds[2][k]*ds[3][l]*synthi[m]

    facts = []
    for i in range(npar):
        if nspecs[i] > 1: fact = params_grid[i][iparams[i][1]] - params_grid[i][iparams[i][0]]
        else: fact = 1
        facts.append(fact)

    synth_interp /= reduce(lambda x, y: x*y, facts)
    wave = np.arange(start, sstop, step)

    return wave, synth_interp


def main():
    wvl_b, flux_b = read_interp_synth(teff=5106, logg=2.312, feh=-1.72, alphafe=0.12, data_path='/raid/gridie/')
    wvl_r, flux_r = read_interp_synth(teff=5106, logg=2.312, feh=-1.72, alphafe=0.12, data_path='/raid/grid7/',start=6300., sstop=9100.)
    print(wvl_r)
    wvl = np.concatenate((wvl_b,wvl_r))
    flux = np.concatenate((flux_b,flux_r))
    outputfile =  '/home/seandrew/synthesize_spectra/510_2.312_-1.72_0.12.csv' 
    with open(outputfile, 'w+') as file:
        for datapoint in range(len(wvl)):
            wavelength_ = str(wvl[datapoint])
            flux_ = str(flux[datapoint])
            file.write(wavelength_ + ',' + flux_ + '\n')

if __name__== "__main__":
    main()
