
# synthesize_spectra.py
# Runs MOOG to generate spectra with specified input parameters
# Outputs synthetic spectrum
#
# - synthesize_spectra: runs MOOG for each set of specified parameters (calls createPar), outputs spectrum
#
# Created 9 Jan 2020
###################################################################

import os
import sys
import glob
import numpy as np
import math
from run_moog import createPar, runMoog
import pandas
import subprocess
from read_synth import read_interp_synth
from smooth_gauss_wrapper import smooth_gauss_wrapper

def synthesize_spectra(galaxy):
	"""synthesizes spectra for stars in a given galaxy
	Inputs:
    galaxy-- specifies one of the six dwarf galaxies (filename + directory where spectra will be saved)
    """
	inputfilename = '/home/seandrew/synthesize_spectra/' + galaxy + '/' + galaxy + 'Params.csv'
	#inputfile contains list of temperatures, logg, fe, and alpha
	data = pandas.read_csv(inputfilename, delimiter=',')

	outputdirectory = '/home/seandrew/synthesize_spectra/' + galaxy + '/'
	#directory for output spectra
	for index, row in data.iterrows():
		logtemp = float(row['Logtemp'])
		temp = math.pow(10,logtemp)
		print(temp)
		logg = row['logg']
		fe = row['[Fe/H]']
		alpha = row['[a/Fe]']
		outputfile = outputdirectory + str(temp) + '_' + str(logg) + '_' + str(fe) + '_' + str(alpha) + '.csv'
		print(outputfile)

		wvl_b, flux_b = read_interp_synth(teff=temp, logg=logg, feh=fe, alphafe=alpha, data_path='/raid/gridie/')
		wvl_r, flux_r = read_interp_synth(teff=temp, logg=logg, feh=fe, alphafe=alpha, data_path='/raid/grid7/',start=6300., sstop=9100.)
		wvl = np.concatenate((wvl_b,wvl_r))
		flux = np.concatenate((flux_b,flux_r))
		with open(outputfile, 'w+') as file:
			for datapoint in range(len(wvl)):
				wavelength_ = str(wvl[datapoint])
				flux_ = str(flux[datapoint])
				file.write(wavelength_ + ',' + flux_ + '\n')
		'''spectrum = runMoog(temp, logg, fe, alpha)
		with open(outputfile, 'w+') as file:
			file.write('wavelength, flux \n')
			for datapoint in range(len(spectrum[0][0])):
				wavelength = str(spectrum[0][1][datapoint])
				flux = str(spectrum[0][0][datapoint])
				file.write(wavelength + ',' + flux + '\n')'''

def main():
    synthesize_spectra('Sculptor')
    synthesize_spectra('Sextans')
    synthesize_spectra('Fornax')
    synthesize_spectra('UrsaMinor')
    synthesize_spectra('Draco')
    synthesize_spectra('NGC6822')
    #synthesize_spectra('BootesI')

if __name__== "__main__":
    main()
