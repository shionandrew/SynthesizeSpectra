# run_moog.py
# Outputs synthetic spectrum
#
# - createPar: for a given *.atm file and linelist, output *.par file
# - runMoog: runs MOOG for linelist (calls createPar)
#
# Created 4 Jan 18
###################################################################

import os
import sys
import glob
import numpy as np
import math
from interp_atmosphere import checkFile, getAtm, writeAtm
import subprocess
import pandas
import tempfile
import shutil

def createPar(name, atmfile='', linelist='', directory=''):
	"""Create *.par file using *.atm file and linelist."""

	# Open linelist and get wavelength range to synthesize spectrum
	wavelengths = np.genfromtxt(linelist, skip_header=1, usecols=0)
	wavelengthrange = [ math.floor(wavelengths[0]),math.ceil(wavelengths[-1]) ]

	# Define filename
	filestr = directory + name + '.par'

	# Check if file already exists
	exists, readytowrite = checkFile(filestr)
	if readytowrite:

		# Outfile names:
		out1 = '\''+directory+'/'+name+'.out1\''
		out2 = '\''+directory+'/'+name+'.out2\''

		# If file exists, open file
		with open(filestr, 'w+') as file:

			# Print lines of .par file
			file.write('synth'+'\n')
			file.write('terminal       '+'\'x11\''+'\n')
			file.write('standard_out   '+out1+'\n')
			file.write('summary_out    '+out2+'\n')
			file.write('model_in       '+'\''+atmfile+'\''+'\n')
			file.write('lines_in       '+'\''+linelist+'\''+'\n')
			file.write('strong        0'+'\n')
			file.write('atmosphere    1'+'\n')
			file.write('molecules     1'+'\n')
			file.write('damping       1'+'\n')
			file.write('trudamp       0'+'\n')
			file.write('lines         1'+'\n')
			file.write('flux/int      0'+'\n')
			file.write('plot          0'+'\n')
			file.write('synlimits'+'\n')
			file.write('  '+'{0:.3f}'.format(wavelengthrange[0])+' '+'{0:.3f}'.format(wavelengthrange[1])+'  0.02  1.00'+'\n')
			file.write('obspectrum    0')

	return filestr, wavelengthrange

def runMoog(temp, logg, fe, alpha, directory='/home/seandrew/', elements=None, abunds=None, solar=None, lines=''):
	"""Run MOOG for each Mn linelist and splice spectra.
	Inputs:
	temp 	 -- effective temperature (K)
	logg 	 -- surface gravity
	fe 		 -- [Fe/H]
	alpha 	 -- [alpha/Fe]
	Keywords:
	dir 	 -- directory to write MOOG output to [default = '/home/seandrew/moogspectra/']
	elements -- list of atomic numbers of elements to add to the *.atm file
	abunds 	 -- list of elemental abundances corresponding to list of elements
	lines    -- list of wavelengths for nearly all elements
	Outputs:
	spectrum -- spliced synthetic spectrum
	"""

	spectrum  = []

	# Define temporary directory to store tempfiles
	tempdir = tempfile.mkdtemp()

	# Create identifying filename (including all parameters + linelist used)
	name = getAtm(temp, logg, fe, alpha, directory='') # Add all parameters to name

	# Create *.atm file (for use with each linelist)
	atmfile = writeAtm(temp, logg, fe, alpha, elements=elements, abunds=abunds, solar=solar, dir= tempdir)

	# Create *.par file
	parname = name
	parfile, wavelengthrange = createPar(parname, atmfile, '/raid/caltech/rb_zr/alllines', directory= tempdir)
	# Run MOOG
	p = subprocess.Popen(['MOOG', parfile], cwd='/home/seandrew/moog17scat/', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# Wait for MOOG to finish running
	p.communicate()


	# Create arrays of wavelengths and fluxes
	outfile = tempdir+'/'+parname+'.out2'

	wavelength = np.linspace(wavelengthrange[0],wavelengthrange[1],math.ceil((wavelengthrange[1]-wavelengthrange[0])/0.02), endpoint=True)
	data = pandas.read_csv(outfile, skiprows=[0,1,-1], delimiter=' ').to_numpy()
	flux = data[~np.isnan(data)][:-1]
	spectrum.append([1.-flux, wavelength])

	# Output synthetic spectrum in a format that continuum_div functions will understand (list of arrays)

	# Clean out the temporary directory
	shutil.rmtree(tempdir)

	return spectrum

