#!/usr/bin/python2.7
 
import numpy as np#para cargar numpy
from pylab import *#cargar mtplotlib
import os#cargar herramientas de archivos
import re # use regular patterns
import sys, getopt # system commands
import string # string functions4
import math
#from scipy import linalg as scipylinalg
#from scipy import stats
#from scipy import ndimage
#from scipy import signal as scipysignal
import scipy
from pylab import *
from numpy import nan
import random
import subprocess
import pandas as pd
import glob
from numpy.linalg import inv
from scipy.stats import chi2
from scipy import stats

import timeit

start = timeit.default_timer()

# open catalog
catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_RVMelnikOnly.csv'
#catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_RVMelnikOnly_SELECTED.csv'


df = pd.read_csv(catalog_name)

largo = len(df)


i = 0
while i < largo:

	if (i+1)%200==0:
		print '\n\t', i+1, '/', largo	# mean values

	mean_pmraval = df['pmra'].iloc[i]
	mean_pmdecval = df['pmdec'].iloc[i]
	mean_plxval = df['parallax'].iloc[i]
	mean_rvval = df['HRV_melnik'].iloc[i]

	# sd's
	sd_pmraval = df['pmra_error'].iloc[i]
	sd_pmdecval = df['pmdec_error'].iloc[i]
	sd_plxval = df['parallax_error'].iloc[i]
	sd_rvval = df['e_HRV_melnik'].iloc[i]

	## gaussian distribution
	# value aroung the mean value, in standard deviations
	N = 1

	a = np.random.normal(mean_pmraval, N*sd_pmraval)
	b = np.random.normal(mean_pmdecval, N*sd_pmdecval)
	c = np.random.normal(mean_plxval, N*sd_plxval)
	d = np.random.normal(mean_rvval, N*sd_rvval)

	#print mean_rvval, sd_rvval
	#print d
	#print abs(mean_rvval - d)
	#print '--------------------------------------------------------\n'

	# offsets
	off_pmra  = 0.1 # mas/yr
	off_pmdec = 0.1 # mas/yr
	off_plx = 0.1 # mas
	off_rv = 5 # km sec


	# add values
	#if df['pmra'].iloc[i] != 0 and df['pmra_error'].iloc[i]!=0:
	#	df['pmra_error'].iloc[i]     = abs(mean_pmraval - a) + off_pmra
	#if df['pmdec'].iloc[i] != 0 and df['pmdec_error'].iloc[i]!=0:
	#	df['pmdec_error'].iloc[i]    = abs(mean_pmdecval - b) + off_pmdec
	#if df['parallax'].iloc[i] != 0 and df['parallax_error'].iloc[i]!=0:
	#	df['parallax_error'].iloc[i] = abs(mean_plxval - c) + off_plx
	#if df['HRV_melnik'].iloc[i] != 0 and df['e_HRV_melnik'].iloc[i]!=0:
	#	df['e_HRV_melnik'].iloc[i]   = abs(mean_rvval - d) + off_rv

	i = i+1

# save catalog with artificial errors

df.to_csv( './'+catalog_name.replace('.csv','_ErrsTweak.csv'), index=False )

# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #