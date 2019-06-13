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
from astropy.io.votable import parse
from scipy.spatial import distance
from astropy.io import ascii

import timeit

start = timeit.default_timer()

save = True
show = True


df_highprob = pd.read_csv( 'highProbTest_MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior_wProbs_wNames' )

#or glob.glob
oc_list = os.listdir("./Clusters")
#oc_list = glob.glob("./Clusters")
#print oc_list
print 'Number of clusters: ', len(oc_list), '\n'

for f in oc_list :


	f = 'NGC_6087'
	#f = 'NGC_1901'
	
	#df = pd.read_csv( 'Clusters/%s/%s_stars_noRVs.csv'%( f , f ) )	
	df = pd.read_csv( 'Clusters/%s/%s_stars_noRVs_wAbsMag_bestIso.csv'%( f , f ) )

	# ra, dec, g, bp, rp

	ra = df['ra']
	dec = df['dec']

	G = df['g_abs']
	BP = df['bp_abs']
	RP = df['rp_abs']

	dfiso = pd.read_csv( 'Clusters/%s/bestIsochrone.csv'%( f ) )	

	G_iso = dfiso['G']
	BP_iso = dfiso['G_BP']
	RP_iso = dfiso['G_RP']

	Z = dfiso['Z'].iloc[0]
	logAge = dfiso['log_age_yr'].iloc[0]
	EBV = dfiso['bestFit_E_B_V'].iloc[0]
	dist = dfiso['bestFit_dist'].iloc[0]

	# CLUSTER DATA?
	subdf_highprob = df_highprob.loc[df_highprob['Cluster_name'] == f]

	r1 = subdf_highprob['r1'].iloc[0]  # degrees
	if r1 <= 1.:
		search_rad = 1.5*r1
	elif  1. < r1 and r1 <= 2.:
		search_rad = r1     # for the biggest clusters
	elif  2. < r1 and r1 <= 3.:
		search_rad = r1/1.5     # for the biggest clusters	
	else:
		search_rad = r1/2.5            # there are clusters with sizes up to 6.

	print '----------------------------------------------------------------------------'
	print 'Cluster: %s\nr1: %s \tsearch_rad: %s\n\t\t\t Best-Fit Parameters: \t Z = %s \t logAge = %s \t E(B-V) = %s \t dist = %i\n'%(f, r1, search_rad, Z, logAge, EBV, dist)
	print '----------------------------------------------------------------------------'

	#plt.figure()
	#plt.axes()

	#f, axarr = plt.subplots(2, sharex=True)
	fig, ax = plt.subplots(1, figsize=(10, 6), sharex=True) # largo, ancho

	ax = plt.gca()

	plt.plot(ra, dec, 'o', color='red', label='OC stars', markersize=10)
	#plt.plot(per2, amp2, 'h', color='blue', label='RRc', markersize=10)

	plt.axis([min(ra)-0.4, max(ra)+0.2, min(dec)-0.05, max(dec)+0.05])

	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)

	plt.xlabel('RA [deg]', fontsize=20)
	plt.ylabel('DEC [deg]', fontsize=20)

	#ax.legend(['A simple line'], loc = 1)
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	ax.text(0.95, 0.95, '%s\n  search rad: %.2f deg'%(f.replace('_',' '), search_rad), transform=ax.transAxes, fontsize=10, horizontalalignment='right', verticalalignment='top', bbox=props) # 0,0 is lower-left and 1,1 is upper-right

	ax.invert_xaxis()
	#plt.xlabel('BP - RP', fontsize=20)
	#plt.ylabel('G', fontsize=20)

	#plt.legend(fontsize=18)

	plt.tight_layout()
	if save == True:
		plt.savefig('Clusters/%s/%s_RADEC.eps'%(f,f))
		plt.savefig('Clusters/%s/%s_RADEC.png'%(f,f))
	if show == True:
		plt.show()



	fig2, ax2 = plt.subplots(1, figsize=(10, 6), sharex=True) # largo, ancho
	plt.plot(BP-RP, G, 'o', color='black', label='OC stars', markersize=5)
	plt.plot(BP_iso - RP_iso, G_iso, '-', color='red', label='Isochrone', linewidth = 4)

	#ax.text(0.95, 0.95, '%s\n Z: %.4f \n logAge: %.2f \n E(B-V): %'%(f.replace('_',' '), search_rad), transform=ax.transAxes, fontsize=10, horizontalalignment='right', verticalalignment='top', bbox=props) # 0,0 is lower-left and 1,1 is upper-right

	ax2.invert_yaxis()
	plt.show()



	break






# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #