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

AbsAp = 'Abs'
#AbsAp = 'Ap'

df_highprob = pd.read_csv( 'highProbTest_MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior_wProbs_wNames' )

#or glob.glob
oc_list = os.listdir("./")
#oc_list = glob.glob("./Clusters")
#print oc_list
print 'Number of clusters: ', len(oc_list), '\n'

for f in oc_list :

	f = 'NGC_6087'
	#f = 'NGC_1901'
	
	# ------------------------------------------------------------------------------------------------------------------------


	dfiso = pd.read_csv( 'bestIsochrone.csv' )	

	G_iso = dfiso['G']
	BP_iso = dfiso['G_BP']
	RP_iso = dfiso['G_RP']

	Z = dfiso['Z'].iloc[0]
	logAge = dfiso['log_age_yr'].iloc[0]
	EBV = dfiso['bestFit_E_B_V'].iloc[0]
	dist = dfiso['bestFit_dist'].iloc[0]

	g_iso  = dfiso['G'] + 0.85926*3.1*EBV - 5. + 5.*np.log10(dist)
	bp_iso = dfiso['G_BP'] + 1.06794*3.1*EBV - 5. + 5.*np.log10(dist)
	rp_iso = dfiso['G_RP'] + 0.65199*3.1*EBV - 5. + 5.*np.log10(dist)

	# ------------------------------------------------------------------------------------------------------------------------

	#df = pd.read_csv( 'Clusters/%s/%s_stars_noRVs.csv'%( f , f ) )	
	df = pd.read_csv( '%s_stars_noRVs_wAbsMag_bestIso.csv'%( f ) )

	# ra, dec, g, bp, rp

	ra = df['ra']
	dec = df['dec']

	G = df['g_abs']
	BP = df['bp_abs']
	RP = df['rp_abs']

	g = df['phot_g_mean_mag'] 
	bp = df['phot_bp_mean_mag'] 
	rp = df['phot_rp_mean_mag']

	# ------------------------------------------------------------------------------------------------------------------------
	# background
	df_bg = pd.read_csv( 'tap10_ngc6087.csv' )

	# ra, dec, g, bp, rp

	ra_bg = df_bg['ra']
	dec_bg = df_bg['dec']

	G_bg = df_bg['phot_g_mean_mag'] - 0.85926*3.1*EBV + 5. - 5.*np.log10(dist)
	BP_bg = df_bg['phot_bp_mean_mag'] - 1.06794*3.1*EBV + 5. - 5.*np.log10(dist)
	RP_bg = df_bg['phot_rp_mean_mag'] - 0.65199*3.1*EBV + 5. - 5.*np.log10(dist)

	g_bg = df_bg['phot_g_mean_mag'] 
	bp_bg = df_bg['phot_bp_mean_mag'] 
	rp_bg = df_bg['phot_rp_mean_mag'] 

	# ------------------------------------------------------------------------------------------------------------------------

	# CLUSTER DATA?

	# 	only works for 1 cepheid
	subdf_highprob_prev = df_highprob.loc[df_highprob['Cluster_name'] == f]

	subdf_highprob = df.loc[df['source_id'] == subdf_highprob_prev['source_id_1'].iloc[0]]

	ra_cep  = subdf_highprob['ra']
	dec_cep = subdf_highprob['dec']

	#G_cep  = subdf_highprob['phot_g_mean_mag'] - 0.85926*3.1*EBV + 5. - 5.*np.log10(dist)
	#BP_cep = subdf_highprob['phot_bp_mean_mag'] - 1.06794*3.1*EBV + 5. - 5.*np.log10(dist)
	#RP_cep = subdf_highprob['phot_rp_mean_mag'] - 0.65199*3.1*EBV + 5. - 5.*np.log10(dist)
	G_cep  = subdf_highprob['g_abs']
	BP_cep = subdf_highprob['bp_abs'] 
	RP_cep = subdf_highprob['rp_abs']

	g_cep = subdf_highprob['phot_g_mean_mag'] 
	bp_cep = subdf_highprob['phot_bp_mean_mag'] 
	rp_cep = subdf_highprob['phot_rp_mean_mag'] 

	# ------------------------------------------------------------------------------------------------------------------------

	r1 = subdf_highprob_prev['r1'].iloc[0]  # degrees
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

	color_OC = 'red'
	color_cep = 'blue'
	color_bg = 'grey'
	color_iso = 'black'

	##f, axarr = plt.subplots(2, sharex=True)
	fig, ax = plt.subplots(2, figsize=(12.5, 6), sharex=True) # largo, ancho


	#fig = plt.figure()
	#ax1 = fig.add_subplot(121)
	#ax2 = fig.add_subplot(122)

	#ax = plt.gca()

	plt.subplot(1, 2, 1)

	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)

	plt.plot(ra_bg, dec_bg, '.', color=color_bg, label='Background', markersize=3, alpha = 0.15)
	plt.plot(ra, dec, 'o', color=color_OC, label='OC stars', markersize=5)
	plt.plot(ra_cep, dec_cep, '*', color=color_cep, label='S Nor', markersize=12)

	plt.grid(color='grey', linestyle='-', linewidth=0.5)

	#plt.axis([min(ra)-0.4, max(ra)+0.2, min(dec)-0.05, max(dec)+0.05])
	plt.axis([max(ra)+0.2, min(ra)-0.2, min(dec)-0.2, max(dec)+0.2])
	

	plt.xlabel('RA [deg]', fontsize=17)
	plt.ylabel('DEC [deg]', fontsize=17)

	#ax.legend(['A simple line'], loc = 1)
	#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	#plt.text(0.95, 0.95, '%s\n  search rad: %.2f deg'%(f.replace('_',' '), search_rad), transform=ax.transAxes, fontsize=10, horizontalalignment='right', verticalalignment='top', bbox=props) # 0,0 is lower-left and 1,1 is upper-right

	#plt.invert_xaxis()
	#plt.xlabel('BP - RP', fontsize=20)
	#plt.ylabel('G', fontsize=20)

	#plt.legend(fontsize=18)

	#plt.tight_layout()
	#if save == True:
	#	plt.savefig('%s_RADEC.eps'%(f))
	#	plt.savefig('%s_RADEC.png'%(f))
	#if show == True:
	#	plt.show()


	plt.subplot(1, 2, 2)
	#fig2, ax2 = plt.subplots(1, figsize=(10, 6), sharex=True) # largo, ancho

	if AbsAp == 'Abs':
		#plt.plot(BP_bg-RP_bg, RP_bg, '.', color=color_bg, label='Background', markersize=3, alpha = 0.15)
		plt.plot(BP-RP, RP, 'o', color=color_OC, label='OC stars', markersize=5)
		plt.plot(BP_iso - RP_iso, RP_iso, '-', color=color_iso, label='Isochrone', linewidth = 4, linestyle = '-.')
		plt.plot(BP_cep-RP_cep, RP_cep, '*', color=color_cep, label='S Nor', markersize=12)
		plt.axis([-1.34, 3.85, 8.5, -6.5])
		plt.xlabel('BP - RP', fontsize=17)
		plt.ylabel('RP', fontsize=17)

	if AbsAp == 'Ap':
		plt.plot(bp_bg-rp_bg, rp_bg, '.', color=color_bg, label='Background', markersize=3, alpha = 0.15)
		plt.plot(bp-rp, rp, 'o', color=color_OC, label='OC stars', markersize=5)
		plt.plot(bp_iso - rp_iso, rp_iso, '-', color=color_iso, label='Isochrone', linewidth = 4, linestyle = '-.')
		plt.plot(bp_cep-rp_cep, rp_cep, '*', color=color_cep, label='S Nor', markersize=12)
		plt.axis([-0.8, 3.4, 18, 4])
		plt.xlabel('bp - rp', fontsize=17)
		plt.ylabel('rp', fontsize=17)

	plt.grid(color='grey', linestyle='-', linewidth=0.5)

	#plt.axis([min(BP-RP)-0.2, max(BP-RP)+0.2, min(RP)-0.2, max(RP)+0.2])
	#plt.axis([-1, 4.5, -5., 8.5])
	
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)

	



	fig.subplots_adjust(wspace=0.18, left=0.10, right=0.97, top=0.91)
	#ax.text(0.95, 0.95, '%s\n Z: %.4f \n logAge: %.2f \n E(B-V): %'%(f.replace('_',' '), search_rad), transform=ax.transAxes, fontsize=10, horizontalalignment='right', verticalalignment='top', bbox=props) # 0,0 is lower-left and 1,1 is upper-right

	#plt.invert_yaxis()
	#if save == True:
	#	plt.savefig('%s_CMD.eps'%(f))
	#	plt.savefig('%s_CMD.png'%(f))
	#if show == True:
	#	plt.show()

	if save == True:
		if AbsAp == 'Abs':
			plt.savefig('%s_radecCMD.pdf'%(f), format='pdf', dpi=1000)
			os.system('pdf2ps %s_radecCMD.pdf %s_radecCMD.ps'%(f, f))
			plt.savefig('%s_radecCMD.png'%(f))
		if AbsAp == 'Ap':
			plt.savefig('%s_radeccmd.pdf'%(f), format='pdf', dpi=1000)
			os.system('pdf2ps %s_radeccmd.pdf %s_radeccmd.ps'%(f, f))
			plt.savefig('%s_radeccmd.png'%(f))	
	if show == True:
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