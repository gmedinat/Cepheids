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

# IDEA: llenar con nans los campos que tengan X = 0 and e_X = 0.
# luego crear un c distinto dependiendo de que valores esten disponibles.

save = False
show = True
verbose = False

#catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_dupRemoved_wPrior.csv'
catalog_name = 'test_membership.csv'

#cols = ['id','coord_ra','coord_dec', 'object', 'visit', 'ccdnum', 'MJDObs', 'base_SdssShape_psf_xx', 'base_SdssShape_psf_yy', 'AIRMASS', 'SEEING', 'EXPTIME', 'base_PsfFlux_flux', 'base_PsfFlux_fluxSigma']
#df = pd.read_csv(indir+'/RR%s/'%file_n+archivo, usecols=cols)
df = pd.read_csv(catalog_name)

#df = df[['id','coord_ra','coord_dec', 'RA_deg', 'DEC_deg', 'object', 'visit', 'ccdnum', 'MJDObs', 'base_SdssShape_psf_xx', 'base_SdssShape_psf_yy', 'AIRMASS', 'SEEING', 'EXPTIME', 'base_PsfFlux_flux', 'base_PsfFlux_fluxSigma']]#1

# Covariance
def cov(x, y):
    xbar, ybar = x.mean(), y.mean()
    return np.sum((x - xbar)*(y - ybar))/(len(x) - 1)

# number of degrees of freedom


df['e_Dist'] = df['Dist']*0.2 # 20% of the distance assumed
df['parallax_cluster'] = 1000./df['Dist']
df['e_parallax_cluster'] = 1000.*df['e_Dist']/(df['Dist']*df['Dist']) 


# add columns with differences
# delta                 cluster   -   cepheid
df['parallax_diff'] = df['parallax_cluster'] - df['parallax']
df['vr_diff']       = df['RV_2'] - df['radial_velocity']
df['pmRA_diff']     = df['pmRA_2'] - df['pmra']
df['pmDEC_diff']    = df['pmDE_2'] - df['pmdec']
#df['FeH_diff']      = df['[]'] - df['mag_2']
#df['age_diff']      = df['Age_D02'] - df['mag_2'] # or logt_K13


n_nans = 0

largo = len(df)

i = 0
while i < largo:

	#print '\n\t', i+1, '/', largo

	subdf = df.iloc[i]

	x_df = subdf[['parallax_diff','vr_diff','pmRA_diff','pmDEC_diff']]
	
	x = x_df.as_matrix()
	x = x.reshape(len(x), 1)
	#xT = np.transpose(x)
	xT = x.T

	if verbose == True:	

		print xT.shape, xT
		print x.shape, x

	# covariance matrix cluster    parallax   RV   pmRA  pmDEC	
	M1 = np.matrix([[subdf['e_parallax_cluster']**2, 0, 0, 0], 
				[0, subdf['e_RV']**2, 0, 0],
				[0, 0, subdf['e_pmRA']**2, 0],
				[0, 0, 0, subdf['e_pmDE']**2]]) 

	if verbose == True:	
		print 'M1', M1.shape
		print M1, '\n'


	# cepheids:
	# parallax_pmra_corr
	# parallax_pmdec_corr
	# pmra_pmdec_corr

	# rv par ?
	# rv pmra ?
	# rv pmdec ?

	cov_par_pmra   = subdf['parallax_pmra_corr']*subdf['parallax_error']*subdf['pmra_error']
	cov_par_rv     = 0.
	cov_par_pmdec  = subdf['parallax_pmdec_corr']*subdf['parallax_error']*subdf['pmdec_error']
	cov_rv_pmra    = 0.
	cov_rv_pmdec   = 0.
	cov_pmra_pmdec = subdf['pmra_pmdec_corr']*subdf['pmra_error']*subdf['pmdec_error']

	# covariance matrix cepheid
	M2 = np.matrix([[subdf['parallax_error']**2, cov_par_rv, cov_par_pmra, cov_par_pmdec], 
				[cov_par_rv, subdf['radial_velocity_error']**2, cov_rv_pmra, cov_rv_pmdec],
				[cov_par_pmra, cov_rv_pmra, subdf['pmra_error']**2, cov_pmra_pmdec],
				[cov_par_pmdec, cov_rv_pmdec, cov_pmra_pmdec, subdf['pmdec_error']**2]])

	if verbose == True:	
		print 'M2', M2.shape
		print M2, '\n'

	# sum of covariance matrices
	Matrix = M1+M2

	Matrix = np.nan_to_num(Matrix) # ESTA LINEA ES POLEMICA! REMOVER EVENTUALMENTE

	if verbose == True:	
		print 'Matrix', Matrix.shape
		print Matrix, '\n'

	# inverse of the sum of covariance matrices

	MatrixInv = Matrix.getI()
	#MatrixInv = inv(Matrix)

	if verbose == True:	

		print 'InvMatrix', MatrixInv.shape
		print MatrixInv, '\n'

		print Matrix * MatrixInv


	# ------------------------------------------------------------------------------------------- #

	#c = xT.dot(x) # 1 x n  dot  n x 1 = scalar
	ca = xT.dot(MatrixInv.dot(x))  
	c = np.squeeze(np.asarray(ca))
	if verbose == True:	
		print ca, c

	if i == 0:
		c_distribution = c
	else: 
		c_distribution = np.append(c_distribution, c)

	del subdf

	#https://stackoverflow.com/questions/11725115/p-value-from-chi-sq-test-statistic-in-python
	#print c
	#print stats.chi2.pdf(float(c), 4)
	#print 1 - stats.chi2.cdf(float(c), 4)


	i = i+1

	if pd.isnull(c) == True:
		print c
		n_nans = n_nans+1
	#	sys.exit(0)

	#if i == 500:
	#	break

del df

#print c_distribution, c_distribution.shape
c_distribution = c_distribution.reshape(len(c_distribution), 1)
#print c_distribution, c_distribution.shape

c_distribution = c_distribution[~pd.isnull(c_distribution)] # remove nans
c_distribution = c_distribution.reshape(len(c_distribution), 1)
#c_distribution = c_distribution[~np.isnan(c_distribution)] # remove nans
#c_distribution = c_distribution[np.isfinite(c_distribution)] # remove nans

print 'number of nans / total: ', n_nans,'/',largo, '\t\t min, max: ', c_distribution.min(), c_distribution.max() 

fig, ax = plt.subplots()
#plt.hist(c_distribution, bins=30)
plt.hist(c_distribution, bins=30, density=True)
#ax.set_xlim(0.1, 1)
#ax.set_ylim(0, 53)
#ax.set_xticks([16,17,18,19,20,21,22])
#ax.set_yticks([5,10,15,20,25,30,35,40,45,50])


#print chi2.ppf(0.01, df), chi2.ppf(0.99999999, df)
#x1 = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.9999999, df), 100)
#x2 = np.linspace(chi2.ppf(0.01, df, scale=10), chi2.ppf(0.9999999999, df, scale=10), 100)
#print x1, x2
#ax.plot(c_distribution, chi2.pdf(c_distribution, 4), 'b-', lw=5, alpha=0.6, label='x chi2 pdf')

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.xlabel('c', fontsize=20)
plt.ylabel('N', fontsize=20)

plt.tight_layout()

if save == True:
	#plt.savefig('histogram_values.eps')
	#plt.savefig('histogram_values.png')
	plt.savefig('histogram_c_values.png')

if show == True:
	plt.show()



