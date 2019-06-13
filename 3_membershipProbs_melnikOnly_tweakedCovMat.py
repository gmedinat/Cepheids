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

# IDEA: llenar con nans los campos que tengan X = 0 and e_X = 0.
# luego crear un c distinto dependiendo de que valores esten disponibles.

save = False
show = False
verbose = False

#print chi2.sf(0.51, 2)*100.
#print stats.chi2.cdf(0.51, 2), 1 - stats.chi2.cdf(0.51, 2)
#sys.exit(0)

#catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior'
#catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_05-04-19_dupRemoved_wPrior.csv'

#catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_RVMelnikOnly.csv'
#catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_RVMelnikOnly_ErrsTweak.csv'
#catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_RVMelnikOnly_SELECTED_ErrsTweak.csv'
catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_27-05-19_dupRem_wPrior_RVMelnikOnly_ErrsTweak.csv'

#catalog_name = 'test_membership.csv'

#cols = ['id','coord_ra','coord_dec', 'object', 'visit', 'ccdnum', 'MJDObs', 'base_SdssShape_psf_xx', 'base_SdssShape_psf_yy', 'AIRMASS', 'SEEING', 'EXPTIME', 'base_PsfFlux_flux', 'base_PsfFlux_fluxSigma']
#df = pd.read_csv(indir+'/RR%s/'%file_n+archivo, usecols=cols)
df = pd.read_csv(catalog_name)

#df = df[['id','coord_ra','coord_dec', 'RA_deg', 'DEC_deg', 'object', 'visit', 'ccdnum', 'MJDObs', 'base_SdssShape_psf_xx', 'base_SdssShape_psf_yy', 'AIRMASS', 'SEEING', 'EXPTIME', 'base_PsfFlux_flux', 'base_PsfFlux_fluxSigma']]#1

# Covariance
def cov(x, y):
    xbar, ybar = x.mean(), y.mean()
    return np.sum((x - xbar)*(y - ybar))/(len(x) - 1)

# number of degrees of freedom


#df['e_Dist'] = df['Dist']*0.2 # 20% of the distance assumed
#df['parallax_cluster'] = 1000./df['Dist']
#df['e_parallax_cluster'] = 1000.*df['e_Dist']/(df['Dist']*df['Dist']) 

parallax_gaia_zp = -0.03

# add columns with differences
# delta                 cluster   -   cepheid
df['parallax_diff'] = df['parallax_cluster_withCG'] - (df['parallax']+parallax_gaia_zp)   # 	, 	#cep_par 		 = cep_par - parallax_gaia_zp
df['vr_diff']       = df['RV_withCarrera19_CARRERA'] - df['HRV_melnik']
df['pmRA_diff']     = df['pmRA_withCG'] - df['pmra']
df['pmDEC_diff']    = df['pmDEC_withCG'] - df['pmdec']
#df['FeH_diff']      = df['[]'] - df['mag_2']
#df['age_diff']      = df['Age_D02'] - df['mag_2'] # or logt_K13

df['dfreed']    = 0.
df['c']      = 0.
df['p_val']   = 0.
df['P_B_A'] = 0.

df['cluster_cepheid'] = ""

n_nans = 0

largo = len(df)

i_0 = 0
i_1 = 0
i_2 = 0
i_3 = 0
i_4 = 0

i = 0
#i = 6200
#i = 2397
#i = 13421
while i < largo:
	#print i
	if (i+1)%200==0:
		print '\n\t', i+1, '/', largo

	subdf = df.iloc[i]

	#x_df = subdf[['parallax_diff','vr_diff','pmRA_diff','pmDEC_diff']]
	
	#x = x_df.as_matrix()
	#x = x.reshape(len(x), 1)
	##########xT = np.transpose(x)
	#xT = x.T

	cep_par    = subdf['parallax'] 
	ecep_par   = subdf['parallax_error']
	cep_rv     = subdf['HRV_melnik']
	ecep_rv    = subdf['e_HRV_melnik']
	cep_pmra   = subdf['pmra']
	ecep_pmra  = subdf['pmra_error']
	cep_pmdec  = subdf['pmdec']
	ecep_pmdec = subdf['pmdec_error']
	
	oc_par    = subdf['parallax_cluster_withCG']
	eoc_par   = subdf['e_parallax_cluster_withCG']
	#oc_rv     = subdf['RV']
	#eoc_rv    = subdf['e_RV']
	oc_rv     = subdf['RV_withCarrera19_CARRERA']
	eoc_rv    = subdf['e_RV_withCarrera19_CARRERA']
	oc_pmra   = subdf['pmRA_withCG']
	eoc_pmra  = subdf['e_pmRA_withCG']
	oc_pmdec  = subdf['pmDEC_withCG']
	eoc_pmdec = subdf['e_pmDEC_withCG']

	cep_par 		 = cep_par - parallax_gaia_zp

	sigma_pm_systematic = 0.035

	ruwe      = subdf['ruwe']	
	if math.isnan(ruwe):
		ruwe = 22 # the maximum value within the catalog is 21.76

	if ruwe <= 1.4:
		ruwe = 1.0

	#print '\n----------------------------------------------------------------------------------'
	#print 'Cepheids: par, e_par, rv, e_rv, pmra, e_pmra, pmdec, e_pmdec --> ', cep_par, ecep_par, cep_rv, ecep_rv, cep_pmra, ecep_pmra, cep_pmdec, ecep_pmdec
	#print 'OC: par, e_par, rv, e_rv, pmra, e_pmra, pmdec, e_pmdec --> ', oc_par, eoc_par, oc_rv, eoc_rv, oc_pmra, eoc_pmra, oc_pmdec, eoc_pmdec

	#print (cep_par == 0), (ecep_par == 0)
	#print (cep_par, ecep_par == (0,0))
	#print (cep_par, ecep_par == (0,0)) or (oc_par,eoc_par == (0,0))

	if (cep_par == 0 and ecep_par == 0) or (oc_par == 0 and eoc_par == 0): # as==0
		
		if (cep_rv == 0 and ecep_rv == 0) or (oc_rv == 0 and eoc_rv == 0):	# bs==0	 

			if (cep_pmra == 0 and ecep_pmra == 0) or (oc_pmra == 0 and eoc_pmra == 0): # cs==0
		
				if (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0
					dof = 0

					p_val = 1
					P_B_A = 1. - p_val
					i_0 = i_0+1

					
				else:
					dof = 1 # only d !=0
				
					x = subdf[['pmDEC_diff']]
	
					# covariance matrix cluster    pmDEC	
					M1 = subdf['e_pmDEC_withCG']**2

					# covariance matrix cepheid
					M2 = (ruwe**2)*subdf['pmdec_error']**2 + sigma_pm_systematic**2

					# sum of covariance matrices
					Matrix = M1+M2

					# inverse of the sum of covariance matrices

					MatrixInv = 1./Matrix

					ca = x*MatrixInv*x				
					c = np.squeeze(np.asarray(ca))	

					if i_1 == 0:
						c_distribution_1dof = c
						i_1 = i_1+1
					else: 
						c_distribution_1dof = np.append(c_distribution_1dof, c)
						i_1 = i_1+1

					df.dfreed.iloc[i]    = dof
					df.c.iloc[i]    = c
					#df.p(c).iloc[i]    = 0.
					df.P_B_A.iloc[i]  = chi2.sf(float(c), dof) 
					df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof) 

					 

			elif (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0			
				dof = 1 # only c!=0

				x = subdf[['pmRA_diff']]

				# covariance matrix cluster    pmRA	
				M1 = subdf['e_pmRA_withCG']**2
				# covariance matrix cepheid
				M2 = (ruwe**2)*subdf['pmra_error']**2 + sigma_pm_systematic**2

				# sum of covariance matrices
				Matrix = M1+M2

				# inverse of the sum of covariance matrices

				MatrixInv = 1./Matrix

				ca = x*MatrixInv*x				
				c = np.squeeze(np.asarray(ca))	

				if i_1 == 0:
					c_distribution_1dof = c
					i_1 = i_1+1
				else: 
					c_distribution_1dof = np.append(c_distribution_1dof, c)
					i_1 = i_1+1

				df.dfreed.iloc[i]    = dof
				df.c.iloc[i]    = c
				#df.p(c).iloc[i]    = 0.
				df.P_B_A.iloc[i]  = chi2.sf(float(c), dof)  
				df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)


			else:
				dof = 2 # c and d != 0

				x_df = subdf[['pmRA_diff','pmDEC_diff']]
	
				x = x_df.as_matrix()
				x = x.reshape(len(x), 1)
				xT = x.T

				# covariance matrix cluster   pmRA  pmDEC	
				M1 = np.matrix([[subdf['e_pmRA_withCG']**2, 0],
						[0, subdf['e_pmDEC_withCG']**2]]) 

				cov_pmra_pmdec = subdf['pmra_pmdec_corr']*subdf['pmra_error']*subdf['pmdec_error']

				# covariance matrix cepheid
				M2 = (ruwe**2)*np.matrix([[subdf['pmra_error']**2,   cov_pmra_pmdec],
						                [cov_pmra_pmdec,             subdf['pmdec_error']**2]])
				M2 = M2 + np.matrix([[sigma_pm_systematic**2,   0],
						                [0,             sigma_pm_systematic**2]])


				# sum of covariance matrices
				Matrix = M1+M2

				# inverse of the sum of covariance matrices

				MatrixInv = Matrix.getI()

				ca = xT.dot(MatrixInv.dot(x))  
				c = np.squeeze(np.asarray(ca))

				if i_2 == 0:
					c_distribution_2dof = c
					i_2 = i_2+1
				else: 
					c_distribution_2dof = np.append(c_distribution_2dof, c)
					i_2 = i_2+1

				df.dfreed.iloc[i]    = dof
				df.c.iloc[i]    = c
				#df.p(c).iloc[i]    = 0.
				df.P_B_A.iloc[i]  = chi2.sf(float(c), dof)  
				df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)


		elif (cep_pmra == 0 and ecep_pmra == 0) or (oc_pmra == 0 and eoc_pmra == 0):
			
			if (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0

				dof = 1 # only b != 0

				x = subdf[['rv_diff']]

				# covariance matrix cluster    RV	
				M1 = subdf['e_RV_withCarrera19_CARRERA']**2
				# covariance matrix cepheid
				M2 = (ruwe**2)*subdf['e_HRV_melnik']**2

				# sum of covariance matrices
				Matrix = M1+M2

				# inverse of the sum of covariance matrices
				MatrixInv = 1./Matrix

				ca = x*MatrixInv*x				
				c = np.squeeze(np.asarray(ca))	

				if i_1 == 0:
					c_distribution_1dof = c
					i_1 = i_1+1
				else: 
					c_distribution_1dof = np.append(c_distribution_1dof, c)
					i_1 = i_1+1

				df.dfreed.iloc[i]    = dof
				df.c.iloc[i]    = c
				#df.p(c).iloc[i]    = 0.
				df.P_B_A.iloc[i]  = chi2.sf(float(c), dof)  
				df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)

			else:
				dof = 2 # b and d != 0	

				x_df = subdf[['vr_diff','pmDEC_diff']]
	
				x = x_df.as_matrix()
				x = x.reshape(len(x), 1)
				xT = x.T

				# covariance matrix cluster    RV   pmDEC	
				M1 = np.matrix([[subdf['e_RV_withCarrera19_CARRERA']**2, 0],
						[0, subdf['e_pmDEC_withCG']**2]]) 

				cov_rv_pmdec   = 0.

				# covariance matrix cepheid
				M2 = (ruwe**2)*np.matrix([[subdf['e_HRV_melnik']**2, cov_rv_pmdec],
						[cov_rv_pmdec, subdf['pmdec_error']**2]])
				M2 = M2 + np.matrix([[0,   0],
						                [0,             sigma_pm_systematic**2]])

				# sum of covariance matrices
				Matrix = M1+M2
		
				# inverse of the sum of covariance matrices
		
				MatrixInv = Matrix.getI()
		
				ca = xT.dot(MatrixInv.dot(x))  
				c = np.squeeze(np.asarray(ca))
		
				if i_2 == 0:
					c_distribution_2dof = c
					i_2 = i_2+1
				else: 
					c_distribution_2dof = np.append(c_distribution_2dof, c)
					i_2 = i_2+1

				df.dfreed.iloc[i]    = dof
				df.c.iloc[i]    = c
				#df.p(c).iloc[i]    = 0.
				df.P_B_A.iloc[i]  = chi2.sf(float(c), dof) 
				df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)
		 

		elif (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0
			dof = 2 # b and c != 0
		
			x_df = subdf[['vr_diff','pmRA_diff']]
	
			x = x_df.as_matrix()
			x = x.reshape(len(x), 1)
			xT = x.T

			# covariance matrix cluster    RV   pmRA 
			M1 = np.matrix([[subdf['e_RV_withCarrera19_CARRERA']**2, 0],
					[0, subdf['e_pmRA_withCG']**2]]) 

			cov_rv_pmra    = 0.
				
			# covariance matrix cepheid
			M2 = (ruwe**2)*np.matrix([[subdf['e_HRV_melnik']**2, cov_rv_pmra],
					[cov_rv_pmra, subdf['pmra_error']**2]])

			M2 = M2 + np.matrix([[0,   0],
						                [0,       sigma_pm_systematic**2]])			

			# sum of covariance matrices
			Matrix = M1+M2
	
			# inverse of the sum of covariance matrices
	
			MatrixInv = Matrix.getI()
	
			ca = xT.dot(MatrixInv.dot(x))  
			c = np.squeeze(np.asarray(ca))
	
			if i_2 == 0:
				c_distribution_2dof = c
				i_2 = i_2+1
			else: 
				c_distribution_2dof = np.append(c_distribution_2dof, c)
				i_2 = i_2+1
		
			df.dfreed.iloc[i]    = dof
			df.c.iloc[i]    = c
			#df.p(c).iloc[i]    = 0.
			df.P_B_A.iloc[i]  = chi2.sf(float(c), dof) 			 
			df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)


		else:
			dof = 3 # b, c and d != 0

			x_df = subdf[['vr_diff','pmRA_diff','pmDEC_diff']]
	
			x = x_df.as_matrix()
			x = x.reshape(len(x), 1)
			xT = x.T


			# covariance matrix cluster    RV   pmRA  pmDEC	
			M1 = np.matrix([[subdf['e_RV_withCarrera19_CARRERA']**2, 0, 0],
					[0, subdf['e_pmRA_withCG']**2, 0],
					[0, 0, subdf['e_pmDEC_withCG']**2]]) 

			cov_rv_pmra    = 0.
			cov_rv_pmdec   = 0.
			cov_pmra_pmdec = subdf['pmra_pmdec_corr']*subdf['pmra_error']*subdf['pmdec_error']
	
			# covariance matrix cepheid
			M2 = (ruwe**2)*np.matrix([[subdf['e_HRV_melnik']**2, cov_rv_pmra, cov_rv_pmdec],
					[cov_rv_pmra, subdf['pmra_error']**2, cov_pmra_pmdec],
					[cov_rv_pmdec, cov_pmra_pmdec, subdf['pmdec_error']**2]])

			M2 = M2 + np.matrix([[0, 0, 0],
					[0, sigma_pm_systematic**2, 0],
					[0,0, sigma_pm_systematic**2]])

			# sum of covariance matrices
			Matrix = M1+M2
	
			# inverse of the sum of covariance matrices
	
			MatrixInv = Matrix.getI()
	
			ca = xT.dot(MatrixInv.dot(x))  
			c = np.squeeze(np.asarray(ca))

			if i_3 == 0:
				c_distribution_3dof = c
				i_3 = i_3+1
			else: 
				c_distribution_3dof = np.append(c_distribution_3dof, c)
				i_3 = i_3+1

			df.dfreed.iloc[i]    = dof
			df.c.iloc[i]    = c
			#df.p(c).iloc[i]    = 0.
			df.P_B_A.iloc[i]  = chi2.sf(float(c), dof)  
			df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)


	elif (cep_rv == 0 and ecep_rv == 0) or (oc_rv == 0 and eoc_rv == 0): # bs==0
		
		if (cep_pmra == 0 and ecep_pmra == 0) or (oc_pmra == 0 and eoc_pmra == 0):

			if (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0
				dof = 1 # only a != 0

				x = subdf[['parallax_diff']]

				# covariance matrix cluster    parallax	
				M1 = subdf['e_parallax_cluster_withCG']**2

				# covariance matrix cepheid
				M2 = (ruwe**2)*subdf['parallax_error']**2

				# sum of covariance matrices
				Matrix = M1+M2
				# inverse of the sum of covariance matrices

				MatrixInv = 1./Matrix

				ca = x*MatrixInv*x				
				c = np.squeeze(np.asarray(ca))	

				if i_1 == 0:
					c_distribution_1dof = c
					i_1 = i_1+1
				else: 
					c_distribution_1dof = np.append(c_distribution_1dof, c)
					i_1 = i_1+1

				df.dfreed.iloc[i]    = dof
				df.c.iloc[i]    = c
				#df.p(c).iloc[i]    = 0.
				df.P_B_A.iloc[i]  = chi2.sf(float(c), dof)  
				df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)

			else:
				dof = 2 # a and d != 0


				x_df = subdf[['parallax_diff','pmDEC_diff']]
	
				x = x_df.as_matrix()
				x = x.reshape(len(x), 1)
				xT = x.T

				# covariance matrix cluster    parallax   pmDEC	
				M1 = np.matrix([[subdf['e_parallax_cluster_withCG']**2, 0], 
						[0, subdf['e_pmDEC_withCG']**2]]) 

				cov_par_pmdec  = subdf['parallax_pmdec_corr']*subdf['parallax_error']*subdf['pmdec_error']
		
				# covariance matrix cepheid
				M2 = (ruwe**2)*np.matrix([[subdf['parallax_error']**2, cov_par_pmdec], 
						[cov_par_pmdec, subdf['pmdec_error']**2]])
				M2 = M2 + np.matrix([[0, 0],
					[0, sigma_pm_systematic**2]])

				# sum of covariance matrices
				Matrix = M1+M2

				# inverse of the sum of covariance matrices

				MatrixInv = Matrix.getI()

				ca = xT.dot(MatrixInv.dot(x))  
				c = np.squeeze(np.asarray(ca))
		
				if i_2 == 0:
					c_distribution_2dof = c
					i_2 = i_2+1
				else: 
					c_distribution_2dof = np.append(c_distribution_2dof, c)
					i_2 = i_2+1
	
				df.dfreed.iloc[i]    = dof
				df.c.iloc[i]    = c
				#df.p(c).iloc[i]    = 0.
				df.P_B_A.iloc[i]  = chi2.sf(float(c), dof) 			
				df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)	 


		elif (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0
			dof = 2 # a and c != 0


			x_df = subdf[['parallax_diff','pmRA_diff']]
	
			x = x_df.as_matrix()
			x = x.reshape(len(x), 1)
			xT = x.T

			# covariance matrix cluster    parallax   pmRA 
			M1 = np.matrix([[subdf['e_parallax_cluster_withCG']**2, 0], 
					[0, subdf['e_pmRA_withCG']**2]]) 
	
			cov_par_pmra   = subdf['parallax_pmra_corr']*subdf['parallax_error']*subdf['pmra_error']
	
			# covariance matrix cepheid
			M2 = (ruwe**2)*np.matrix([[subdf['parallax_error']**2, cov_par_pmra], 
					[cov_par_pmra, subdf['pmra_error']**2]])

			M2 = M2 + np.matrix([[0, 0],
					[0, sigma_pm_systematic**2]])
	
			# sum of covariance matrices
			Matrix = M1+M2
	
			# inverse of the sum of covariance matrices
	
			MatrixInv = Matrix.getI()
	
			ca = xT.dot(MatrixInv.dot(x))  
			c = np.squeeze(np.asarray(ca))
	
			if i_2 == 0:
				c_distribution_2dof = c
				i_2 = i_2+1
			else: 
				c_distribution_2dof = np.append(c_distribution_2dof, c)
				i_2 = i_2+1

			df.dfreed.iloc[i]    = dof
			df.c.iloc[i]    = c
			#df.p(c).iloc[i]    = 0.
			df.P_B_A.iloc[i]  = chi2.sf(float(c), dof) 
			df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)
			 

		else:
			dof = 3 # a, c and d != 0

			x_df = subdf[['parallax_diff','pmRA_diff','pmDEC_diff']]
	
			x = x_df.as_matrix()
			x = x.reshape(len(x), 1)
			xT = x.T

			# covariance matrix cluster    parallax  pmRA  pmDEC	
			M1 = np.matrix([[subdf['e_parallax_cluster_withCG']**2, 0, 0], 
					[0, subdf['e_pmRA_withCG']**2, 0],
					[0, 0, subdf['e_pmDEC_withCG']**2]]) 

			cov_par_pmra   = subdf['parallax_pmra_corr']*subdf['parallax_error']*subdf['pmra_error']
			cov_par_pmdec  = subdf['parallax_pmdec_corr']*subdf['parallax_error']*subdf['pmdec_error']
			cov_pmra_pmdec = subdf['pmra_pmdec_corr']*subdf['pmra_error']*subdf['pmdec_error']

			# covariance matrix cepheid
			M2 = (ruwe**2)*np.matrix([[subdf['parallax_error']**2, cov_par_pmra, cov_par_pmdec], 
					[cov_par_pmra, subdf['pmra_error']**2, cov_pmra_pmdec],
					[cov_par_pmdec, cov_pmra_pmdec, subdf['pmdec_error']**2]])

			M2 = M2 + np.matrix([[0, 0, 0],
					[0, sigma_pm_systematic**2, 0],
					[0,0, sigma_pm_systematic**2]])

			M2 = M2 + np.matrix([	[subdf['parallax']*0. + 0, 0, 0], 
							 		[0, subdf['pmra']*0. + 0, 0],
							 		[0, 0, subdf['pmdec']*0. + 0]		])

			# sum of covariance matrices
			Matrix = M1+M2

			# inverse of the sum of covariance matrices
	
			#print subdf['Cluster_name'], subdf['source_id_1']
			#print M1
			#print M2

			MatrixInv = Matrix.getI()
	
			ca = xT.dot(MatrixInv.dot(x))  
			c = np.squeeze(np.asarray(ca))

			if i_3 == 0:
				c_distribution_3dof = c
				i_3 = i_3+1
			else: 
				c_distribution_3dof = np.append(c_distribution_3dof, c)
				i_3 = i_3+1

			df.dfreed.iloc[i]    = dof
			df.c.iloc[i]    = c
			#df.p(c).iloc[i]    = 0.
			df.P_B_A.iloc[i]  = chi2.sf(float(c), dof)
			df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)

			#print subdf[['Cluster_name']], subdf[['Var_Name']]
			#print '\n', M1, '\n'
			#print M2, '\n'
			#print Matrix, '\n'
			#print MatrixInv, '\n'
			#print x, '\n'
			#print c, '\n'
			#print 'pba = ', chi2.sf(float(c), dof) 
			#print 'pval = ', 1. - chi2.sf(float(c), dof)
			#sys.exit(0)


	elif (cep_pmra == 0 and ecep_pmra == 0) or (oc_pmra == 0 and eoc_pmra == 0):
		if (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0
			dof = 2 # a and b != 0

			x_df = subdf[['parallax_diff','vr_diff']]
	
			x = x_df.as_matrix()
			x = x.reshape(len(x), 1)
			xT = x.T

			# covariance matrix cluster    parallax   RV
			M1 = np.matrix([[subdf['e_parallax_cluster_withCG']**2, 0], 
					[0, subdf['e_RV_withCarrera19_CARRERA']**2]]) 

			cov_par_rv     = 0.

			# covariance matrix cepheid
			M2 = (ruwe**2)*np.matrix([[subdf['parallax_error']**2, cov_par_rv], 
					[cov_par_rv, subdf['e_HRV_melnik']**2]])
	
	
			# sum of covariance matrices
			Matrix = M1+M2
	
			# inverse of the sum of covariance matrices
	
			MatrixInv = Matrix.getI()
	
			ca = xT.dot(MatrixInv.dot(x))  
			c = np.squeeze(np.asarray(ca))
	
			if i_2 == 0:
				c_distribution_2dof = c
				i_2 = i_2+1
			else: 
				c_distribution_2dof = np.append(c_distribution_2dof, c)
				i_2 = i_2+1

			df.dfreed.iloc[i]    = dof
			df.c.iloc[i]    = c
			#df.p(c).iloc[i]    = 0.
			df.P_B_A.iloc[i]  = chi2.sf(float(c), dof) 
			df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)


		else:
			dof = 3 # a, b and d != 0

			x_df = subdf[['parallax_diff','vr_diff','pmDEC_diff']]
	
			x = x_df.as_matrix()
			x = x.reshape(len(x), 1)
			xT = x.T

			# covariance matrix cluster    parallax   RV  pmDEC	
			M1 = np.matrix([[subdf['e_parallax_cluster_withCG']**2, 0, 0], 
					[0, subdf['e_RV_withCarrera19_CARRERA']**2, 0],
					[0, 0, subdf['e_pmDEC_withCG']**2]]) 

			cov_par_rv     = 0.
			cov_par_pmdec  = subdf['parallax_pmdec_corr']*subdf['parallax_error']*subdf['pmdec_error']
			cov_rv_pmdec   = 0.

			# covariance matrix cepheid
			M2 = (ruwe**2)*np.matrix([[subdf['parallax_error']**2, cov_par_rv, cov_par_pmdec], 
					[cov_par_rv, subdf['e_HRV_melnik']**2, cov_rv_pmdec],
					[cov_par_pmdec, cov_rv_pmdec, subdf['pmdec_error']**2]])

			M2 = M2 + np.matrix([[0, 0, 0],
					[0, 0, 0],
					[0,0, sigma_pm_systematic**2]])

			# sum of covariance matrices
			Matrix = M1+M2
	
			# inverse of the sum of covariance matrices

			MatrixInv = Matrix.getI()

			ca = xT.dot(MatrixInv.dot(x))  
			c = np.squeeze(np.asarray(ca))

			if i_3 == 0:
				c_distribution_3dof = c
				i_3 = i_3+1
			else: 
				c_distribution_3dof = np.append(c_distribution_3dof, c)
				i_3 = i_3+1

			df.dfreed.iloc[i]    = dof
			df.c.iloc[i]    = c
			#df.p(c).iloc[i]    = 0.
			df.P_B_A.iloc[i]  = chi2.sf(float(c), dof) 
			df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)
			 

	elif (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0
		dof = 3 # a, b and c != 0

		x_df = subdf[['parallax_diff','vr_diff','pmRA_diff']]
	
		x = x_df.as_matrix()
		x = x.reshape(len(x), 1)
		xT = x.T

		# covariance matrix cluster    parallax   RV   pmRA
		M1 = np.matrix([[subdf['e_parallax_cluster_withCG']**2, 0, 0], 
				[0, subdf['e_RV_withCarrera19_CARRERA']**2, 0],
				[0, 0, subdf['e_pmRA_withCG']**2]]) 

		cov_par_pmra   = subdf['parallax_pmra_corr']*subdf['parallax_error']*subdf['pmra_error']
		cov_par_rv     = 0.
		cov_rv_pmra    = 0.

		# covariance matrix cepheid
		M2 = (ruwe**2)*np.matrix([[subdf['parallax_error']**2, cov_par_rv, cov_par_pmra], 
				[cov_par_rv, subdf['e_HRV_melnik']**2, cov_rv_pmra],
				[cov_par_pmra, cov_rv_pmra, subdf['pmra_error']**2]])

		M2 = M2 + np.matrix([[0, 0, 0],
					[0, 0, 0],
					[0, 0, sigma_pm_systematic**2]])

		# sum of covariance matrices
		Matrix = M1+M2

		# inverse of the sum of covariance matrices

		MatrixInv = Matrix.getI()

		ca = xT.dot(MatrixInv.dot(x))  
		c = np.squeeze(np.asarray(ca))

		if i_3 == 0:
			c_distribution_3dof = c
			i_3 = i_3+1
		else: 
			c_distribution_3dof = np.append(c_distribution_3dof, c)
			i_3 = i_3+1
		
		df.dfreed.iloc[i]    = dof
		df.c.iloc[i]    = c
		#df.p(c).iloc[i]    = 0.
		df.P_B_A.iloc[i]  = chi2.sf(float(c), dof) 
		df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)


	else: # everything is != 0
		dof = 4 # a,b,c and d != 0

		x_df = subdf[['parallax_diff','vr_diff','pmRA_diff','pmDEC_diff']]
	
		#print x_df

		x = x_df.as_matrix()
		x = x.reshape(len(x), 1)
		xT = x.T

		# covariance matrix cluster    parallax   RV   pmRA  pmDEC	
		M1 = np.matrix([[subdf['e_parallax_cluster_withCG']**2, 0, 0, 0], 
				[0, subdf['e_RV_withCarrera19_CARRERA']**2, 0, 0],
				[0, 0, subdf['e_pmRA_withCG']**2, 0],
				[0, 0, 0, subdf['e_pmDEC_withCG']**2]]) 

		cov_par_pmra   = subdf['parallax_pmra_corr']*subdf['parallax_error']*subdf['pmra_error']
		cov_par_rv     = 0.
		cov_par_pmdec  = subdf['parallax_pmdec_corr']*subdf['parallax_error']*subdf['pmdec_error']
		cov_rv_pmra    = 0.
		cov_rv_pmdec   = 0.
		cov_pmra_pmdec = subdf['pmra_pmdec_corr']*subdf['pmra_error']*subdf['pmdec_error']

		# covariance matrix cepheid
		M2 = (ruwe**2)*np.matrix([[subdf['parallax_error']**2, cov_par_rv, cov_par_pmra, cov_par_pmdec], 
				[cov_par_rv, subdf['e_HRV_melnik']**2, cov_rv_pmra, cov_rv_pmdec],
				[cov_par_pmra, cov_rv_pmra, subdf['pmra_error']**2, cov_pmra_pmdec],
				[cov_par_pmdec, cov_rv_pmdec, cov_pmra_pmdec, subdf['pmdec_error']**2]])

		M2 = M2 + np.matrix([[0,0,0,0], 
				[0,0,0,0],
				[0,0, sigma_pm_systematic**2, 0],
				[0,0,0, sigma_pm_systematic**2]])

		#np.set_printoptions(formatter={"float_kind": lambda x: "%.4g" % x})
		#print '\n cluster values: \n', 	subdf[['parallax_cluster_withCG', 'RV_withCarrera19_CARRERA', 'pmRA_withCG', 'pmDEC_withCG']]
		#print '\n cepheid values: \n', subdf[['parallax','HRV_melnik','pmra','pmdec']]
		#print 'M2 without major tweak: \n', M2, '\n', '\n'

		M2 = np.matrix([[abs(subdf['parallax']*0.15) + 0., 0, 0, 0], 
							 [0, abs(subdf['HRV_melnik']*0.2) + 0, 0, 0],
							 [0, 0, abs(subdf['pmra']*0.15) + 0.0, 0],
							 [0, 0, 0, abs(subdf['pmdec']*0.15) + 0.0]])

		# sum of covariance matrices
		Matrix = M1+M2

		# inverse of the sum of covariance matrices

		MatrixInv = Matrix.getI()

		ca = xT.dot(MatrixInv.dot(x))  
		c = np.squeeze(np.asarray(ca))

		if i_4 == 0:
			c_distribution_4dof = c
			i_4 = i_4+1
		else: 
			c_distribution_4dof = np.append(c_distribution_4dof, c)
			i_4 = i_4+1
	
		df.dfreed.iloc[i]    = dof
		df.c.iloc[i]    = c
		#df.p(c).iloc[i]    = 0.
		df.P_B_A.iloc[i]  = chi2.sf(float(c), dof) 
		df.p_val.iloc[i]  = 1. - chi2.sf(float(c), dof)

		#print subdf[['Cluster_name']], subdf[['Var_Name']]
		#print '\n', M1, '\n'
		#print 'M2 with major tweak: \n', M2, '\n'
		#print Matrix, '\n'
		#print MatrixInv, '\n'
		#print x, '\n'
		#print c, '\n'
		#print 'pba = ', chi2.sf(float(c), dof) 
		#print 'pval = ', 1. - chi2.sf(float(c), dof)
		#sys.exit(0)

	if pd.isnull(subdf['Var_Name']):  # == ""
		#print 'entro1'
		#df.cluster_cepheid.iloc[i]  = subdf['Cluster_name'] + '_' + '%.3f_%.3f'%(subdf['ra'], subdf['dec'])
		df.cluster_cepheid.iloc[i]  = subdf['Cluster_name'] + '_' 
	#if subdf['Var_Name'] != "": # != ""
	else:
		#print 'entro2'
		df.cluster_cepheid.iloc[i]  = subdf['Cluster_name'] + '_' + subdf['Var_Name'].replace(" ","_")

	del subdf

	i = i+1

df['P_A*p_val'] = df['P_A']*df['p_val']
df['P_A*P_B_A'] = df['P_A']*df['P_B_A']

df.to_csv( './'+catalog_name.replace('.csv','_covMatTweak_wProbs.csv'), index=False )

del df



#print "1 degree of freedom:", c_distribution_1dof.shape,'\n', c_distribution_1dof
print "2 degree of freedom:", c_distribution_2dof.shape,'\n', c_distribution_2dof
#print "3 degree of freedom:", c_distribution_3dof.shape,'\n', c_distribution_3dof
#print "4 degree of freedom:", c_distribution_4dof.shape,'\n', c_distribution_4dof

np.savetxt( "./c_distribution_1dof.dat", c_distribution_1dof, fmt='%s',delimiter="   ", header='c')
np.savetxt( "./c_distribution_2dof.dat", c_distribution_2dof, fmt='%s',delimiter="   ", header='c')
np.savetxt( "./c_distribution_3dof.dat", c_distribution_3dof, fmt='%s',delimiter="   ", header='c')
np.savetxt( "./c_distribution_4dof.dat", c_distribution_4dof, fmt='%s',delimiter="   ", header='c')

print '\n\ni_0 = ', i_0
print 'i_1 = ', i_1
print 'i_2 = ', i_2
print 'i_3 = ', i_3
print 'i_4 = ', i_4, '\n'

# en este punto se deberian tener los sgtes arreglos:
# c_distribution_0dof, c_distribution_1dof, c_distribution_2dof, c_distribution_3dof, c_distribution_4dof  
# que distribuyen como
# chi2 con 0, 1, 2, 3 y 4 grados de libertad, respectivamente.


c_distribution_1dof = c_distribution_1dof.reshape(len(c_distribution_1dof), 1)
c_distribution_2dof = c_distribution_2dof.reshape(len(c_distribution_2dof), 1)
c_distribution_3dof = c_distribution_3dof.reshape(len(c_distribution_3dof), 1)
c_distribution_4dof = c_distribution_4dof.reshape(len(c_distribution_4dof), 1)

# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #

fig, ax = plt.subplots()
plt.hist(c_distribution_1dof, bins=30, density=True)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.xlabel('c1', fontsize=20)
plt.ylabel('N', fontsize=20)

plt.tight_layout()

if save == True:
	#plt.savefig('histogram_values_1dof.eps')
	#plt.savefig('histogram_values_1dof.png')
	plt.savefig('histogram_c_values_1dof.png')

#if show == True:
#	plt.show()

# ---------------------------------------------------------------------------------- #

fig, ax = plt.subplots()
plt.hist(c_distribution_2dof, bins=30, density=True)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.xlabel('c2', fontsize=20)
plt.ylabel('N', fontsize=20)

plt.tight_layout()

if save == True:
	#plt.savefig('histogram_values_2dof.eps')
	#plt.savefig('histogram_values_2dof.png')
	plt.savefig('histogram_c_values_2dof.png')

#if show == True:
#	plt.show()

# ---------------------------------------------------------------------------------- #

fig, ax = plt.subplots()
plt.hist(c_distribution_3dof, bins=30, density=True)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.xlabel('c3', fontsize=20)
plt.ylabel('N', fontsize=20)

plt.tight_layout()

if save == True:
	#plt.savefig('histogram_values_3dof.eps')
	#plt.savefig('histogram_values_3dof.png')
	plt.savefig('histogram_c_values_3dof.png')

#if show == True:
#	plt.show()

# ---------------------------------------------------------------------------------- #

fig, ax = plt.subplots()
#plt.hist(c_distribution_4dof, bins=60, density=True)
plt.hist(c_distribution_4dof, bins=60, density=False)

c4_mean = np.mean(c_distribution_4dof)
print c4_mean
#x4 = np.linspace(chi2.ppf(0.01, 4, scale=10, loc=100), chi2.ppf(0.9999, df, scale=10, loc=100), 100)
#ax.plot(x4, chi2.pdf(x4, df, scale=10, loc=100), 'g-', lw=5, alpha=0.6, label='Scale=10, loc=100, df=%i'%df)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.xlabel('c4', fontsize=20)
plt.ylabel('N', fontsize=20)

plt.tight_layout()

if save == True:
	#plt.savefig('histogram_values_4dof.eps')
	#plt.savefig('histogram_values_4dof.png')
	plt.savefig('histogram_c_values_4dof.png')

fig, ax = plt.subplots()
x4 = np.linspace(chi2.ppf(0.01, 4, loc=c4_mean), chi2.ppf(0.9999, 4, loc=c4_mean), 10000)
plt.plot(x4, chi2.pdf(x4, 4, loc=c4_mean), 'g-', lw=5, alpha=0.6, label='loc=%f, df=%i'%(c4_mean,4))

if show == True:
	plt.show()








