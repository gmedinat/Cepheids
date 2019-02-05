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

catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior'
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

i_0 = 0
i_1 = 0
i_2 = 0
i_3 = 0
i_4 = 0

i = 0
while i < largo:

	#print '\n\t', i+1, '/', largo

	subdf = df.iloc[i]

	#x_df = subdf[['parallax_diff','vr_diff','pmRA_diff','pmDEC_diff']]
	
	#x = x_df.as_matrix()
	#x = x.reshape(len(x), 1)
	##########xT = np.transpose(x)
	#xT = x.T

	cep_par    = subdf['parallax'] 
	ecep_par   = subdf['parallax_error']
	cep_rv     = subdf['radial_velocity']
	ecep_rv    = subdf['radial_velocity_error']
	cep_pmra   = subdf['pmra']
	ecep_pmra  = subdf['pmra_error']
	cep_pmdec  = subdf['pmdec']
	ecep_pmdec = subdf['pmdec_error']
	
	oc_par    = subdf['parallax_cluster']
	eoc_par   = subdf['e_parallax_cluster']
	oc_rv     = subdf['RV_2']
	eoc_rv    = subdf['e_RV']
	oc_pmra   = subdf['pmRA_2']
	eoc_pmra  = subdf['e_pmRA']
	oc_pmdec  = subdf['pmDE_2']
	eoc_pmdec = subdf['e_pmDE']

	print '\n----------------------------------------------------------------------------------'
	print 'Cepheids: par, e_par, rv, e_rv, pmra, e_pmra, pmdec, e_pmdec --> ', cep_par, ecep_par, cep_rv, ecep_rv, cep_pmra, ecep_pmra, cep_pmdec, ecep_pmdec
	print 'OC: par, e_par, rv, e_rv, pmra, e_pmra, pmdec, e_pmdec --> ', oc_par, eoc_par, oc_rv, eoc_rv, oc_pmra, eoc_pmra, oc_pmdec, eoc_pmdec

	print (cep_par == 0), (ecep_par == 0)
	print (cep_par, ecep_par == (0,0))
	print (cep_par, ecep_par == (0,0)) or (oc_par,eoc_par == (0,0))

	if (cep_par == 0 and ecep_par == 0) or (oc_par == 0 and eoc_par == 0): # as==0
		
		if (cep_rv == 0 and ecep_rv == 0) or (oc_rv == 0 and eoc_rv == 0):	# bs==0	 

			if (cep_pmra == 0 and ecep_pmra == 0) or (oc_pmra == 0 and eoc_pmra == 0): # cs==0
		
				if (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0
					dof = 0

					P_B_A = 1
					i_0 = i_0+1
					
				else:
					dof = 1 # only d !=0
				
					x = subdf[['pmDEC_diff']]
	
					# covariance matrix cluster    pmDEC	
					M1 = subdf['e_pmDE']**2

					# covariance matrix cepheid
					M2 = subdf['pmdec_error']**2

					# sum of covariance matrices
					Matrix = M1+M2

					# inverse of the sum of covariance matrices

					MatrixInv = 1./Matrix

					c = x*MatrixInv*x					

					if i_1 == 0:
						c_distribution_1dof = c
						i_1 = i_1+1
					else: 
						c_distribution_1dof = np.append(c_distribution_1dof, c)
						i_1 = i_1+1


			elif (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0			
				dof = 1 # only c!=0

				x = subdf[['pmRA_diff']]

				# covariance matrix cluster    pmRA	
				M1 = subdf['e_pmRA']**2
				# covariance matrix cepheid
				M2 = subdf['pmra_error']**2

				# sum of covariance matrices
				Matrix = M1+M2

				# inverse of the sum of covariance matrices

				MatrixInv = 1./Matrix

				c = x*MatrixInv*x					

				if i_1 == 0:
					c_distribution_1dof = c
					i_1 = i_1+1
				else: 
					c_distribution_1dof = np.append(c_distribution_1dof, c)
					i_1 = i_1+1


			else:
				dof = 2 # c and d != 0

				x_df = subdf[['pmRA_diff','pmDEC_diff']]
	
				x = x_df.as_matrix()
				x = x.reshape(len(x), 1)
				xT = x.T

				# covariance matrix cluster   pmRA  pmDEC	
				M1 = np.matrix([[subdf['e_pmRA']**2, 0],
						[0, subdf['e_pmDE']**2]]) 

				cov_pmra_pmdec = subdf['pmra_pmdec_corr']*subdf['pmra_error']*subdf['pmdec_error']

				# covariance matrix cepheid
				M2 = np.matrix([[subdf['pmra_error']**2, cov_pmra_pmdec],
						[cov_pmra_pmdec, subdf['pmdec_error']**2]])


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




		elif (cep_pmra == 0 and ecep_pmra == 0) or (oc_pmra == 0 and eoc_pmra == 0):
			
			if (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0

				dof = 1 # only b != 0

				x = subdf[['rv_diff']]

				# covariance matrix cluster    RV	
				M1 = subdf['e_RV']**2
				# covariance matrix cepheid
				M2 = subdf['radial_velocity_error']**2

				# sum of covariance matrices
				Matrix = M1+M2

				# inverse of the sum of covariance matrices
				MatrixInv = 1./Matrix

				c = x*MatrixInv*x					
				if i_1 == 0:
					c_distribution_1dof = c
					i_1 = i_1+1
				else: 
					c_distribution_1dof = np.append(c_distribution_1dof, c)
					i_1 = i_1+1

			else:
				dof = 2 # b and d != 0	

				x_df = subdf[['vr_diff','pmDEC_diff']]
	
				x = x_df.as_matrix()
				x = x.reshape(len(x), 1)
				xT = x.T

				# covariance matrix cluster    RV   pmDEC	
				M1 = np.matrix([[subdf['e_RV']**2, 0],
						[0, subdf['e_pmDE']**2]]) 

				cov_rv_pmdec   = 0.

				# covariance matrix cepheid
				M2 = np.matrix([[subdf['radial_velocity_error']**2, cov_rv_pmdec],
						[cov_rv_pmdec, subdf['pmdec_error']**2]])

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





		 

		elif (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0
			dof = 2 # b and c != 0
		
			x_df = subdf[['vr_diff','pmRA_diff']]
	
			x = x_df.as_matrix()
			x = x.reshape(len(x), 1)
			xT = x.T

			# covariance matrix cluster    RV   pmRA 
			M1 = np.matrix([[subdf['e_RV']**2, 0],
					[0, subdf['e_pmRA']**2]]) 

			cov_rv_pmra    = 0.
				
			# covariance matrix cepheid
			M2 = np.matrix([[subdf['radial_velocity_error']**2, cov_rv_pmra],
					[cov_rv_pmra, subdf['pmra_error']**2]])

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

		else:
			dof = 3 # b, c and d != 0

			x_df = subdf[['vr_diff','pmRA_diff','pmDEC_diff']]
	
			x = x_df.as_matrix()
			x = x.reshape(len(x), 1)
			xT = x.T


			# covariance matrix cluster    RV   pmRA  pmDEC	
			M1 = np.matrix([[subdf['e_RV']**2, 0, 0],
					[0, subdf['e_pmRA']**2, 0],
					[0, 0, subdf['e_pmDE']**2]]) 

			cov_rv_pmra    = 0.
			cov_rv_pmdec   = 0.
			cov_pmra_pmdec = subdf['pmra_pmdec_corr']*subdf['pmra_error']*subdf['pmdec_error']
	
			# covariance matrix cepheid
			M2 = np.matrix([[subdf['radial_velocity_error']**2, cov_rv_pmra, cov_rv_pmdec],
					[cov_rv_pmra, subdf['pmra_error']**2, cov_pmra_pmdec],
					[cov_rv_pmdec, cov_pmra_pmdec, subdf['pmdec_error']**2]])

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



	elif (cep_rv == 0 and ecep_rv == 0) or (oc_rv == 0 and eoc_rv == 0): # bs==0
		
		if (cep_pmra == 0 and ecep_pmra == 0) or (oc_pmra == 0 and eoc_pmra == 0):

			if (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0
				dof = 1 # only a != 0

				x = subdf[['parallax_diff']]

				# covariance matrix cluster    parallax	
				M1 = subdf['e_parallax_cluster']**2

				# covariance matrix cepheid
				M2 = subdf['parallax_error']**2

				# sum of covariance matrices
				Matrix = M1+M2
				# inverse of the sum of covariance matrices

				MatrixInv = 1./Matrix

				c = x*MatrixInv*x					

				if i_1 == 0:
					c_distribution_1dof = c
					i_1 = i_1+1
				else: 
					c_distribution_1dof = np.append(c_distribution_1dof, c)
					i_1 = i_1+1




			else:
				dof = 2 # a and d != 0


				x_df = subdf[['parallax_diff','pmDEC_diff']]
	
				x = x_df.as_matrix()
				x = x.reshape(len(x), 1)
				xT = x.T

				# covariance matrix cluster    parallax   pmDEC	
				M1 = np.matrix([[subdf['e_parallax_cluster']**2, 0], 
						[0, subdf['e_pmDE']**2]]) 

				cov_par_pmdec  = subdf['parallax_pmdec_corr']*subdf['parallax_error']*subdf['pmdec_error']
		
				# covariance matrix cepheid
				M2 = np.matrix([[subdf['parallax_error']**2, cov_par_pmdec], 
						[cov_par_pmdec, subdf['pmdec_error']**2]])

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


		elif (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0
			dof = 2 # a and c != 0


			x_df = subdf[['parallax_diff','pmRA_diff']]
	
			x = x_df.as_matrix()
			x = x.reshape(len(x), 1)
			xT = x.T

			# covariance matrix cluster    parallax   pmRA 
			M1 = np.matrix([[subdf['e_parallax_cluster']**2, 0], 
					[0, subdf['e_pmRA']**2]]) 
	
			cov_par_pmra   = subdf['parallax_pmra_corr']*subdf['parallax_error']*subdf['pmra_error']
	
			# covariance matrix cepheid
			M2 = np.matrix([[subdf['parallax_error']**2, cov_par_pmra], 
					[cov_par_pmra, subdf['pmra_error']**2]])
	
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

		else:
			dof = 3 # a, c and d != 0

			x_df = subdf[['parallax_diff','pmRA_diff','pmDEC_diff']]
	
			x = x_df.as_matrix()
			x = x.reshape(len(x), 1)
			xT = x.T

			# covariance matrix cluster    parallax  pmRA  pmDEC	
			M1 = np.matrix([[subdf['e_parallax_cluster']**2, 0, 0], 
					[0, subdf['e_pmRA']**2, 0],
					[0, 0, subdf['e_pmDE']**2]]) 

			cov_par_pmra   = subdf['parallax_pmra_corr']*subdf['parallax_error']*subdf['pmra_error']
			cov_par_pmdec  = subdf['parallax_pmdec_corr']*subdf['parallax_error']*subdf['pmdec_error']
			cov_pmra_pmdec = subdf['pmra_pmdec_corr']*subdf['pmra_error']*subdf['pmdec_error']

			# covariance matrix cepheid
			M2 = np.matrix([[subdf['parallax_error']**2, cov_par_pmra, cov_par_pmdec], 
					[cov_par_pmra, subdf['pmra_error']**2, cov_pmra_pmdec],
					[cov_par_pmdec, cov_pmra_pmdec, subdf['pmdec_error']**2]])

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

	

	elif (cep_pmra == 0 and ecep_pmra == 0) or (oc_pmra == 0 and eoc_pmra == 0):
		if (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0
			dof = 2 # a and b != 0

			x_df = subdf[['parallax_diff','vr_diff']]
	
			x = x_df.as_matrix()
			x = x.reshape(len(x), 1)
			xT = x.T

			# covariance matrix cluster    parallax   RV
			M1 = np.matrix([[subdf['e_parallax_cluster']**2, 0], 
					[0, subdf['e_RV']**2]]) 

			cov_par_rv     = 0.

			# covariance matrix cepheid
			M2 = np.matrix([[subdf['parallax_error']**2, cov_par_rv], 
					[cov_par_rv, subdf['radial_velocity_error']**2]])
	
	
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



		else:
			dof = 3 # a, b and d != 0

			x_df = subdf[['parallax_diff','vr_diff','pmDEC_diff']]
	
			x = x_df.as_matrix()
			x = x.reshape(len(x), 1)
			xT = x.T

			# covariance matrix cluster    parallax   RV  pmDEC	
			M1 = np.matrix([[subdf['e_parallax_cluster']**2, 0, 0], 
					[0, subdf['e_RV']**2, 0],
					[0, 0, subdf['e_pmDE']**2]]) 

			cov_par_rv     = 0.
			cov_par_pmdec  = subdf['parallax_pmdec_corr']*subdf['parallax_error']*subdf['pmdec_error']
			cov_rv_pmdec   = 0.

			# covariance matrix cepheid
			M2 = np.matrix([[subdf['parallax_error']**2, cov_par_rv, cov_par_pmdec], 
					[cov_par_rv, subdf['radial_velocity_error']**2, cov_rv_pmdec],
					[cov_par_pmdec, cov_rv_pmdec, subdf['pmdec_error']**2]])

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


	elif (cep_pmdec == 0 and ecep_pmdec == 0) or (oc_pmdec == 0 and eoc_pmdec == 0): # ds == 0
		dof = 3 # a, b and c != 0

		x_df = subdf[['parallax_diff','vr_diff','pmRA_diff']]
	
		x = x_df.as_matrix()
		x = x.reshape(len(x), 1)
		xT = x.T

		# covariance matrix cluster    parallax   RV   pmRA
		M1 = np.matrix([[subdf['e_parallax_cluster']**2, 0, 0], 
				[0, subdf['e_RV']**2, 0],
				[0, 0, subdf['e_pmRA']**2]]) 

		cov_par_pmra   = subdf['parallax_pmra_corr']*subdf['parallax_error']*subdf['pmra_error']
		cov_par_rv     = 0.
		cov_rv_pmra    = 0.

		# covariance matrix cepheid
		M2 = np.matrix([[subdf['parallax_error']**2, cov_par_rv, cov_par_pmra], 
				[cov_par_rv, subdf['radial_velocity_error']**2, cov_rv_pmra],
				[cov_par_pmra, cov_rv_pmra, subdf['pmra_error']**2]])

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



	else: # everything is != 0
		dof = 4 # a,b,c and d != 0

		x_df = subdf[['parallax_diff','vr_diff','pmRA_diff','pmDEC_diff']]
	
		x = x_df.as_matrix()
		x = x.reshape(len(x), 1)
		xT = x.T

		# covariance matrix cluster    parallax   RV   pmRA  pmDEC	
		M1 = np.matrix([[subdf['e_parallax_cluster']**2, 0, 0, 0], 
				[0, subdf['e_RV']**2, 0, 0],
				[0, 0, subdf['e_pmRA']**2, 0],
				[0, 0, 0, subdf['e_pmDE']**2]]) 

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


	del subdf

	i = i+1


del df





print '\n\ni_0 = ', i_0
print 'i_1 = ', i_1
print 'i_2 = ', i_2
print 'i_3 = ', i_3
print 'i_4 = ', i_4, '\n'



# en este punto se deberian tener los sgtes arreglos:
# c_distribution_0dof, c_distribution_1dof, c_distribution_2dof, c_distribution_3dof, c_distribution_4dof  
# que distribuyen como
# chi2 con 0, 1, 2, 3 y 4 grados de libertad, respectivamente.




# ---------------------------------------------------------------------------------- #






