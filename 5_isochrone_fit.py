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

#highprob_catalog_name = 'highProbTest_MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior_wProbs_wNames'
#df_highprob = pd.read_csv(highprob_catalog_name)

highprob_catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_wProbs.csv'
highprob_catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_27-05-19_dupRem_wPrior_RVMelnikOnly_ErrsTweak_covMatTweak_wProbs.csv'

df_highprob_0 = pd.read_csv(highprob_catalog_name)
#df_highprob = df_highprob_0.loc[df_highprob_0['p_val'] >= 0.05]
df_highprob = df_highprob_0.loc[df_highprob_0['P_B_A'] >= 0.05]

def closest_point(point, points):
    """ Find closest point from a list of points. """
    #return points[cdist([point], points).argmin()]
    return distance.cdist([point], [points], 'euclidean').argmin()

############## INCLUDE MAGNITUDE UNCERTAINTIES!!!!!!!!!!!!!
############# FILTER STARS BY MAGNITUDE, SOMEHOW!!!! 

verbose = False

# first try with 1 OC...
#cluster_name = 'Majaess_225'   # OJO CON LOS CLUSTERS QUE TIENEN DISTANCIAS = 0
#cluster_name = 'NGC_1901'
#cluster_name = 'NGC_6087'
#cluster_name = 'Alessi_24'
#cluster_name = 'Hogg_7'
#cluster_name = 'IC_2391'
#cluster_name = 'ASCC_63'
cluster_name = 'NGC_7235'
#cluster_name = 'Westerlund_2'

catalog_name = './Clusters/%s/%s_stars_noRVs.csv'%(cluster_name,cluster_name)
df = pd.read_csv(catalog_name)

index_ofthiscluster = df_highprob[df_highprob['Cluster_name']==cluster_name].index.values.astype(int)[0] # check where in the original highprob catalog is this cluster

EB_V_guess = df_highprob_0['E(B-V)'].iloc[index_ofthiscluster] # unique index
dist_guess = df_highprob_0['Dist'].iloc[index_ofthiscluster]

# this is just to record what is the closest point in the isochrone. Is gonna be changed in the next loop
df['closest_G']    = 0
df['closest_G_BP'] = 0
df['closest_G_RP'] = 0

g  = df['phot_g_mean_mag']
bp = df['phot_bp_mean_mag']
rp = df['phot_rp_mean_mag']

g_bp =  df['phot_g_mean_mag'] - df['phot_bp_mean_mag']
bp_rp =  df['phot_bp_mean_mag'] - df['phot_rp_mean_mag']

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------
# [M/H]=log(Z/Zo), with Zo=0.019   LOG 10 ?
# or  [M/H]=log(Z/X)-log(Z/X)o, with (Z/X)o=0.0207 and Y=0.2485+1.78Z for PARSEC tracks.
# linea creadora, centrada en valores conocidos
# only consider metallicities fe/h > -0.6   (typical of OCs?)
# if df_highprob['[Fe/H]'].iloc[index_ofthiscluster] != null   !!!!!

#Zo = 0.019
Zo = 0.0152 # updated solar content

# look for the closest value to that one in a list of available metallicities
# select range of metallicities to be checked
FeH_guess_ini = df_highprob_0['FeH'].iloc[index_ofthiscluster]
#Z_guess = 10**(df_highprob['[Fe/H]'].iloc[index_ofthiscluster])*Zo # IS IT LOG 10, RIGHT??????????????????!!!!!!!!!!!!!!!!!

if float(FeH_guess_ini) == 0 or FeH_guess_ini is None: # no info is available. Assume solar met. as initial guess
	Z_guess_ini = 0.0152
else:
	Z_guess_ini = 10.**df_highprob_0['FeH'].iloc[index_ofthiscluster]*Zo

#approximate Z_guess_ini to a valid value. (0.006, 0.007, 0.008, ..., 0.027, 0.028, 0.029) Round
Z_guess = round(Z_guess_ini,3)

if (Z_guess - 0.005) <= 0.006: # given by the limits in the isochrones I downloaded
	initial_Z = 0.006
	final_Z   = Z_guess + 0.005 	
elif (Z_guess + 0.005) >=  0.029: # given by the limits in the isochrones I downloaded
	initial_Z = Z_guess - 0.005
	final_Z   = 0.029
else:
	initial_Z = Z_guess - 0.005
	final_Z   = Z_guess + 0.005	
print FeH_guess_ini
print Z_guess_ini
print Z_guess
print initial_Z, final_Z

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------

min_redchi2 = 100000000000

initial_age = 710 # logAge x 100
final_age = 740 # logAge x 100

age_guess_ini = df_highprob_0['Age_D02'].iloc[index_ofthiscluster]
if age_guess_ini == 0:
	age_guess_ini = df_highprob_0['logt_K13'].iloc[index_ofthiscluster]

# note that, for example, 7.95 rounds to 8.0
if round(age_guess_ini,1) <= age_guess_ini:  #     7.9    7.9    8.0       7.9   7.92    8.0         8.0  8.0 8.1
	age_guess = int( round(round(age_guess_ini,1)*100))

else: # if age_guess < round(age_guess)     7.9   7.97  8.0    
	age_guess = int( round(round(age_guess_ini,1) - 0.1)*100)

if age_guess == 0:
	initial_age = 660	
	final_age = 700
	# FLAG IT!
elif age_guess <= 660:
	initial_age = 660	
	final_age = 700	
elif (age_guess - 40) <= 660:
	initial_age = 660	
	final_age = age_guess + 40
	# FLAG IT! 
elif age_guess >= 890:
	initial_age = 850	
	final_age = 890
	# FLAG IT!
elif (age_guess + 40) >= 890:
	initial_age = age_guess - 40	
	final_age = 890
	# FLAG IT!
else:
	initial_age = age_guess - 40  	
	final_age = age_guess + 40 

print age_guess_ini
print age_guess
print initial_age, final_age
sys.exit(0)	

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------

current_initial_age = initial_age
current_final_age = current_initial_age + 10

while current_initial_age < final_age: # this is to open the files. Each file has different metallicities, different ages in the range (current_initial_age,current_final_age)

	#output685-710.dat
	set_name = './Isochrones/step10/output%s-%s.dat'%(current_initial_age,current_final_age)
	table = ascii.read(set_name)
	#table = Table.read('testRR1_all_g.fits', format='fits')

	#Z_range_min = 0.0152
	#Z_range_max = 0.0152

	df_set = table.to_pandas()
	#df_set.rename(columns={'col1': 'Z', 'col2': 'log_age_yr', 'col3': 'M_ini','col4': 'M_act','col5': 'LogL/Lo','col6': 'LogTe','col7': 'LogG','col8': 'mbol','col9': 'G','col10': 'G_BP','col11': 'G_RP','col12': 'int_IMF','col13': 'stage'}, inplace=True)


	# output new interface: # Zini logAge Mini int_IMF  Mass   logL    logTe  logg  label   McoreTP C_O  period0 period1 pmode  
	#Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z 	 mbolmag  Gmag    G_BPmag  G_RPmag  B_Tmag  V_Tmag  Jmag    Hmag    Ksmag
	df_set.rename(columns={'col1':'Zini','col2':'log_age_yr','col3':'Mini','col4':'int_IMF','col5':'Mass','col6':'logL','col7':'logTe','col8':'logg','col9':'label','col10':'McoreTP','col11':'C_O','col12':'period0','col13':'period1','col14':'pmode', 'col15':'Mloss','col16':'tau1m','col17':'X','col18':'Y','col19':'Xc','col20':'Xn','col21':'Xo','col22':'Cexcess','col23':'Z','col24':'mbolmag','col25':'G','col26':'G_BP','col27':'G_RP','col28':'B_Tmag','col29':'V_Tmag','col30':'Jmag','col31':'Hmag','col32':'Ksmag'}, inplace=True)

	df_set = df_set[['Zini','log_age_yr','Mini','Mass','logL','logTe','logg','label','Mloss','Z', 'G','G_BP','G_RP']]#1

	#age_limit = 8.6989  # 500 Myr
	#age_limit = 8.8129  # 650 Myr
	#df_set = df_set.loc[ df_set['log_age_yr'] < age_limit ] # only consider isochrones with ages < log(500 Myr) = 8.6989   or ages < log(650 Myr) = 8.8129

	#df_set = pd.read_csv(set_name)

	diff_met = df_set.Zini.unique()
	diff_met = np.sort(diff_met)
	#print diff_met
	diff_met = diff_met[diff_met >= initial_Z]
	#print diff_met	
	diff_met = diff_met[diff_met <= final_Z]
	#print diff_met

	met_range_min = diff_met[0]
	met_range_max = diff_met[-1]

	len_diff_met = len(diff_met)

	index_met = 0

	while index_met < len_diff_met: # con esto se recorre el set completo de isocronas

		df_set_met = df_set.loc[df_set['Zini'] == diff_met[index_met]]
		#df_set_met = df_set.loc[df_set['Zini'] == 0.018]

		largo_set_met = len(df_set_met['Zini']) # numero de puntos con metalicidad especifica

		largo_oc = len(g) # numero de estrellas del cluster
		largo_set = len(df_set_met) # set de isocronas, con distintas edades, e igual metalicidad diff_met[index_met]

		#sys.exit(0)

		# ALERTAAAAAAAAAAAAAAAAAAAAAAAAAAA!!!!!
		# DETECTAR CUANTOS VALORES DISTINTOS DE LOG AGE HAY.
		#SELECCIONAR PUNTOS CON EDAD ESPECIFICA
		#ESA ES UNA ISOCHRONA, SIN SHIFT EN COLOR EXCESS Y DIST
		diff_ages = df_set_met.log_age_yr.unique()
		diff_ages = np.sort(diff_ages)

		age_range_min = diff_ages[0]
		age_range_max = diff_ages[-1]

		len_diff_ages = len(diff_ages)

		if verbose == False: # print at least some information
				print '\n'
				print '------------------------------------------------------------------------------------------------------------------------------'

		print '\t Cluster: ', cluster_name, '\tZ_range: [', initial_Z, ', ', final_Z, ']\t logA_range: [', initial_age, ', ', final_age, ']'
		print '\t Metallicity: ', df_set_met['Zini'].iloc[0] , "\t edades: ", diff_ages, "\t edades totales: ", len_diff_ages


		index_age = 0
		while index_age < len_diff_ages: # con esto se recorre el set completo de isocronas con metalicidad constante

			df_set_age = df_set_met.loc[df_set_met['log_age_yr'] == diff_ages[index_age]]

			largo_set_age = len(df_set_age['Zini']) # numero de puntos con edad especifica

			if verbose == True:
				print 'largo_oc, largo_set, largo_set_age: \t', largo_oc, largo_set, largo_set_age, '\n'
			if verbose == False: # print at least some information
				#if index_age == 0:
				#	print '\n'
				#	print '------------------------------------------------------------------------------------------------------------------------------'
				print '\tZ = %s \tlogAge = %s'%(df_set_age['Zini'].iloc[0], df_set_age['log_age_yr'].iloc[0]) 

			#i = 0
			#while i < largo_set_age: # con esto se recorre el set de isocronas con metalicidad y edad constante
		
			#	print i, ' \ ' , largo_set_age 

			#if len(index_ofthiscluster) > 1:
			#	EB_V_guess = df_highprob['E(B-V)'].iloc[index_ofthiscluster[0]]  # if the cluster is in the catalog more than once
			#else:
			#	EB_V_guess = df_highprob['E(B-V)'].iloc[index_ofthiscluster] # unique index # NOW THIS IS ABOVE!, OUTSIDE THE LOOPS
			#print EB_V_guess
		
			time_before_ebv = timeit.default_timer()

			# crear grid de E(B-V)
			# linea creadora, centrada en valores conocidos
			gridsize_ebv = 0.1
			EB_V_array = np.arange( EB_V_guess - 0.5, EB_V_guess + 0.5, gridsize_ebv) # start, end, step

			EB_V_range_min = EB_V_guess - 0.5
			EB_V_range_max = EB_V_guess + 0.5

			largo_grid_ebv = len(EB_V_array)

			# recorrerlo con j
			j = 0
			while j < largo_grid_ebv: # con esto se recorre el set de isocronas con metalicidad, edad y color excess constante
				
				gridsize_dist = 10.
			
				if dist_guess == 0:
					dist_array = np.arange(50, 400., gridsize_dist )   # we will use later log10(dist), so dist has to be != 0
					dist_range_min = 50
					dist_range_max = 500

				elif dist_guess - 300. <= 0:
					dist_array = np.arange(50, dist_guess + 200., gridsize_dist )
					dist_range_min = 50
					dist_range_max = dist_guess + 200.
				else:
					dist_array = np.arange(dist_guess - 200., dist_guess + 200., gridsize_dist )
					dist_range_min = dist_guess - 200.
					dist_range_max = dist_guess + 200.

				if verbose == True:
					print cluster_name, '\tZ_range: [%.4f, %.4f] \t logAge_range: [%.2f, %.2f] \t E(B-V)_range: [%.2f, %.2f] \t dist_range: [%i, %i]'%( met_range_min, met_range_max, age_range_min, age_range_max, EB_V_range_min, EB_V_range_max, dist_range_min, dist_range_max) 
					print '\t\t\t    Z    log_age E(B-V)  dist '

				#   	A_G/A_V = 0.85926 A_BP/A_V = 1.06794   A_RP/A_V = 0.65199
				# These values are for a G2V star, using Cardelli_etal89 + ODonnell94 extinction curve with Rv=3.1.
				df['g_ebv_corrected'] = df['phot_g_mean_mag'] - 0.85926*3.1*EB_V_array[j]
				df['bp_ebv_corrected'] = df['phot_bp_mean_mag'] - 1.06794*3.1*EB_V_array[j]
				df['rp_ebv_corrected'] = df['phot_rp_mean_mag'] - 0.65199*3.1*EB_V_array[j]

				largo_grid_dist = len(dist_array)

				time_before_dists = timeit.default_timer()

				# recorrerlo con k
				k = 0
				while k < largo_grid_dist: # con esto se recorre el set de isocronas con metalicidad, edad, color excess y distancia constante

					df['g_abs'] = df['g_ebv_corrected'] + 5. - 5.*np.log10(dist_array[k])
					df['bp_abs'] = df['bp_ebv_corrected'] + 5. - 5.*np.log10(dist_array[k])
					df['rp_abs'] = df['rp_ebv_corrected'] + 5. - 5.*np.log10(dist_array[k])

					g_abs = df['g_abs']
					bp_abs = df['bp_abs']
					rp_abs = df['rp_abs']

					df_set_age = df_set_age_0.loc[df_set_age_0['G'] < (max(g_abs) + 2) ].copy()
					# The labels are: 0=PMS, 1=MS, 2=SGB, 3=RGB, (4,5,6)=different stages of CHEB, 7=EAGB, 8=TPAGB. 
					#df_set_age = df_set_age_0.loc[df_set_age_0['label'] <= 1 ].copy()
					df_set_age = df_set_age.reset_index(drop=True)


					# mag err ^2 = (1.086 err flux/ flux) ^2 + sigma_zero point ^2 
					# assuming sigma_zero_point = 0:
					# or 1.09*flux/flux_err from http://dc.zah.uni-heidelberg.de/tableinfo/gaia.dr2light?tapinfo=True#note-ruwe

					df['e_g']  = 1.09*df['phot_g_mean_flux_error']/df['phot_g_mean_flux']
					df['e_bp'] = 1.09*df['phot_bp_mean_flux_error']/df['phot_bp_mean_flux']
					df['e_rp'] = 1.09*df['phot_rp_mean_flux_error']/df['phot_rp_mean_flux']

					#e_g_abs  = 0.03
					#e_bp_abs = 0.03
					#e_rp_abs = 0.03
					df['e_g_abs']  = df['e_g']
					df['e_bp_abs'] = df['e_bp']
					df['e_rp_abs'] = df['e_rp']

					e_g_abs = df['e_g_abs']
					e_bp_abs = df['e_bp_abs']
					e_rp_abs = df['e_rp_abs']

					sum_red_chi2 = 0

					contador_nans = 0

					#print k, '/', largo_grid_dist

					index_oc = 0
					while index_oc < largo_oc:

						#point = [(df['bp_abs'].iloc[index_oc], df['rp_abs'].iloc[index_oc])]
						point = [(df['bp_abs'].iloc[index_oc] - df['rp_abs'].iloc[index_oc], df['g_abs'].iloc[index_oc])]
					
						G = df_set_age['G'].values
						G_BP = df_set_age['G_BP'].values
						G_RP = df_set_age['G_RP'].values
						#points = tuple(zip(G_BP, G_RP))
						points = tuple(zip(G_BP-G_RP, G))
				
						dimensional_dist_array = distance.cdist(point, points, 'euclidean')
						index_min_dist = dimensional_dist_array[0].argmin()

						closest_g_val = df_set_age['G'].iloc[index_min_dist]
						closest_bp_val = df_set_age['G_BP'].iloc[index_min_dist]
						closest_rp_val = df_set_age['G_RP'].iloc[index_min_dist]
						closest_bp_rp_val = closest_bp_val - closest_rp_val

						#print index_oc
						#print 'antes'
						# just to record what is the closest point in the isochrone
						#df['closest_G'].iloc[index_oc] = closest_g_val
						#df['closest_G_BP'].iloc[index_oc] = closest_bp_val
						#df['closest_G_RP'].iloc[index_oc] = closest_rp_val
						#df.closest_G.iloc[index_oc] = closest_g_val
						#df.closest_G_BP.iloc[index_oc] = closest_bp_val  
						#df.closest_G_RP.iloc[index_oc] = closest_rp_val
						#print 'despues'

						# compute chi2

						#red_chi2 = ( ( g_abs[index_oc] - df_set_age['G'].iloc[index_min_dist] )/e_g_abs[index_oc]**2 + ( bp_abs[index_oc] - df_set_age['G_BP'].iloc[index_min_dist] )/e_bp_abs[index_oc]**2 ) / 2. 
						#red_chi2 = ( ( ( g_abs[index_oc] - df_set_age['G'].iloc[index_min_dist] )/e_g_abs )**2. + ( ( bp_abs[index_oc] - df_set_age['G_BP'].iloc[index_min_dist] )/e_bp_abs )**2. ) / 2. 
						
						#red_chi2 = ( ( ( bp_abs[index_oc] - closest_bp_val )/e_bp_abs[index_oc] )**2. + ( ( rp_abs[index_oc] - closest_rp_val )/e_rp_abs[index_oc] )**2. ) / 2. 
						red_chi2 = ( ( ( (bp_abs[index_oc] - rp_abs[index_oc]) - closest_bp_rp_val )/sqrt(e_bp_abs[index_oc]**2 + e_rp_abs[index_oc]**2) )**2. + ( ( g_abs[index_oc] - closest_g_val )/e_g_abs[index_oc] )**2. ) / 2. 

						if math.isnan(red_chi2):
							red_chi2 = 0 # no sumara a sum_red_chi2
							contador_nans = contador_nans+1

						sum_red_chi2 = sum_red_chi2 + red_chi2

						index_oc = index_oc+1

					# compute the overall/average chi2 for the OC stars w/respect to this isochrone
					avg_redchi2 = sum_red_chi2 / (largo_oc - contador_nans) 

					if avg_redchi2 < min_redchi2:
						min_redchi2 = avg_redchi2

						met_min = df_set_age['Zini'].iloc[0]
						age_min = df_set_age['log_age_yr'].iloc[0]
						ebv_min = EB_V_array[j]
						dist_min = dist_array[k]

						avg_redchi2_min = avg_redchi2

					if verbose == True:
						print '\t\t\t %.4f   %.2f   %.2f   %i  '%(df_set['Zini'].iloc[0], df_set_age['log_age_yr'].iloc[0], EB_V_array[j], dist_array[k]), '\t avg_chi2: ', avg_redchi2

					# go check the next dist
					k = k+1	


				if verbose == True:
					print '------------------------------------------------------------------------------------------------------------------------------------------------------------'

				# ---------------------------------------------------------------------------------- #

				time_after_dists = timeit.default_timer()
				total_time_dists = time_after_dists - time_before_dists
				mins_time_dists, secs_time_dists = divmod(total_time_dists, 60)
				hours_time_dists, mins_time_dists = divmod(mins_time_dists, 60)

				if verbose == True:
					sys.stdout.write("Time spent on this grid: %d:%d:%.2d\n\n" % (hours_time_dists, mins_time_dists, secs_time_dists))
				# ---------------------------------------------------------------------------------- #

				# go check the next ebv
				j = j+1


			time_after_ebv = timeit.default_timer()
			total_time_ebv = time_after_ebv - time_before_ebv
			mins_time_ebv, secs_time_ebv = divmod(total_time_ebv, 60)
			hours_time_ebv, mins_time_ebv = divmod(mins_time_ebv, 60)

			if verbose == False:
				sys.stdout.write("\t\t\t\t\t(Time spent on this grid: %d:%d:%.2d)\n" % (hours_time_ebv, mins_time_ebv, secs_time_ebv))

			# go check the next age
			index_age = index_age + 1

		#break # delete the line df_set_met = df_set.loc[df_set['Zini'] == 0.018] and uncomment the one before

		# go check next metallicity
		index_met = index_met + 1

	current_initial_age = current_final_age
	current_final_age = current_initial_age + 10

#print '\n', cluster_name, '\n--------------- MINIMOS --------------- \t LITERATURE\nZ = %.4f ([Fe/H] = %.4f) \t\t Z_lit = %.4f ([Fe/H] = %4.f) \nlogAge = %.2f \t\t\t\t logAge_lit = %.2f (Dias), %.2f (Kharchenko) \nE(B-V) = %.2f \t\t\t\t E(B-V)_lit = %.2f \ndist = %i \t\t\t\t dist_lit = %i \n\n avg_redchi2 = %f \n---------------------------------------'%(met_min, np.log10(met_min/Zo), 10.**df_highprob['[Fe/H]'].iloc[index_ofthiscluster]*Zo, df_highprob['[Fe/H]'].iloc[index_ofthiscluster], age_min, df_highprob['Age_D02'].iloc[index_ofthiscluster], df_highprob['logt_K13'].iloc[index_ofthiscluster], ebv_min, df_highprob['E(B-V)'].iloc[index_ofthiscluster], dist_min, df_highprob['Dist'].iloc[index_ofthiscluster], avg_redchi2_min)
print '\n', cluster_name, '\n--------------- MINIMOS --------------- \t LITERATURE\nZ = %.4f ([Fe/H] = %.4f) \t\t Z_lit = %.4f ([Fe/H] = %4.f) \nlogAge = %.2f \t\t\t\t logAge_lit = %.2f (Dias), %.2f (Kharchenko) \nE(B-V) = %.2f \t\t\t\t E(B-V)_lit = %.2f \ndist = %i \t\t\t\t dist_lit = %i \n\n avg_redchi2 = %f \n---------------------------------------'%(met_min, np.log10(met_min/Zo), 10.**df_highprob_0['FeH'].iloc[index_ofthiscluster]*Zo, df_highprob_0['FeH'].iloc[index_ofthiscluster], age_min, df_highprob_0['Age_D02'].iloc[index_ofthiscluster], df_highprob_0['logt_K13'].iloc[index_ofthiscluster], ebv_min, df_highprob_0['E(B-V)'].iloc[index_ofthiscluster], dist_min, df_highprob_0['Dist'].iloc[index_ofthiscluster], avg_redchi2_min)

# save table with stars of the OC, adding absolute magnitudes according to the best E(B-V) and dist
df['g_ebv_corrected'] = df['phot_g_mean_mag'] - 0.85926*3.1*ebv_min # ebv_min
df['bp_ebv_corrected'] = df['phot_bp_mean_mag'] - 1.06794*3.1*ebv_min
df['rp_ebv_corrected'] = df['phot_rp_mean_mag'] - 0.65199*3.1*ebv_min
df['g_abs'] = df['g_ebv_corrected'] + 5. - 5.*np.log10(dist_min)
df['bp_abs'] = df['bp_ebv_corrected'] + 5. - 5.*np.log10(dist_min)
df['rp_abs'] = df['rp_ebv_corrected'] + 5. - 5.*np.log10(dist_min)
df.to_csv( catalog_name.replace('.csv','_wAbsMag_bestIso.csv'), index=False )


# recover isochrone data
# into dataframe
# save dataframe in the directory
#bestset_name = './Isochrones/0.0152_isochs_full' 


# note that, for example, 7.95 rounds to 8.0

if round(age_min,1) <= age_min:  #     7.9    7.9    8.0       7.9   7.92    8.0         8.0  8.0 8.1
	print 'age_min = ', age_min, '\tFile where this is found: ',round(age_min,1), round(age_min,1) + 0.1
	bestset_name = './Isochrones/step10/output%s-%s.dat' % ( int( round( round(age_min,1)*100) ) , int( round( (round(age_min,1) + 0.1)*100) ) ) 

else: # if age_min < round(age_min)     7.9   7.97  8.0    
	print 'age_min = ', age_min, '\tFile where this is found: ', round(age_min,1) - 0.1, round(age_min,1)	
	bestset_name = './Isochrones/step10/output%s-%s.dat' % ( int( round( (round(age_min,1) - 0.1)*100 ) ), int( round( round(age_min,1)*100) ) ) 

#bestset_name = './Isochrones/step10/output%s_%s.dat' % ( , ) 
print bestset_name

table_temp = ascii.read(bestset_name)

df_temp = table_temp.to_pandas()
#df_temp.rename(columns={'col1': 'Z', 'col2': 'log_age_yr', 'col3': 'M_ini','col4': 'M_act','col5': 'LogL/Lo','col6': 'LogTe','col7': 'LogG','col8': 'mbol','col9': 'G','col10': 'G_BP','col11': 'G_RP','col12': 'int_IMF','col13': 'stage'}, inplace=True)
df_temp.rename(columns={'col1':'Zini','col2':'log_age_yr','col3':'Mini','col4':'int_IMF','col5':'Mass','col6':'logL','col7':'logTe','col8':'logg','col9':'label','col10':'McoreTP','col11':'C_O','col12':'period0','col13':'period1','col14':'pmode', 'col15':'Mloss','col16':'tau1m','col17':'X','col18':'Y','col19':'Xc','col20':'Xn','col21':'Xo','col22':'Cexcess','col23':'Z','col24':'mbolmag','col25':'G','col26':'G_BP','col27':'G_RP','col28':'B_Tmag','col29':'V_Tmag','col30':'Jmag','col31':'Hmag','col32':'Ksmag'}, inplace=True)
df_save_prev = df_temp.loc[df_temp['log_age_yr'] == age_min]
df_save = df_save_prev.loc[df_save_prev['Zini'] == met_min]
df_save['bestFit_E_B_V'] = ebv_min
df_save['bestFit_dist'] = dist_min

bestisoc_name = './Clusters/%s/%s_bestIsochrone.csv'%(cluster_name,cluster_name)
df_save.to_csv( bestisoc_name, index=False )


del(df)
del(df_set)
del(df_set_age)

del df_set_age_0

del(df_temp)
del(df_save)

del(df_highprob)
del df_highprob_0

# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d\n\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #