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

highprob_catalog_name = 'highProbTest_MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior_wProbs_wNames'
df_highprob = pd.read_csv(highprob_catalog_name)

def closest_point(point, points):
    """ Find closest point from a list of points. """
    #return points[cdist([point], points).argmin()]
    return distance.cdist([point], [points], 'euclidean').argmin()

############## INCLUDE MAGNITUDE UNCERTAINTIES!!!!!!!!!!!!!
############# FILTER STARS BY MAGNITUDE, SOMEHOW!!!! 

# first try with 1 OC...
#cluster_name = 'Majaess_225'   # OJO CON LOS CLUSTERS QUE TIENEN DISTANCIAS = 0
#cluster_name = 'NGC_1901'
#cluster_name = 'NGC_6087'
#cluster_name = 'Alessi_24'
#cluster_name = 'Hogg_7'
cluster_name = 'IC_2391'
catalog_name = './Clusters/%s/%s_stars_noRVs.csv'%(cluster_name,cluster_name)
df = pd.read_csv(catalog_name)

index_ofthiscluster = df_highprob[df_highprob['Cluster_name']==cluster_name].index.values.astype(int)[0] # check where in the original highprob catalog is this cluster

g  = df['phot_g_mean_mag']
bp = df['phot_bp_mean_mag']
rp = df['phot_rp_mean_mag']

g_bp =  df['phot_g_mean_mag'] - df['phot_bp_mean_mag']
bp_rp =  df['phot_bp_mean_mag'] - df['phot_rp_mean_mag']

#e_g =
#e_bp = 
#e_ rp =

#e_g_bp =
#e_bp_rp =

# [M/H]=log(Z/Zo), with Zo=0.019   LOG 10 ?
# crear grid de E(B-V)
# linea creadora, centrada en valores conocidos

# only consider metallicities fe/h > -0.6   (typical of OCs?)

# if df_highprob['[Fe/H]'].iloc[index_ofthiscluster] != null   !!!!!

Zo = 0.019
#Z_guess = 10**(df_highprob['[Fe/H]'].iloc[index_ofthiscluster])*Zo # IS IT LOG 10, RIGHT??????????????????!!!!!!!!!!!!!!!!!

# look for the closest value to that one in a list of available metallicities
# select range of metallicities to be checked

# open isochrone set for given metallicity. Each row of this file has a specific age, and a specific mass.
set_name = './Isochrones/0.0152_isochs_full'
table = ascii.read(set_name)
#table = Table.read('testRR1_all_g.fits', format='fits')

Z_range_min = 0.0152
Z_range_max = 0.0152

df_set = table.to_pandas()
#df_set.rename(columns={'col1': 'Z', 'col2': 'log_age_yr', 'col3': 'M_ini','col4': 'M_act','col5': 'LogL/Lo','col6': 'LogTe','col7': 'LogG','col8': 'mbol','col9': 'G','col10': 'G_BP','col11': 'G_RP','col12': 'int_IMF','col13': 'stage'}, inplace=True)


# output new interface: # Zini logAge Mini int_IMF  Mass   logL    logTe  logg  label   McoreTP C_O  period0 period1 pmode  
#Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z 	 mbolmag  Gmag    G_BPmag  G_RPmag  B_Tmag  V_Tmag  Jmag    Hmag    Ksmag
df_set.rename(columns={'col1':'Zini','col2':'logAge','col3':'Mini','col4':'int_IMF','col5':'Mass','col6':'logL','col7':'logTe','col8':'logg','col9':'label','col10':'McoreTP','col11':'C_O','col12':'period0','col13':'period1','col14':'pmode', 'col15':'Mloss','col16':'tau1m','col17':'X','col18':'Y','col19':'Xc','col20':'Xn','col21':'Xo','col22':'Cexcess','col23':'Z','col24':'mbolmag','col25':'G','col26':'G_BP','col27':'G_RP','col28':'B_Tmag','col29':'V_Tmag','col30':'Jmag','col31':'Hmag','col32':'Ksmag'}, inplace=True)
#df_set.rename(columns={'col1': 'Zini', 'col2': 'log_age_yr', 'col3': 'M_ini','col4': 'M_act','col5': 'LogL/Lo','col6': 'LogTe','col7': 'LogG',

age_limit = 8.6989  # 500 Myr
#age_limit = 8.8129  # 650 Myr
df_set = df_set.loc[ df_set['log_age_yr'] < age_limit ] # only consider isochrones with ages < log(500 Myr) = 8.6989   or ages < log(650 Myr) = 8.8129

#df_set = pd.read_csv(set_name)

largo_oc = len(g) # numero de estrellas del cluster
largo_set = len(df_set) # set de isocronas, con distintas edades

#sys.exit(0)

# ALERTAAAAAAAAAAAAAAAAAAAAAAAAAAA!!!!!
# DETECTAR CUANTOS VALORES DISTINTOS DE LOG AGE HAY.
#SELECCIONAR PUNTOS CON EDAD ESPECIFICA
#ESA ES UNA ISOCHRONA, SIN SHIFT EN COLOR EXCESS Y DIST
diff_ages = df_set.log_age_yr.unique()
diff_ages = np.sort(diff_ages)


age_range_min = diff_ages[0]
age_range_max = diff_ages[-1]

len_diff_ages = len(diff_ages)

print '\t Metallicity: ', df_set['Z'].iloc[0] , "\t edades: ", diff_ages, "\t edades totales: ", len_diff_ages

min_redchi2 = 10000

index_age = 0
while index_age < len_diff_ages: # con esto se recorre el set completo de isocronas con metalicidad constante

	df_set_age = df_set.loc[df_set['log_age_yr'] == diff_ages[index_age]]

	largo_set_age = len(df_set_age['Z']) # numero de puntos con edad especifica

	print 'larog_oc, largo_set, largo_set_age: \t', largo_oc, largo_set, largo_set_age, '\n'

	#i = 0
	#while i < largo_set_age: # con esto se recorre el set de isocronas con metalicidad y edad constante
	
	#	print i, ' \ ' , largo_set_age 

	#if len(index_ofthiscluster) > 1:
	#	EB_V_guess = df_highprob['E(B-V)'].iloc[index_ofthiscluster[0]]  # if the cluster is in the catalog more than once
	#else:
	#	EB_V_guess = df_highprob['E(B-V)'].iloc[index_ofthiscluster] # unique index
	EB_V_guess = df_highprob['E(B-V)'].iloc[index_ofthiscluster] # unique index
	#print EB_V_guess
	
	# crear grid de E(B-V)
	# linea creadora, centrada en valores conocidos
	gridsize_ebv = 0.01
	EB_V_array = np.arange( EB_V_guess - 0.5, EB_V_guess + 0.5, gridsize_ebv) # start, end, step

	EB_V_range_min = EB_V_guess - 0.5
	EB_V_range_max = EB_V_guess + 0.5

	largo_grid_ebv = len(EB_V_array)

	# recorrerlo con j
	j = 0
	while j < largo_grid_ebv: # con esto se recorre el set de isocronas con metalicidad, edad y color excess constante

		# crear grid dist
		# linea creadora, centrada en valores conocidos
		dist_guess = df_highprob['Dist'].iloc[index_ofthiscluster]

		#gridsize_dist = 10.
		gridsize_dist = 10.
		
		if dist_guess == 0:
			dist_array = np.arange(10, 400., gridsize_dist )   # we will use later log10(dist), so dist has to be != 0
			dist_range_min = 10
			dist_range_max = 500

		elif dist_guess - 300. <= 0:
			dist_array = np.arange(10, dist_guess + 200., gridsize_dist )
			dist_range_min = 10
			dist_range_max = dist_guess + 200.
		else:
			dist_array = np.arange(dist_guess - 200., dist_guess + 200., gridsize_dist )
			dist_range_min = dist_guess - 200.
			dist_range_max = dist_guess + 200.

		#print '\n',
		#print index_ofthiscluster 
		#print EB_V_guess, EB_V_range_min, EB_V_range_max
		#print dist_guess, dist_range_min, dist_range_max
		print cluster_name, '\tZ_range: [%.4f, %.4f] \t logAge_range: [%.2f, %.2f] \t E(B-V)_range: [%.2f, %.2f] \t dist_range: [%i, %i]'%( Z_range_min, Z_range_max, age_range_min, age_range_max, EB_V_range_min, EB_V_range_max, dist_range_min, dist_range_max) 
		print '\t\t\t    Z    log_age E(B-V)  dist '
		#sys.exit(0)

		# hacer correcion de magnitudes por ebv
		#extinction = # EB_V_array[j]

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

			# aqui tenemos una isochrona con valores especificos de met, edad, ebv y dist

			# hacer correcion de magnitudes corregidas por ebv, por dist

			df['g_abs'] = df['g_ebv_corrected'] + 5. - 5.*np.log10(dist_array[k])
			df['bp_abs'] = df['bp_ebv_corrected'] + 5. - 5.*np.log10(dist_array[k])
			df['rp_abs'] = df['rp_ebv_corrected'] + 5. - 5.*np.log10(dist_array[k])

			g_abs = df['g_abs']
			bp_abs = df['bp_abs']
			rp_abs = df['rp_abs']

			e_g_abs = 0.03
			e_bp_abs = 0.03
			e_rp_abs = 0.03

			sum_red_chi2 = 0

			contador_nans = 0

			index_oc = 0
			while index_oc < largo_oc:

				# for a specific star in the cluster, df[index_oc], look for the closest match in the current isochrone	
				
				#min_dimensional_dist = 10000
				#ii = 0
				#while ii < largo_set_age:

					##dimensional_dist = np.sqrt( ( df['g_abs'].iloc[index_oc] - df_set_age['G'].iloc[ii] )**2. + ( df['bp_abs'].iloc[index_oc] - df_set_age['G_BP'].iloc[ii] )**2. )
					#dimensional_dist = np.sqrt( ( df['bp_abs'].iloc[index_oc] - df_set_age['G_BP'].iloc[ii] )**2. + ( df['rp_abs'].iloc[index_oc] - df_set_age['G_RP'].iloc[ii] )**2. )

					#if dimensional_dist < min_dimensional_dist:
					#	min_dimensional_dist = dimensional_dist
					#	index_min_dist = ii

					#ii = ii+1

				point = [(df['bp_abs'].iloc[index_oc], df['rp_abs'].iloc[index_oc])]
				
				G_BP = df_set_age['G_BP'].values
				G_RP = df_set_age['G_RP'].values
				points = tuple(zip(G_BP, G_RP))
				
				dimensional_dist_array = distance.cdist(point, points, 'euclidean')
				index_min_dist = dimensional_dist_array[0].argmin()

				#print 'Cluster Name:', cluster_name, '\tindex OC star: ', index_oc, 'index isochrone mass: ', index_min_dist
				#print '\t(%f, %f) ---> (%f, %f)\t    best match in isochrone: (%f, %f)'%(g[index_oc], bp[index_oc], g_abs[index_oc], bp_abs[index_oc], df_set_age['G'].iloc[index_min_dist], df_set_age['G_BP'].iloc[index_min_dist] )

				# compute chi2

				#red_chi2 = ( ( g_abs[index_oc] - df_set_age['G'].iloc[index_min_dist] )/e_g_abs[index_oc]**2 + ( bp_abs[index_oc] - df_set_age['G_BP'].iloc[index_min_dist] )/e_bp_abs[index_oc]**2 ) / 2. 
				#red_chi2 = ( ( ( g_abs[index_oc] - df_set_age['G'].iloc[index_min_dist] )/e_g_abs )**2. + ( ( bp_abs[index_oc] - df_set_age['G_BP'].iloc[index_min_dist] )/e_bp_abs )**2. ) / 2. 
				red_chi2 = ( ( ( bp_abs[index_oc] - df_set_age['G_BP'].iloc[index_min_dist] )/e_g_abs )**2. + ( ( rp_abs[index_oc] - df_set_age['G_RP'].iloc[index_min_dist] )/e_bp_abs )**2. ) / 2. 

				if math.isnan(red_chi2):
					red_chi2 = 0 # no sumara a sum_red_chi2
					contador_nans = contador_nans+1
					#print '\t\t ', g_abs[index_oc], df_set_age['G'].iloc[index_min_dist], (g_abs[index_oc] - df_set_age['G'].iloc[index_min_dist]), bp_abs[index_oc], df_set_age['G_BP'].iloc[index_min_dist], ( bp_abs[index_oc] - df_set_age['G_BP'].iloc[index_min_dist] )
					#print index_oc, bp[index_oc], np.log10(dist_array[k]), bp_abs[index_oc]
					#sys.exit(0)	

				sum_red_chi2 = sum_red_chi2 + red_chi2

				index_oc = index_oc+1

			# compute the overall/average chi2 for the OC stars w/respect to this isochrone
			avg_redchi2 = sum_red_chi2 / (largo_oc - contador_nans) 

			if avg_redchi2 < min_redchi2:
				min_redchi2 = avg_redchi2

				met_min = df_set['Z'].iloc[0]
				age_min = df_set_age['log_age_yr'].iloc[0]
				ebv_min = EB_V_array[j]
				dist_min = dist_array[k]

				avg_redchi2_min = avg_redchi2

			#print '\tZ, logAge, E(B-V), dist: %.4f   %.2f   %.2f   %i  '%(df_set['Z'].iloc[0], df_set_age['log_age_yr'].iloc[0], EB_V_array[j], dist_array[k]), '\t avg_chi2: ', avg_redchi2
			print '\t\t\t %.4f   %.2f   %.2f   %i  '%(df_set['Z'].iloc[0], df_set_age['log_age_yr'].iloc[0], EB_V_array[j], dist_array[k]), '\t avg_chi2: ', avg_redchi2

			#sys.exit(0)

			# go check the next dist
			k = k+1	

		print '------------------------------------------------------------------------------------------------------------------------------------------------------------'

		# ---------------------------------------------------------------------------------- #

		time_after_dists = timeit.default_timer()
		total_time_dists = time_after_dists - time_before_dists
		mins_time_dists, secs_time_dists = divmod(total_time_dists, 60)
		hours_time_dists, mins_time_dists = divmod(mins_time_dists, 60)

		sys.stdout.write("Time spent on this grid: %d:%d:%.2d\n\n" % (hours_time_dists, mins_time_dists, secs_time_dists))
		# ---------------------------------------------------------------------------------- #

		# go check the next ebv
		j = j+1

	#i = i+1

	# go check the next age
	index_age = index_age + 1

print '\n', cluster_name, '\n--------------- MINIMOS --------------- \t LITERATURE\nZ = %.4f ([Fe/H] = %.4f) \t\t Z_lit = %.4f ([Fe/H] = %4.f) \nlogAge = %.2f \t\t\t\t logAge_lit = %.2f (Dias), %.2f (Kharchenko) \nE(B-V) = %.2f \t\t\t\t E(B-V)_lit = %.2f \ndist = %i \t\t\t\t dist_lit = %i \n\n avg_redchi2 = %f \n---------------------------------------'%(met_min, np.log10(met_min/Zo), 10.**df_highprob['[Fe/H]'].iloc[index_ofthiscluster]*Zo, df_highprob['[Fe/H]'].iloc[index_ofthiscluster], age_min, df_highprob['Age_D02'].iloc[index_ofthiscluster], df_highprob['logt_K13'].iloc[index_ofthiscluster], ebv_min, df_highprob['E(B-V)'].iloc[index_ofthiscluster], dist_min, df_highprob['Dist'].iloc[index_ofthiscluster], avg_redchi2_min)
#sys.exit(0)


# save table with stars of the OC, adding absolute magnitudes according to the best E(B-V) and dist
df['g_ebv_corrected'] = g # ebv_min
df['bp_ebv_corrected'] = bp 
df['rp_ebv_corrected'] = rp 
df['g_abs'] = df['g_ebv_corrected'] + 5. - 5.*np.log10(dist_min)
df['bp_abs'] = df['bp_ebv_corrected'] + 5. - 5.*np.log10(dist_min)
df['rp_abs'] = df['rp_ebv_corrected'] + 5. - 5.*np.log10(dist_min)
df.to_csv( catalog_name.replace('.csv','_wAbsMag_bestIso.csv'), index=False )


# recover isochrone data
# into dataframe
# save dataframe in the directory
bestset_name = './Isochrones/0.0152_isochs_full'
table_temp = ascii.read(bestset_name)

df_temp = table_temp.to_pandas()
df_temp.rename(columns={'col1': 'Z', 'col2': 'log_age_yr', 'col3': 'M_ini','col4': 'M_act','col5': 'LogL/Lo','col6': 'LogTe','col7': 'LogG','col8': 'mbol','col9': 'G','col10': 'G_BP','col11': 'G_RP','col12': 'int_IMF','col13': 'stage'}, inplace=True)
df_save = df_temp.loc[df_temp['log_age_yr'] == age_min]
df_save['bestFit_E_B_V'] = ebv_min
df_save['bestFit_dist'] = dist_min

bestisoc_name = './Clusters/%s/bestIsochrone.csv'%(cluster_name)
df_save.to_csv( bestisoc_name, index=False )


del(df)
del(df_set)
del(df_set_age)
del(df_temp)
del(df_save)

del(df_highprob)

# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d\n\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #