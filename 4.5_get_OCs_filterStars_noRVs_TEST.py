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

import timeit

reload(sys)
sys.setdefaultencoding('utf8')

def votable_to_pandas(votable_file):
	#votable = parse(votable_file)
	votable = votable_file
	table = votable.get_first_table().to_table(use_names_over_ids=True)
	return table.to_pandas()


start = timeit.default_timer()

sigma_num = 5.
signif_err = 1.5

# be careful with the signif_err value.
# if its too big (>3), objects with small proper motions and parallaxes belonging to distant
# clusters (which have small pm's and parallaxes) will be rejected ---> not enough stars for isochrone fitting

# if its too small (~1), distant clusters will have more stars, but closer ones will be more "contamined" by stars with poor measurements

# 5
# 3

catalog_name = 'highProbTest_MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior_wProbs_wNames'
#catalog_name = 'ascc64_combo.csv'
#catalog_name = 'trumpler3_combo.csv'
#catalog_name = 'ngc6229_combo.csv'
catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_wProbs.csv'
catalog_name = 'gaiaV_combo_test'

df_ini = pd.read_csv(catalog_name)

largo_total = len(df_ini)

# posterior mahalanobis dist > cierto nivel
#df = df_ini.loc[df_ini['p_val'] >= 0.05]
df = df_ini.loc[df_ini['P_B_A'] >= 0.05]


#df = pd.read_csv(catalog_name)

largo_df = len(df)

lista_clusters = list()

diff_name = df.Cluster_name.unique()
diff_name = np.sort(diff_name)

largo = len(diff_name)
#largo = len(diff_name)


clusters_number = 0
print '\n'
i = 0
while i < largo:

	#subdf = df.iloc[i]
	subdf = df.loc[df['Cluster_name'] == diff_name[i]]
	subdf = subdf.iloc[0]

	Cluster_Name = subdf['Cluster_name']

	RA_cluster    = subdf['RA_deg']
	DEC_cluster    = subdf['DEC_deg']

	pmRA_cluster  = subdf['pmRA_withCG']
	epmRA_cluster = subdf['e_pmRA_withCG']
	pmDEC_cluster  = subdf['pmDEC_withCG']
	epmDEC_cluster = subdf['e_pmDEC_withCG']

	parallax_cluster = subdf['parallax_cluster_withCG']
	eparallax_cluster = subdf['e_parallax_cluster_withCG']

	RV_cluster = subdf['RV_withCarrera19_CARRERA']
	eRV_cluster = subdf['e_RV_withCarrera19_CARRERA']

	r1 = subdf['r1']  # degrees
	if r1 <= 1.:
		search_rad = 1.5*r1
	elif  1. < r1 and r1 <= 2.:
		search_rad = r1     # for the biggest clusters
	elif  2. < r1 and r1 <= 3.:
		search_rad = r1/1.5     # for the biggest clusters	
	else:
		search_rad = r1/2.5            # there are clusters with sizes up to 6.

	print Cluster_Name, '   %s/%s'%(i+1, largo), '\t', RA_cluster, ', ', DEC_cluster, '\t r1: ', r1, '\t search radius: ', search_rad,' deg'

	if RV_cluster == 0 and eRV_cluster == 0:
		print '\n\tpmra = %.2f'%pmRA_cluster, ' +- %.2f'%epmRA_cluster, '\tpmdec = %.2f'%pmDEC_cluster, ' +- %.2f'%epmDEC_cluster, '\tpm_abs = %.2f'%sqrt(pmRA_cluster**2 + pmDEC_cluster**2), ' +- %.2f'%sqrt(epmRA_cluster**2 + epmDEC_cluster**2), '\tparallax = %.2f'%parallax_cluster, ' +- %.2f'%eparallax_cluster
	else:
		print '\n\tpmra = %.2f'%pmRA_cluster, ' +- %.2f'%epmRA_cluster, '\tpmdec = %.2f'%pmDEC_cluster, ' +- %.2f'%epmDEC_cluster, '\tpm_abs = %.2f'%sqrt(pmRA_cluster**2 + pmDEC_cluster**2), ' +- %.2f'%sqrt(epmRA_cluster**2 + epmDEC_cluster**2), '\tparallax = %.2f'%parallax_cluster, ' +- %.2f'%eparallax_cluster, '\tRV = %.2f'%RV_cluster, ' +- %.2f'%eRV_cluster

	indir = 'Clusters/%s'%Cluster_Name

	file_name = indir+"/%s_stars.vot"%Cluster_Name


	if (not Cluster_Name in str(lista_clusters)):

		lista_clusters.append(Cluster_Name)
		clusters_number = clusters_number+1

		catalogFile_name = file_name.replace('stars.vot', 'stars_noConstraints.vot')
		print '\tfrom ', catalogFile_name, '   to ', file_name.replace('.vot','_noRVs_TEST.csv')
		print '\treading' , catalogFile_name, '...\n'
		votable_file = parse(catalogFile_name)

		df_OC_0 = votable_to_pandas(votable_file)
		df_OC   = df_OC_0.loc[df_OC_0['ruwe'] <= 1.4] 

		A = (pmRA_cluster != 0 and epmRA_cluster != 0)
		B = (pmDEC_cluster != 0 and epmDEC_cluster != 0)
		C = (parallax_cluster != 0 and eparallax_cluster !=0)

		if A == True and B == True and C == True: 
			print '\nA and B and C', '\t\t\t\t\t Lenght without cuts: ', len(df_OC_0)
	
			df_OC = df_OC.loc[ (np.sqrt(df_OC['pmra']**2. + df_OC['pmdec']**2.) / np.sqrt(df_OC['pmra_error']**2. + df_OC['pmra_error']**2.) > signif_err) & (df_OC['parallax'] / df_OC['parallax_error'] > signif_err) ]

			pmAbs_cluster = sqrt(pmRA_cluster**2 + pmDEC_cluster**2)
			epmAbs_cluster = sqrt(epmRA_cluster**2 + epmDEC_cluster**2)

			print '\t Cuts: '
			print '\t %f <= pmra <= %f '%( pmRA_cluster-sigma_num*epmRA_cluster, pmRA_cluster+sigma_num*epmRA_cluster )
			print '\t %f <= pmdec <= %f '%( pmDEC_cluster-sigma_num*epmDEC_cluster, pmDEC_cluster+sigma_num*epmDEC_cluster )
			print '\t %f <= pmAbs <= %f '%( pmAbs_cluster-sigma_num*epmAbs_cluster, pmAbs_cluster+sigma_num*epmAbs_cluster)
			print '\t %f <= parallax <= %f '%( parallax_cluster-sigma_num*eparallax_cluster, parallax_cluster+sigma_num*eparallax_cluster )


			#df1 = df_OC.loc[(df_OC['pmra'] >= pmRA_cluster-sigma_num*epmRA_cluster) & (df_OC['pmra'] <= pmRA_cluster+sigma_num*epmRA_cluster)]	
			#df2 = df1.loc[(df1['pmdec'] >= pmDEC_cluster-sigma_num*epmDEC_cluster) & (df1['pmdec'] <= pmDEC_cluster+sigma_num*epmDEC_cluster)]
			df2 = df_OC.loc[(sqrt(df_OC['pmra']**2 + df_OC['pmdec']**2)  >= pmAbs_cluster-sigma_num*epmAbs_cluster) & (sqrt(df_OC['pmra']**2 + df_OC['pmdec']**2)  <=  pmAbs_cluster+sigma_num*epmAbs_cluster)]
			df3 = df2.loc[(df2['parallax'] >= parallax_cluster-sigma_num*eparallax_cluster) & (df2['parallax'] <= parallax_cluster+sigma_num*eparallax_cluster)]
			df3.to_csv( file_name.replace('.vot','_noRVs_TEST.csv'), index=False )

			print '\t\t\t\t\t\t Lenght with quality cut: ', len(df_OC), '\t Lenght with quality+first cut: ', len(df2), '\t Lenght with all cuts: ',len(df3)

			#del(df1) 
			del(df2)
			del(df3)

		elif A == True and B == True and C == False:
			print 'A and B not C'
			
			df_OC = df_OC.loc[ (np.sqrt(df_OC['pmra']**2. + df_OC['pmdec']**2.) / np.sqrt(df_OC['pmra_error']**2. + df_OC['pmra_error']**2.) > signif_err) ]

			pmAbs_cluster = sqrt(pmRA_cluster**2 + pmDEC_cluster**2)
			epmAbs_cluster = sqrt(epmRA_cluster**2 + epmDEC_cluster**2)

			print '\t Cuts: '
			print '\t %f <= pmra <= %f '%( pmRA_cluster-sigma_num*epmRA_cluster, pmRA_cluster+sigma_num*epmRA_cluster )
			print '\t %f <= pmdec <= %f '%( pmDEC_cluster-sigma_num*epmDEC_cluster, pmDEC_cluster+sigma_num*epmDEC_cluster )
			print '\t %f <= pmAbs <= %f '%( pmAbs_cluster-sigma_num*epmAbs_cluster, pmAbs_cluster+sigma_num*epmAbs_cluster)


			#df1 = df_OC.loc[(df_OC['pmra'] >= pmRA_cluster-sigma_num*epmRA_cluster) & (df_OC['pmra'] <= pmRA_cluster+sigma_num*epmRA_cluster)]	
			#df2 = df1.loc[(df1['pmdec'] >= pmDEC_cluster-sigma_num*epmDEC_cluster) & (df1['pmdec'] <= pmDEC_cluster+sigma_num*epmDEC_cluster)]
			df2 = df_OC.loc[(sqrt(df_OC['pmra']**2 + df_OC['pmdec']**2)  >= pmAbs_cluster-sigma_num*epmAbs_cluster) & (sqrt(df_OC['pmra']**2 + df_OC['pmdec']**2)  <=  pmAbs_cluster+sigma_num*epmAbs_cluster)]

			df2.to_csv( file_name.replace('.vot','_noRVs_TEST.csv'), index=False )

			#del(df1)
			del(df2)

		elif A == True and B == False and C == True:
			print 'A not B and C'

			df_OC = df_OC.loc[ (df_OC['pmra'] / df_OC['pmra_error'] > signif_err) & (df_OC['parallax'] / df_OC['parallax_error'] > signif_err) ]

			print '\t Cuts: '
			print '\t %f <= pmra <= %f '%( pmRA_cluster-sigma_num*epmRA_cluster, pmRA_cluster+sigma_num*epmRA_cluster )
			print '\t %f <= parallax <= %f '%( parallax_cluster-sigma_num*eparallax_cluster, parallax_cluster+sigma_num*eparallax_cluster )
			
			df1 = df_OC.loc[(df_OC['pmra'] >= pmRA_cluster-sigma_num*epmRA_cluster) & (df_OC['pmra'] <= pmRA_cluster+sigma_num*epmRA_cluster)]	
			df2 = df1.loc[(df1['parallax'] >= parallax_cluster-sigma_num*eparallax_cluster) & (df1['parallax'] <= parallax_cluster+sigma_num*eparallax_cluster)]
			df2.to_csv( file_name.replace('.vot','_noRVs_TEST.csv'), index=False )

			del(df1)
			del(df2)
			
		elif A == True and B == False and C == False:
			print 'A not B not C'

			df_OC = df_OC.loc[ (df_OC['pmra'] / df_OC['pmra_error'] > signif_err) ]

			print '\t Cuts: '
			print '\t %f <= pmra <= %f '%( pmRA_cluster-sigma_num*epmRA_cluster, pmRA_cluster+sigma_num*epmRA_cluster )
			
			df1 = df_OC.loc[(df_OC['pmra'] >= pmRA_cluster-sigma_num*epmRA_cluster) & (df_OC['pmra'] <= pmRA_cluster+sigma_num*epmRA_cluster)]	
			df1.to_csv( file_name.replace('.vot','_noRVs_TEST.csv'), index=False )

			del(df1)
			
		elif A == False and B == True and C == True:
			print 'not A and B and C'

			df_OC = df_OC.loc[ (df_OC['pmdec'] / df_OC['pmdec_error'] > signif_err) & (df_OC['parallax'] / df_OC['parallax_error'] > signif_err) ]			
			
			print '\t Cuts: '
			print '\t %f <= pmdec <= %f '%( pmDEC_cluster-sigma_num*epmDEC_cluster, pmDEC_cluster+sigma_num*epmDEC_cluster )
			print '\t %f <= parallax <= %f '%( parallax_cluster-sigma_num*eparallax_cluster, parallax_cluster+sigma_num*eparallax_cluster )

			df1 = df_OC.loc[(df_OC['pmdec'] >= pmDEC_cluster-sigma_num*epmDEC_cluster) & (df_OC['pmdec'] <= pmDEC_cluster+sigma_num*epmDEC_cluster)]
			df2 = df1.loc[(df1['parallax'] >= parallax_cluster-sigma_num*eparallax_cluster) & (df1['parallax'] <= parallax_cluster+sigma_num*eparallax_cluster)]
			df2.to_csv( file_name.replace('.vot','_noRVs_TEST.csv'), index=False )

			del(df1)
			del(df2)

		elif A == False and B == True and C == False:
			print 'not A and B not C'

			df_OC = df_OC.loc[ (df_OC['pmdec'] / df_OC['pmdec_error'] > signif_err)  ]

			print '\t Cuts: '
			print '\t %f <= pmdec <= %f '%( pmDEC_cluster-sigma_num*epmDEC_cluster, pmDEC_cluster+sigma_num*epmDEC_cluster )
			
			df1 = df_OC.loc[(df_OC['pmdec'] >= pmDEC_cluster-sigma_num*epmDEC_cluster) & (df_OC['pmdec'] <= pmDEC_cluster+sigma_num*epmDEC_cluster)]
			df1.to_csv( file_name.replace('.vot','_noRVs_TEST.csv'), index=False )

			del(df1)

		elif A == False and B == False and C == True:
			print 'not A not B and C'

			df_OC = df_OC.loc[ (df_OC['parallax'] / df_OC['parallax_error'] > signif_err) ]

			print '\t Cuts: '
			print '\t %f <= parallax <= %f '%( parallax_cluster-sigma_num*eparallax_cluster, parallax_cluster+sigma_num*eparallax_cluster )
				
			df1 = df_OC.loc[(df_OC['parallax'] >= parallax_cluster-sigma_num*eparallax_cluster) & (df_OC['parallax'] <= parallax_cluster+sigma_num*eparallax_cluster)]
			df1.to_csv( file_name.replace('.vot','_noRVs_TEST.csv'), index=False )

			del(df1)	
		
		else: # everything == 0
			print 'not A not B not C'

			print 'No Cuts'
			
			df1 = df_OC			
			df1.to_csv( file_name.replace('.vot','_noRVs_TEST.csv'), index=False )

			del(df1)
		
		del(df_OC)

	#sys.exit(0)	
	print '--------------------------------------------------------------------------------------------------------\n'
	
	i = i+1


print '\nNumber of different clusters = ', clusters_number
# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #
	

