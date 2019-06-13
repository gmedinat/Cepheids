#!/usr/bin/python3.6

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
import pyvo as vo
from astropy.io.votable import writeto

import timeit

start = timeit.default_timer()


#example of successful query:
#wget -O myResultFile.vot --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT TOP 10 * FROM gaiadr2.gaia_source" 'http://gaia.ari.uni-heidelberg.de/tap/sync'


#SELECT TOP 5
#  source_id, tgas.ra, tgas.dec, tm.raj2000,
#tm.dej2000, hmag, e_hmag
#FROM gaiadr1.tgas_source as tgas
#JOIN twomass AS tm
#ON 1=CONTAINS (
#     POINT('ICRS', tm.raj2000, tm.dej2000),
#     CIRCLE('ICRS', tgas.ra, tgas.dec, 1.5/3600))

#SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe
#FROM gaiadr2.gaia_source as gaia
#WHERE 1=CONTAINS (
#     POINT('ICRS', 272.3805, -20.79),
#     CIRCLE('ICRS', gaia.ra, gaia.dec, 1.5/3600))

#SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe
#FROM gaiadr2.gaia_source as gaia
#WHERE 1=CONTAINS(POINT('ICRS', gaia.ra, gaia.dec),
#CIRCLE('ICRS', 272.3805, -20.79, 60./3600. ))

mode = 'wo_errors'
#mode = 'w_errors'
#sigma_num = 5.

#catalog_name = 'highPostTest_MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior_wProbs'
#catalog_name = 'highProbTest_MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior_wProbs_wNames'
catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_wProbs.csv'
catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_RVMelnikOnly_ErrsTweak_covMatTweak_wProbs.csv'
catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_27-05-19_dupRem_wPrior_RVMelnikOnly_ErrsTweak_covMatTweak_wProbs.csv'

tap_service_url = "http://gaia.ari.uni-heidelberg.de/tap"

df_ini = pd.read_csv(catalog_name)

lista_clusters = list()

largo = len(df_ini)

# posterior mahalanobis dist > cierto nivel
#df = df_ini.loc[df_ini['p_val'] >= 0.05]
#df = df_ini.loc[df_ini['P_B_A'] >= 0.05]
df = df_ini.loc[df_ini['P_A*P_B_A'] >= 0.0000001]


len_df = len(df)

diff_name = df.Cluster_name.unique()
diff_name = np.sort(diff_name)

len_diff_name = len(diff_name)


#diff_met = df_set.Zini.unique()
#diff_met = np.sort(diff_met)
#df_set_met = df_set.loc[df_set['Zini'] == diff_met[index_met]]
#print(df)

clusters_number = 0
print('\nlargo total: ', largo, '\t largo high prob: ', len_df, '\t number of different clusters: ', len_diff_name, '\n\n' )
del df_ini

i = 0
#i = 79
#i = 105
#while i < largo:
while i < len_diff_name:


	#if i == 134:
	#	sys.exit(0)
	#if (i+1)%10==0:
	#	print('\n\t', i+1, '/', largo

	#subdf = df.iloc[i]
	subdf = df.loc[df['Cluster_name'] == diff_name[i]]
	#if len(subdf) > 1:
	#print('len(subdf) > 1')
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
	#print(r1)
	#print(r1.iloc[0])
	#sys.exit(0)
	if r1 <= 1.:
		search_rad = 1.5*r1
	elif  1. < r1 and r1 <= 2.:
		search_rad = r1     # for the biggest clusters
	elif  2. < r1 and r1 <= 3.:
		search_rad = r1/1.5     # for the biggest clusters	
	else:
		search_rad = r1/2.5            # there are clusters with sizes up to 6.

	print(Cluster_Name, '   %s/%s'%(i+1, len_diff_name), '\t', RA_cluster, ', ', DEC_cluster, '\tr1: ', r1, '\t search radius: ', search_rad,' deg')

	#if mode == 'w_errors':
	#	indir = 'Clusters/w_constrains/%s'%Cluster_Name
	#if mode == 'wo_errors':
	#	indir = 'Clusters/no_constrains/%s'%Cluster_Name
	indir = 'Clusters/%s'%Cluster_Name
	if os.path.exists(indir) == False:
		os.system('mkdir '+indir)

	if mode == 'w_errors':
		file_name = indir+"/%s_stars.vot"%Cluster_Name
	if mode == 'wo_errors':
		file_name = indir+"/%s_stars_noConstraints.vot"%Cluster_Name

	icrs_string = "'ICRS'"

	# ADD pmra, pmdec, parallax/distance constrains
	if (not Cluster_Name in str(lista_clusters)) and mode == 'w_errors':
		lista_clusters.append(Cluster_Name)
		clusters_number = clusters_number+1

		A = (pmRA_cluster != 0 and epmRA_cluster != 0)
		B = (pmDEC_cluster != 0 and epmDEC_cluster != 0)
		C = (parallax_cluster != 0 and eparallax_cluster !=0)
		D = (RV_cluster != 0 and eRV_cluster !=0)
		
		if A == True and B == True and C == True and D == True:
			print('A and B and C and D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.parallax < %f  AND gaia.parallax > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+sigma_num*epmRA_cluster, pmRA_cluster-sigma_num*epmRA_cluster, pmDEC_cluster+sigma_num*epmDEC_cluster, pmDEC_cluster-sigma_num*epmDEC_cluster, parallax_cluster+sigma_num*eparallax_cluster, parallax_cluster-sigma_num*eparallax_cluster, RV_cluster+sigma_num*eRV_cluster, RV_cluster-sigma_num*eRV_cluster)

		elif A == True and B == True and C == True and D == False:
			print('A and B and C not D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.parallax < %f  AND gaia.parallax > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+sigma_num*epmRA_cluster, pmRA_cluster-sigma_num*epmRA_cluster, pmDEC_cluster+sigma_num*epmDEC_cluster, pmDEC_cluster-sigma_num*epmDEC_cluster, parallax_cluster+sigma_num*eparallax_cluster, parallax_cluster-sigma_num*eparallax_cluster)

		elif A == True and B == True and C == False and D == True:
			print('A and B not C and D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+sigma_num*epmRA_cluster, pmRA_cluster-sigma_num*epmRA_cluster, pmDEC_cluster+sigma_num*epmDEC_cluster, pmDEC_cluster-sigma_num*epmDEC_cluster, RV_cluster+sigma_num*eRV_cluster, RV_cluster-sigma_num*eRV_cluster)

		elif A == True and B == False and C == True and D == True:
			print('A not B and C and D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.parallax < %f  AND gaia.parallax > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+sigma_num*epmRA_cluster, pmRA_cluster-sigma_num*epmRA_cluster, parallax_cluster+sigma_num*eparallax_cluster, parallax_cluster-sigma_num*eparallax_cluster, RV_cluster+sigma_num*eRV_cluster, RV_cluster-sigma_num*eRV_cluster)
			
		elif A == True and B == True and C == False and D == False:
			print('A and B not C not D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.pmdec < %f  AND gaia.pmdec > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+sigma_num*epmRA_cluster, pmRA_cluster-sigma_num*epmRA_cluster, pmDEC_cluster+sigma_num*epmDEC_cluster, pmDEC_cluster-sigma_num*epmDEC_cluster)
			
		elif A == True and B == False and C == True and D == False:
			print('A not B and C not D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.parallax < %f  AND gaia.parallax > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+sigma_num*epmRA_cluster, pmRA_cluster-sigma_num*epmRA_cluster, parallax_cluster+sigma_num*eparallax_cluster, parallax_cluster-sigma_num*eparallax_cluster)

		elif A == True and B == False and C == False and D == True:
			print('A not B not C and D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+sigma_num*epmRA_cluster, pmRA_cluster-sigma_num*epmRA_cluster, RV_cluster+sigma_num*eRV_cluster, RV_cluster-sigma_num*eRV_cluster)

		elif A == False and B == True and C == True and D == True:
			print('not A and B and C and D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.parallax < %f  AND gaia.parallax > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmDEC_cluster+sigma_num*epmDEC_cluster, pmDEC_cluster-sigma_num*epmDEC_cluster, parallax_cluster+sigma_num*eparallax_cluster, parallax_cluster-sigma_num*eparallax_cluster, RV_cluster+sigma_num*eRV_cluster, RV_cluster-sigma_num*eRV_cluster)
			
		elif A == False and B == True and C == True and D == False:
			print('not A and B and C not D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.parallax < %f  AND gaia.parallax > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmDEC_cluster+sigma_num*epmDEC_cluster, pmDEC_cluster-sigma_num*epmDEC_cluster, parallax_cluster+sigma_num*eparallax_cluster, parallax_cluster-sigma_num*eparallax_cluster)
			
		elif A == False and B == True and C == False and D == True:
			print('not A and B not C and D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmDEC_cluster+sigma_num*epmDEC_cluster, pmDEC_cluster-sigma_num*epmDEC_cluster, RV_cluster+sigma_num*eRV_cluster, RV_cluster-sigma_num*eRV_cluster)
		
		elif A == False and B == False and C == True and D == True:
			print('not A not B and C and D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.parallax < %f  AND gaia.parallax > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, parallax_cluster+sigma_num*eparallax_cluster, parallax_cluster-sigma_num*eparallax_cluster, RV_cluster+sigma_num*eRV_cluster, RV_cluster-sigma_num*eRV_cluster)

		elif A == True and B == False and C == False and D == False:
			print('A not B not C not D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+sigma_num*epmRA_cluster, pmRA_cluster-sigma_num*epmRA_cluster)

		elif A == False and B == True and C == False and D == False:
			print('not A and B not C not D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmdec < %f  AND gaia.pmdec > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmDEC_cluster+sigma_num*epmDEC_cluster, pmDEC_cluster-sigma_num*epmDEC_cluster)

		elif A == False and B == False and C == True and D == False:
			print('not A not B and C not D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.parallax < %f  AND gaia.parallax > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, parallax_cluster+sigma_num*eparallax_cluster, parallax_cluster-sigma_num*eparallax_cluster)

			tap_service = vo.dal.TAPService(tap_service_url)
			tap_result = tap_service.run_async(query)
			writeto(tap_result.table, file_name)

		elif A == False and B == False and C == False and D == True:
			print('not A not B not C and D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, RV_cluster+sigma_num*eRV_cluster, RV_cluster-sigma_num*eRV_cluster)

		else: # everything == 0
			print('not A not B not C not D')
			query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f ))"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad)


		print(query)
		tap_service = vo.dal.TAPService(tap_service_url)
		#tap_result = tap_service.run_async(query)
		tap_result = tap_service.search(query, maxrec=1000000) # 1.000.000 rows allowed
		print('\tWriting ', file_name, '...')
		writeto(tap_result.table, file_name)

		#job = tap_service.submit_job(query)
		#job.execution_duration = 7200
		#job.run()
		#print(job.phase, job.quote, job.url)

	#pmRA_cluster  epmRA_cluster pmDEC_cluster  epmDEC_cluster  parallax_cluster 	eparallax_cluster RV_cluster eRV_cluster

	if (not Cluster_Name in str(lista_clusters)) and mode == 'wo_errors':
		lista_clusters.append(Cluster_Name)
		clusters_number = clusters_number+1

		#query = "SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f ))"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad)
		query = "SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, phot_g_mean_flux, phot_g_mean_flux_error, phot_bp_mean_flux, phot_bp_mean_flux_error, phot_rp_mean_flux, phot_rp_mean_flux_error, phot_proc_mode, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f ))"%(icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad)

		if os.path.exists(file_name) :
			
			print('\tit already exists... No need to re-download.')

		else:
			print(query)
			tap_service = vo.dal.TAPService(tap_service_url)
			#tap_result = tap_service.run_async(query)
			tap_result = tap_service.search(query, maxrec=10000000) # 10.000.000 rows allowed
			print('\tWriting ', file_name, '...')
			writeto(tap_result.table, file_name)

		#job = tap_service.submit_job(query)
		#job.execution_duration = 7200
		#job.run()
		#print(job.phase, job.quote, job.url)

#SELECT source_id, ra, dec, l, b, parallax, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, ra_error, dec_error, parallax_error, pmra_error, pmdec_error, bp_g, bp_rp, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, radial_velocity, radial_velocity_error, ruwe FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT('ICRS', gaia.ra, gaia.dec), CIRCLE('ICRS', 272.3805, -20.79, 60./3600. ))

	print('--------------------------------------------------------------------------------------------------------\n')

	#if not Cluster_Name in str(lista_clusters):
	#	print('repetido!'
	#	lista_clusters.append(Cluster_Name)
	#	print("added ", Cluster_Name

	#	clusters_number = clusters_number+1
		#break
	#print(any(Cluster_Name in s for s in lista_clusters) if 'abc' in str(my_list):
	#lista_clusters.append(Cluster_Name)
	#print(str(lista_clusters)

	i = i+1


print(lista_clusters, '\nNumber of different clusters = ', clusters_number)
# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #
