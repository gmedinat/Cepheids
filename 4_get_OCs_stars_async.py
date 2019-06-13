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

#SELECT *
#FROM gaiadr2.gaia_source as gaia
#WHERE 1=CONTAINS (
#     POINT('ICRS', 272.3805, -20.79),
#     CIRCLE('ICRS', gaia.ra, gaia.dec, 1.5/3600))

#SELECT *
#FROM gaiadr2.gaia_source as gaia
#WHERE 1=CONTAINS(POINT('ICRS', gaia.ra, gaia.dec),
#CIRCLE('ICRS', 272.3805, -20.79, 60./3600. ))

#mode = 'wo_errors'
mode = 'w_errors'

#catalog_name = 'highPostTest_MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior_wProbs'
catalog_name = 'highProbTest_MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior_wProbs_wNames'

df = pd.read_csv(catalog_name)

lista_clusters = list()

largo = len(df)

clusters_number = 0
print '\n'
i = 0
while i < largo:

	#if (i+1)%10==0:
	#	print '\n\t', i+1, '/', largo

	subdf = df.iloc[i]

	Cluster_Name = subdf['Cluster_name']

	RA_cluster    = subdf['RA_deg']
	DEC_cluster    = subdf['DEC_deg']

	pmRA_cluster  = subdf['pmRA_withCG']
	epmRA_cluster = subdf['e_pmRA_withCG']
	pmDEC_cluster  = subdf['pmDEC_withCG']
	epmDEC_cluster = subdf['e_pmDEC_withCG']

	parallax_cluster = subdf['parallax_cluster_withCG']
	eparallax_cluster = subdf['e_parallax_cluster_withCG']

	RV_cluster = subdf['RV']
	eRV_cluster = subdf['e_RV']

	r1 = subdf['r1']  # degrees
	search_rad = 1.5*r1

	print Cluster_Name, '\t', RA_cluster, ', ', DEC_cluster, '\t search radius: ', search_rad,' deg'

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

	#os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT TOP 10 * FROM gaiadr2.gaia_source" "http://gaia.ari.uni-heidelberg.de/tap/async" '%file_name)

	icrs_string = "'ICRS'"
	#print "'ICRS'"
	#sys.exit(0)

	# ADD pmra, pmdec, parallax/distance constrains
	if (not Cluster_Name in str(lista_clusters)) and mode == 'w_errors':
        	lista_clusters.append(Cluster_Name)
		clusters_number = clusters_number+1

        # gaia.pmra, gaia.pmdec, gaia.parallax, gaia.radial_velocity

		# this one!
		#os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.parallax < %f  AND gaia.parallax > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster))

		A = (pmRA_cluster != 0 and epmRA_cluster != 0)
		B = (pmDEC_cluster != 0 and epmDEC_cluster != 0)
		C = (parallax_cluster != 0 and eparallax_cluster !=0)
		D = (RV_cluster != 0 and eRV_cluster !=0)

		if A == True and B == True and C == True and D == True:
			print 'A and B and C and D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.parallax < %f  AND gaia.parallax > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.parallax < %f  AND gaia.parallax > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster))



		elif A == True and B == True and C == True and D == False:
			print 'A and B and C not D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.parallax < %f  AND gaia.parallax > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.parallax < %f  AND gaia.parallax > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster))

		elif A == True and B == True and C == False and D == True:
			print 'A and B not C and D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster))


		elif A == True and B == False and C == True and D == True:
			print 'A not B and C and D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.parallax < %f  AND gaia.parallax > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.parallax < %f  AND gaia.parallax > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster))

		elif A == True and B == True and C == False and D == False:
			print 'A and B not C not D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.pmdec < %f  AND gaia.pmdec > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.pmdec < %f  AND gaia.pmdec > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster))

		elif A == True and B == False and C == True and D == False:
			print 'A not B and C not D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.parallax < %f  AND gaia.parallax > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.parallax < %f  AND gaia.parallax > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster))


		elif A == True and B == False and C == False and D == True:
			print 'A not B not C and D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster))


		elif A == False and B == True and C == True and D == True:
			print 'not A and B and C and D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.parallax < %f  AND gaia.parallax > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.parallax < %f  AND gaia.parallax > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster))

		elif A == False and B == True and C == True and D == False:
			print 'not A and B and C not D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.parallax < %f  AND gaia.parallax > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.parallax < %f  AND gaia.parallax > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster))

		elif A == False and B == True and C == False and D == True:
			print 'not A and B not C and D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmdec < %f  AND gaia.pmdec > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster))

		elif A == False and B == False and C == True and D == True:
			print 'not A not B and C and D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.parallax < %f  AND gaia.parallax > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.parallax < %f  AND gaia.parallax > %f AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster))


		elif A == True and B == False and C == False and D == False:
			print 'A not B not C not D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmra < %f  AND gaia.pmra > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmRA_cluster+2*epmRA_cluster, pmRA_cluster-2*epmRA_cluster))

		elif A == False and B == True and C == False and D == False:
			print 'not A and B not C not D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmdec < %f  AND gaia.pmdec > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.pmdec < %f  AND gaia.pmdec > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, pmDEC_cluster+2*epmDEC_cluster, pmDEC_cluster-2*epmDEC_cluster))

		elif A == False and B == False and C == True and D == False:
			print 'not A not B and C not D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.parallax < %f  AND gaia.parallax > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.parallax < %f  AND gaia.parallax > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, parallax_cluster+2*eparallax_cluster, parallax_cluster-2*eparallax_cluster))

		elif A == False and B == False and C == False and D == True:
			print 'not A not B not C and D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f )) AND gaia.radial_velocity < %f  AND gaia.radial_velocity > %f " "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad, RV_cluster+2*eRV_cluster, RV_cluster-2*eRV_cluster))
		else: # everything == 0
			print 'not A not B not C not D'
			print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f ))" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad), '\n'
			os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f ))" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad))

	#pmRA_cluster  epmRA_cluster pmDEC_cluster  epmDEC_cluster  parallax_cluster 	eparallax_cluster RV_cluster eRV_cluster


	if (not Cluster_Name in str(lista_clusters)) and mode == 'wo_errors':
        	lista_clusters.append(Cluster_Name)
		clusters_number = clusters_number+1
		print 'wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f ))" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad), '\n'
		os.system('wget -O %s --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT(%s, gaia.ra, gaia.dec), CIRCLE(%s, %f, %f, %f ))" "http://gaia.ari.uni-heidelberg.de/tap/async" '%(file_name, icrs_string, icrs_string, RA_cluster, DEC_cluster, search_rad))

#SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT('ICRS', gaia.ra, gaia.dec), CIRCLE('ICRS', 272.3805, -20.79, 60./3600. ))

	print '--------------------------------------------------------------------------------------------------------\n'
	
	#if not Cluster_Name in str(lista_clusters):
	#	print 'repetido!'
	#	lista_clusters.append(Cluster_Name)
	#	print "added ", Cluster_Name

	#	clusters_number = clusters_number+1
		#break
	#print any(Cluster_Name in s for s in lista_clusters) if 'abc' in str(my_list):
	#lista_clusters.append(Cluster_Name)    	
	#print str(lista_clusters)


	i = i+1


print lista_clusters, '\nNumber of different clusters = ', clusters_number
# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #
