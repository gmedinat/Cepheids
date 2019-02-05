#!/usr/bin/python2.7

import numpy as np
from numpy import nan
from pylab import *
import math
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import csv

from astropy.table import Table
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy
import timeit

import pandas as pd

start = timeit.default_timer()

def ang_distance(ra1, dec1, ra2, dec2):

    return sqrt( (ra1-ra2)**2 + (dec1-dec2)**2 ) 

def within_rad(ra1, dec1, ra2, dec2, ang_rad):

    return 0

# ITS BETTER TO DO THIS STEP IN TOPCAT. It takes too long matching.    

catalog_cep = './Gaiadr2_cepheids_FULL.csv'
ra_cep_name  = 'ra'
dec_cep_name = 'dec'

catalog_ocs = './Dias2002a_OCs.cvs'
ra_oc_name  = 'RA_deg'
dec_oc_name = 'DEC_deg'
angrad_oc_name = 'apparent_radius_deg'

output = './match-Dias-Gaiafull.csv'

df_cep = pd.read_csv(catalog_cep)
df_ocs = pd.read_csv(catalog_ocs)

len_cep = len(df_cep)

len_ocs = len(df_ocs)

print len_cep, len_ocs

first_detection = False

i = 0
#while i < len_ocs:
while i < 10:

	angular_rad = df_ocs.iloc[i][angrad_oc_name]/60.

	#df_ocs['distancia_deg'] = (((df_ocs[ra_oc_name]-df_cep[ra_cep_name])**2.+(df_ocs[dec_oc_name]-df_cep[dec_cep_name])**2.)**0.5)
	# if ...

	print df_ocs.iloc[i]['Cluster'], '\t \t', i+1, '/', len_ocs 

	ra_oc  = df_ocs.iloc[i][ra_oc_name]
	dec_oc = df_ocs.iloc[i][dec_oc_name]

	#df_cep['distancia_deg'] =

	# new_df = FILTRO

	# delete column

	#c = SkyCoord(ra=ra_oc*u.degree, dec=dec_oc*u.degree)  
	#catalog = SkyCoord(ra=df_cep[ra_cep_name]*u.degree, dec=df_cep[dec_cep_name]*u.degree)  
	
	break

	#j = 0
	#while j < len_cep:

		#print j
		#ra_cep  = df_cep.iloc[j][ra_cep_name]
		#dec_cep = df_cep.iloc[j][dec_cep_name]

		#c1 = SkyCoord(ra_oc*u.deg, dec_oc*u.deg)
		#c2 = SkyCoord(ra_cep*u.deg, dec_cep*u.deg)

		#sep = c1.separation(c2)
		#sep_deg = sep.degree

		#if sep_deg <= 1: #or dist <= 1.5*angular_rad:

		#	df_match_0 = pd.concat([df_ocs.iloc[i], df_cep.iloc[j]], axis=1)
			
		#	if first_detection == False:
		#		df_match   = df_match_0		
		#		first_detection = True	

		#	if first_detection == True:
		#		df_match   = pd.concat([df_match, df_match_0], axis=0) 

		#	print 'Combo! \t', df_cep.iloc[j]['designation']
		#	print '\t\t\t (%f, %f) - - (%f, %f)\t\t sep: %f'% (ra_oc,dec_oc,ra_cep,dec_cep, sep_deg) 
			
	#	j = j+1


	print '-----------------------------------------------------------------------------------'

	#break
	i = i+1
   
#print '\n\n len(df_match) = ', len(df_match)
#df_match.to_csv(output)
 
stop = timeit.default_timer()
total_time = stop-start
# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n" % (hours, mins, secs))



