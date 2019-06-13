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
from pylab import *
from numpy import nan
import random
import subprocess
import timeit
import pandas as pd

start = timeit.default_timer()

indir = './Clusters'

#oc_list = glob.glob('./Clusters')

catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_RVMelnikOnly_ErrsTweak_covMatTweak_wProbs.csv'
catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_27-05-19_dupRem_wPrior_RVMelnikOnly_ErrsTweak_covMatTweak_wProbs.csv'


df_ini = pd.read_csv(catalog_name)
df = df_ini.loc[df_ini['P_A*P_B_A'] >= 0.0000001]

#oc_list = os.listdir(indir)
oc_list = df.Cluster_name.unique()
oc_list = np.sort(oc_list)
#print oc_list
print 'Number of clusters: ', len(oc_list), '\n'

i = 1
for oc in oc_list :

	print oc , '\t\t', i, '/', len(oc_list)

	# Westerlund_2_starsFromCG18.csv 
	fileName = indir+'/'+oc+'/'+oc+'_starsFromCG18.csv'

	# Westerlund_2_stars_noConstraints.vot
	fileGAIA = indir+'/'+oc+'/'+oc+'_stars_noConstraints.vot'

	# Westerlund_2_stars_noRVs.csv
	#fileGAIA = indir+'/'+oc+'/'+oc+'_stars_noRVs.csv'

	if os.path.exists(fileName):
		print '\t\t Exists!'

		outFile = fileName.replace('.csv','_matchGaiaDR2.csv')
		print '\t\t', fileGAIA
		print '\t\t', fileName

		print '\t\tMatching catalogs ...'
		#os.system('java -jar stilts.jar tcopy ifmt=fits in=%s/%s ofmt=ascii  out=%s/%s'%(indir, f, outdir, f.replace('.fits','.dat')) )
		
		os.system('java -jar stilts.jar tmatch2 matcher=sky ifmt1=csv in1=%s ifmt2=votable in2=%s values1="_RAJ2000 _DEJ2000" values2="ra dec" join=1and2 find=best1 params="3" fixcols=dups suffix1="_1" suffix2="" ofmt=csv out=%s'%(fileName,fileGAIA,outFile) )
		#os.system('java -jar stilts.jar tmatch2 matcher=sky ifmt1=csv in1=%s ifmt2=csv in2=%s values1="_RAJ2000 _DEJ2000" values2="ra dec" join=1and2 find=best1 params="3" fixcols=dups suffix1="" suffix2="_2" ofmt=csv out=%s'%(fileName,fileGAIA,outFile) )


		print '\t\t\t\t\t\t\t\t\t%s created!\n\n'%outFile

	else:	
		print '\t\t%s doesnt exist!\n\n'%fileName

	i = i+1		




# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #