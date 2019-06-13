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

inputcatalog = 'match_MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_wProbs_WITH_CG18.csv'
inputcatalog = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_RVMelnikOnly_ErrsTweak_covMatTweak_wProbs.csv'
inputcatalog = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_27-05-19_dupRem_wPrior_RVMelnikOnly_ErrsTweak_covMatTweak_wProbs.csv'


df = pd.read_csv(inputcatalog)
#df_highprob = df.loc[df['p_val'] >= 0.05]
#df_highprob = df.loc[df['P_B_A'] >= 0.05]
df_highprob = df.loc[df['P_A*P_B_A'] >= 0.0000001]


diff_ocs = df_highprob.Cluster_name.unique()
diff_ocs = np.sort(diff_ocs)

#sys.argv[0]

len_diff_ocs = len(diff_ocs)

index_ocs = 0
#index_ocs = 39
#index_ocs = 22
while index_ocs < len_diff_ocs: 

	cluster_name = diff_ocs[index_ocs]

	print 'Current Cluster: ', cluster_name, '\t\t\t\t\t', index_ocs+1,'/',len_diff_ocs

	if not os.path.exists( './Clusters/%s/%s_bestIsochrone_CG18.csv'%(cluster_name,cluster_name) ) :

		semistart = timeit.default_timer()

		os.system( 'python 5_isochrone_fit_CG18data.py %s'%(cluster_name) )

		# ---------------------------------------------------------------------------------- #

		semistop = timeit.default_timer()
		semi_time = semistop - semistart

		# output running time in a nice format.
		semimins, semisecs = divmod(semi_time, 60)
		semihours, semimins = divmod(semimins, 60)

		sys.stdout.write("\n\nPartial running time: %d:%d:%d\n\n" % (semihours, semimins, semisecs))

		# ---------------------------------------------------------------------------------- #

		print '--------------------------------------------------------------------------------------------------------------------------------------------------'
		print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
		print '--------------------------------------------------------------------------------------------------------------------------------------------------'
		print '--------------------------------------------------------------------------------------------------------------------------------------------------'
		print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
		print '--------------------------------------------------------------------------------------------------------------------------------------------------\n\n'

	else:
		print '\t\t\tBest isochrone already computed... We are not doing it again. FOR SURE!!\t'

	index_ocs = index_ocs+1

# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d\n\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #