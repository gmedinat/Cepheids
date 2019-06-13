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
signif_err = 3.5

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
catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19_dupRemoved_wPrior_RVMelnikOnly_ErrsTweak_covMatTweak_wProbs.csv'
catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_27-05-19_dupRem_wPrior_RVMelnikOnly_ErrsTweak_covMatTweak_wProbs.csv'

df_ini = pd.read_csv(catalog_name)

largo_total = len(df_ini)

# posterior mahalanobis dist > cierto nivel
#df = df_ini.loc[df_ini['p_val'] >= 0.05]
#df = df_ini.loc[df_ini['P_B_A'] >= 0.05]
df = df_ini.loc[df_ini['P_A*P_B_A'] >= 0.0000001]


#df = pd.read_csv(catalog_name)

largo_df = len(df)

catalog_CG18 = './papers/cantat-gaudins/CG-2018_1229ocs_memberships.csv'
df_cg        = pd.read_csv(catalog_CG18)

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


	subdf_cg = df_cg.loc[df_cg['Cluster'] == Cluster_Name]

	largo_df_cg = len(subdf_cg)

	if largo_df_cg > 0:
		print i+1, '/', largo, 'Cluster: ', Cluster_Name, '\t\t Cluster in the Cantat-Gaudin et al. 2018 membership database!'
		print 'Number of stars: ', largo_df_cg, '\n'

		name_file = './Clusters/%s/%s_starsFromCG18.csv'%(Cluster_Name, Cluster_Name)

		print 'Saving file ', name_file

		subdf_cg.to_csv( name_file, index=False )

	else:
		print i+1, '/', largo, 'Cluster: ', Cluster_Name, '\t\t No match in the Cantat-Gaudin et al. 2018 membership database!\n'
			
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
	
