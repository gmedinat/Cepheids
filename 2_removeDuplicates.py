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

#catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1.csv'
#atalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_12-04-19.csv'
catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_27-05-19.csv'
df = pd.read_csv('./'+catalog_name)


df_duplicates_removed = df.drop_duplicates(subset=['Cluster', 'RA_deg','DEC_deg','ra','dec'], keep='first')

df_duplicates_removed.to_csv( './'+catalog_name.replace('.csv', '_dupRemoved.csv'), index=False )




