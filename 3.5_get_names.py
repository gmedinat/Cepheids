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
from numpy.linalg import inv
from scipy.stats import chi2
from scipy import stats
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy import coordinates
import astropy.units as u

import timeit

start = timeit.default_timer()

#result_table = Simbad.query_object("m1")
#print result_table


# works only for ICRS coordinates:
#c = coordinates.SkyCoord("05h35m17.3s -05d23m28s", frame='icrs')
#r = 5 * u.arcminute
#result_table = Simbad.query_region(c, radius=r)


#c = coordinates.SkyCoord(10.625, 41.2, frame='icrs', unit='deg')
#r = 5 * u.arcminute
#r = 50 * u.arcsecond
#result_table = Simbad.query_region(c, radius=r)

#print result_table


catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior'

df = pd.read_csv(catalog_name)
largo = len(df)

i = 0
while i < largo:
    ra, dec = df['ra'].iloc[i], df['dec'].iloc[i]
    print ra, dec,

    c = coordinates.SkyCoord(ra,dec, frame='icrs', unit='deg')
    r = 10 * u.arcsecond
    result = Simbad.query_region(c, radius=r)
    print result

    i = i+1






# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #
