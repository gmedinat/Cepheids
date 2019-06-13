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

from isochrones.dartmouth import Dartmouth_Isochrone
#from isochrones.dartmouth import Dartmouth_FastIsochrone
from isochrones.mist import MIST_Isochrone
from isochrones.padova import Padova_Isochrone
#from isochrones.basti import Basti_Isochrone

#from isochrones.dartmouth.tri import write_tri
#from isochrones.padova.tri import write_tri
#write_tri()

import timeit

start = timeit.default_timer()

# useful links:
# http://slittlefair.staff.shef.ac.uk/teaching/phy241/resources/plotting_isochrones.html
# https://isochrones.readthedocs.io/en/latest/
# https://isochrones.readthedocs.io/en/latest/_modules/isochrones/isochrone.html

#catalog_name = "./isochronetest_ngc7062.csv"
#df = pd.read_csv(catalog_name)

#g = df['phot_g_mean_mag']
#bp = df['phot_bp_mean_mag']
#rp = df['phot_rp_mean_mag']

#padoviso = Padova_Isochrone(bands=['B','V'])
#padoviso = Padova_Isochrone()
#mist = MIST_Isochrone()
iso = Dartmouth_Isochrone(bands=['B','V'])

#sys.exit(0)

#print(padoviso.bands)

#print(mist.radius(1.0, 9.7, 0.0)) #M/Msun, log10(age), Fe/H
#print(mist.bands)


#mass, age, feh, distance, AV = (0.95, 9.61, -0.2, 200, 0.2)
#print(mist.mag['g'](mass, age, feh, distance, AV))
#print(iso.bands)

#model = iso.isochrone(9.235) # log10 of the age
#print(type(model))

#model_b = model.B_mag
#model_v = model.V_mag
# calculate B-V for this model
#model_bv = model_b - model_v


#fig, axis = plt.subplots(figsize=(8,6))
#axis.plot(bv,v,'r.',alpha=0.1)
#axis.invert_yaxis()
#axis.set_xlabel('B-V')
#axis.set_ylabel('V')

# now I plot the isochrone
#axis.plot(model_bv,model_v)
#plt.show()


#from isochrones import StarModel
#from isochrones.mist import MIST_Isochrone

#spectroscopic properties (value, uncertainty)
#Teff = (5770, 80)
#logg = (4.44, 0.08)
#feh = (0.00, 0.10)
#mist = MIST_Isochrone()
#model  = StarModel(mist, Teff=Teff, logg=logg, feh=feh)
#model.fit()
#model.corner_physical()
#model.save_hdf('starmodel.h5')





# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #
