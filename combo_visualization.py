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

catalog_name = 'highPostTest_MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior_wProbs'
#catalog_name = 'MatchAndWithin_D02+K13_Gaia_2deg_5r1_wPrior_wProbs'

df = pd.read_csv(catalog_name)

largo = len(df)

fig, ax = plt.subplots(figsize=(17,8))
#ax.set_xlim(-0.65, 1.75)
#ax.set_ylim(-0.65, 1.89)
ax.grid(False)

ax.set_xlabel('RA (deg)', fontsize=20)
ax.set_ylabel('DEC (deg)', fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.tight_layout()
#plt.savefig('color_color_grri_PS.png')
#plt.savefig('color_color_grri_PS.eps')

print '\n'
i = 0
while i < largo:

    print i+1, '/', largo

    ra1 = df['RA_deg'].iloc[i]
    dec1 = df['DEC_deg'].iloc[i]

    ra2 = df['ra'].iloc[i]
    dec2 = df['dec'].iloc[i]

    rect_dist = sqrt( (ra1-ra2)**2 + (dec1-dec2)**2   )

    if rect_dist%360 > 180:
        if ra1 < ra2:
            ramin = ra1
            decmin = dec1
            ramax = ra2
            decmax = dec2
        else:
            ramin = ra2
            decmin = dec2
            ramax = ra1
            decmax= dec1
        m = (decmin-decmax)/( (360+ramin) - ramax )
        b = decmax - ramax*m
        dec_inter = m*360 + b

        ax.plot([ramax, 360], [decmax, dec_inter], color='k', linestyle='-', linewidth=2)
        ax.plot([0, ramin], [dec_inter, decmin], color='k', linestyle='-', linewidth=2)

    else:
        ax.plot([ra1, ra2], [dec1, dec2], color='k', linestyle='-', linewidth=2)
    ax.plot(ra1,dec1, 'bo', markersize=10)
    ax.plot(ra2,dec2, 'ro', markersize=3)

    clustername = '   ' + df['Cluster'].iloc[i]
    ax.annotate(clustername, (ra1, dec1))

    print clustername, '\t', rect_dist, '\t', rect_dist%360


    i = i+1

# ---------------------------------------------------------------------------------- #

stop = timeit.default_timer()
total_time = stop - start

# output running time in a nice format.
mins, secs = divmod(total_time, 60)
hours, mins = divmod(mins, 60)

sys.stdout.write("\n\nTotal running time: %d:%d:%d.\n\n" % (hours, mins, secs))

# ---------------------------------------------------------------------------------- #


plt.show()
