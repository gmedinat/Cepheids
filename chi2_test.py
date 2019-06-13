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
#from scipy.stats import chi2
from scipy import stats

from scipy.stats import chi2
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)

df = 4
mean, var, skew, kurt = chi2.stats(df, moments='mvsk')

#print chi2.ppf(0.01, df), chi2.ppf(0.9999, df)
x1 = np.linspace(chi2.ppf(0.01, df), chi2.ppf(0.9999, df), 100)
x2 = np.linspace(chi2.ppf(0.01, df, scale=10), chi2.ppf(0.9999, df, scale=10), 100)
x3 = np.linspace(chi2.ppf(0.01, df, scale=10, loc=100), chi2.ppf(0.9999, df, scale=10, loc=100), 100)
#print x1, x2
ax.plot(x1, chi2.pdf(x1, df), 'b-', lw=5, alpha=0.6, label='No param, df=%i'%df)
ax.plot(x2, chi2.pdf(x2, df, scale=10), 'r-', lw=5, alpha=0.6, label='Scale=10, df=%i'%df)
ax.plot(x3, chi2.pdf(x3, df, scale=10, loc=100), 'g-', lw=5, alpha=0.6, label='Scale=10, loc=100, df=%i'%df)

#rv = chi2(df)
#ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')

#vals = chi2.ppf([0.001, 0.5, 0.999], df)
#np.allclose([0.001, 0.5, 0.999], chi2.cdf(vals, df))

#r = chi2.rvs(df, size=1000)

#ax.hist(r, density=True, histtype='stepfilled', alpha=0.2)
ax.legend(loc='best', frameon=False)
plt.xlabel('x', fontsize=12)
plt.ylabel('chi2 pdf', fontsize=12)
plt.tight_layout()
#plt.show()

mean, var, skew, kurt = chi2.stats(df, moments='mvsk')
print mean, var, skew, kurt
print chi2.mean(df, loc=0, scale=1), chi2.mean(df, loc=0, scale=10), chi2.mean(df, loc=10, scale=1), chi2.mean(df, loc=10, scale=10), chi2.mean(3, loc=0, scale=10)


print chi2.sf(3.84, 1), chi2.sf(3.84, 1, scale=1, loc=100)
print chi2.sf(3.84, 4), chi2.sf(3.84, 4, scale=1, loc=100)
print chi2.sf(9.49, 4), chi2.sf(9.49, 4, scale=1, loc=100)







