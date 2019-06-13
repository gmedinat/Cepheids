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
import pyvo as vo
from astropy.io.votable import writeto


import timeit

start = timeit.default_timer()

#wget -O Clusters/Majaess_225/Majaess_225_stars.vot --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT('ICRS', gaia.ra, gaia.dec), CIRCLE('ICRS', 354.950000, 61.928056, 0.302481 )) AND gaia.pmra < -1.570000  AND gaia.pmra > -2.450000 AND gaia.pmdec < -2.250000  AND gaia.pmdec > -2.930000 AND gaia.radial_velocity < -1.000000  AND gaia.radial_velocity > -113.000000" "http://gaia.ari.uni-heidelberg.de/tap/sync"  

#wget -O Clusters/NGC_1901/NGC_1901_stars.vot --post-data="REQUEST=doQuery&LANG=ADQL&QUERY=SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT('ICRS', gaia.ra, gaia.dec), CIRCLE('ICRS', 79.545833, -68.450000, 0.277500 )) AND gaia.pmra < 1.400000  AND gaia.pmra > 0.160000 AND gaia.pmdec < 13.590000  AND gaia.pmdec > 12.630000 AND gaia.parallax < 2.730000  AND gaia.parallax > 2.290000 AND gaia.radial_velocity < 4.500000  AND gaia.radial_velocity > -3.500000" "http://gaia.ari.uni-heidelberg.de/tap/sync" 

tap_service_url = "http://gaia.ari.uni-heidelberg.de/tap"
#query = "SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT('ICRS', gaia.ra, gaia.dec), CIRCLE('ICRS', 354.950000, 61.928056, 0.302481 )) AND gaia.pmra < -1.570000  AND gaia.pmra > -2.450000 AND gaia.pmdec < -2.250000  AND gaia.pmdec > -2.930000 AND gaia.radial_velocity < -1.000000  AND gaia.radial_velocity > -113.000000" # Majaess_225
#query = "SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT('ICRS', gaia.ra, gaia.dec), CIRCLE('ICRS', 79.545833, -68.450000, 0.277500 )) AND gaia.pmra < 1.400000  AND gaia.pmra > 0.160000 AND gaia.pmdec < 13.590000  AND gaia.pmdec > 12.630000 AND gaia.parallax < 2.730000  AND gaia.parallax > 2.290000 AND gaia.radial_velocity < 4.500000  AND gaia.radial_velocity > -3.500000" # NGC_1901
query = "SELECT * FROM gaiadr2.gaia_source as gaia WHERE 1=CONTAINS(POINT('ICRS', gaia.ra, gaia.dec), CIRCLE('ICRS', 82.854167, 34.245000, 0.120000 )) AND gaia.pmra < -0.030000  AND gaia.pmra > -0.550000 AND gaia.pmdec < -3.120000  AND gaia.pmdec > -3.840000 AND gaia.parallax < 0.608696  AND gaia.parallax > 0.260870" # NGC_1931
tap_service = vo.dal.TAPService(tap_service_url)
tap_result = tap_service.run_async(query)
writeto(tap_result.table, 'filename')
#tap_result = tap_service.search(query)

print(tap_result.votable)

print(tap_service)
print(tap_result)
print(tap_result.table)
#tap_result.write(filename, format='latex')
#tap_result.table.to_table.write('filename', format='votable')

#from astropy.io.votable import parse_single_table
#table = parse_single_table(tap_result.votable).to_table()
#print(table)
#table.write(filename, format='latex')  

writeto(tap_result.table, 'filename', tabledata_format='TABLEDATA')
writeto(tap_result.table, 'filename.fits', tabledata_format='FITS')
#writeto(tap_result.table, 'filename', tabledata_format='fits')
#writeto(tap_result.table, 'filename', tabledata_format='csv')
#writeto(tap_result.table, 'filename', tabledata_format='ascii')
#tap_result.tablewriteto('filename', format='votable')
#tap_result.table.to_table.write('filename', format='votable')

#tap_result.get_first_table().format = "binary"
#tap_result.to_xml("new_votable.xml", tabledata_format="binary")
#tap_result.to_xml("new_votable.fits", tabledata_format="FITS")

#PYVO
			#>>TAP
			#en vez de search , query 
			#tapserviceurl = "http://gaia.ari.uni-heidelberg.de/tap"
			#tap service = pyvo dal TAPservice (tapserviceurl)
			#tap service .run async (query)


			#otra forma:
			#PHASE=RUN
