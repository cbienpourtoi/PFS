"""
File: lightcone_explorer.py
Copyright: Loic Le Tiran, 2014
Contact: loic.le-tiran@obspm.fr
Licence: GNU GPL v3

Description:
Tool for reading the lightcones from 
http://galformod.mpa-garching.mpg.de/qa/mrobs/pages/surveys/PFS.jsp
"""


import numpy as np
import scipy
#import pyfits
import sys
import os
import glob
import csv
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm
from astropy.io import fits
#from string import upper,lower
#from numpy import recfromcsv
#from numpy.random import normal
#import math
from time import gmtime, strftime


plot_extension = ".png"
plot_directory = "./plots/"
def savemyplot(name):
	fig.savefig(plot_directory+name+plot_extension)
	return


"""
Table of redshift bins, properties, limit magnitudes and selection filters
"""

# Make a table of the spectral axis of the Alhambra catalog
z_bins = [2, 3, 4, 5, 6, 7]
mag_limit = [24., 24.3, 24.5, 24.9, 24.9, 25.3]
selec_filter_1 = ['', 'u', 'g', '', '', '']
selec_filter_2 = ['', 'g', 'r', '', '', '']
selec_filter_3 = ['', 'r', 'i', '', '', '']

#filters = Table([filter_names,filter_leff], names=('Filter', 'Lambda eff (A)'), meta={'name': 'table of the filters'})



"""
Open files
"""

file_number = 1

# Lightcones path
conepath = "./data/lightcones/"
conename = "wmap1_bc03_"+str(file_number).zfill(3)+"_igm1.fits"

hdulist = fits.open(conepath+conename)
#cone = hdulist[1].data
cols = hdulist[1].columns
hdulist.close()

sys.exit()

print "There are "+str(len(cone))+" objects in the cone."

"""
Just playing with the files: z distribution.
"""

'''
fig = plt.figure()
plt.title("Redshift distribution for the lightcone")
plt.xlabel("Apparent Redshift (Z_APP)")
plt.ylabel("#")
plt.hist(cone_z2['Z_APP'], bins=200)
plt.show()
savemyplot("z_dist")
plt.close()
'''



"""
Data selection
"""

mask = cone.field('Z_APP')>2

print strftime("%Y-%m-%d %H:%M:%S", gmtime())
# selecting all objects with z>2 takes 20 min
cone_z2 = cone[mask]

print strftime("%Y-%m-%d %H:%M:%S", gmtime())


#TODO : Make field of view to the right size







