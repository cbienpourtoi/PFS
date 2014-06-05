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
from astropy.table import Table
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

dz = 0.5 # Delta in redshift (ie selection between z-dz and z+dz)


# Make a table of the redshifts, selection filters and limit magnitude associated (limit mag always on the redder)
z_bins =         [      2.,       3.,       4.,       5.,       6.,     7.]
mag_limit =      [     24.,     24.3,     24.5,     24.9,     24.9,   25.3]
selec_filter_1 = [      '', 'SDSS_U', 'SDSS_G', 'SDSS_R', 'SDSS_I',    'Z']
selec_filter_2 = [      '', 'SDSS_G', 'SDSS_R', 'SDSS_I', 'SDSS_Z',    'Y']
selec_filter_3 = ['SDSS_R', 'SDSS_R', 'SDSS_I', 'SDSS_Z',      'Y',    'J']
selection_properties = Table([z_bins, mag_limit, selec_filter_1, selec_filter_2, selec_filter_3], names=('z', 'LimitMag', 'Filter1', 'Filter2', 'Filter3'), meta={'name': 'table of the selection properties'})



"""
Open files
"""

file_number = 1

# Lightcones path
conepath = "./data/lightcones/"
conename = "wmap1_bc03_"+str(file_number).zfill(3)+"_igm1.fits"

hdulist = fits.open(conepath+conename)
cone = hdulist[1].data
cols = hdulist[1].columns
hdulist.close()


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
Redshift and limit magnitude selection
"""

#for i in np.arange(len(selection_properties)):
for i in [1]:
	print "redshift: " + str(selection_properties['z'][i])
	mask_z = np.abs( cone.field('Z_APP') - selection_properties['z'][i] ) < dz
	mask_mag = cone.field(selection_properties['Filter3'][i]) < selection_properties['LimitMag'][i]
	mask = mask_mag #& mask_z
	print strftime("%Y-%m-%d %H:%M:%S", gmtime())

# selecting all objects with z>2 takes 20 min
cone_z2 = cone[mask]

print "Number of candidates : " + str(len(cone_z2))

print strftime("%Y-%m-%d %H:%M:%S", gmtime())

fig = plt.figure()
plt.title("Redshift vs Redder Magnitude")
plt.xlabel("Apparent Redshift (Z_APP)")
plt.ylabel("Magnitude in redder color (?)")
plt.hist2d(cone_z2['Z_APP'], cone_z2[selection_properties['Filter3'][i]], bins=1000)
plt.show()
savemyplot("z_vs_mag")
plt.close()


"""
Color selection
"""


fig = plt.figure()
plt.title("Colors")
plt.xlabel(selection_properties['Filter2'][i] + "-" +selection_properties['Filter3'][i])
plt.ylabel(selection_properties['Filter1'][i] + "-" +selection_properties['Filter2'][i])
#plt.hist2d(cone_z2[selection_properties['Filter2'][i]] - cone_z2[selection_properties['Filter3'][i]], cone_z2[selection_properties['Filter1'][i]] - cone_z2[selection_properties['Filter2'][i]], bins=50, range=([-1.,2.5],[-1.,3.5]))
#plt.plot(cone_z2[selection_properties['Filter2'][i]] - cone_z2[selection_properties['Filter3'][i]], cone_z2[selection_properties['Filter1'][i]] - cone_z2[selection_properties['Filter2'][i]], '.')
plt.scatter(cone_z2[selection_properties['Filter2'][i]] - cone_z2[selection_properties['Filter3'][i]], cone_z2[selection_properties['Filter1'][i]] - cone_z2[selection_properties['Filter2'][i]],c=cone_z2['Z_APP'],vmin=1,vmax=5,cmap=plt.cm.spectral)
cbar = plt.colorbar() 
plt.plot([-1., 0.3], [1., 1.], '-b')
plt.plot([1.4, 1.4], [2.6, 3.5], '-b')
plt.plot([0.3, 1.4], [1.0, 2.6], '-b')

plt.xlim(-1.,2.5) 
plt.ylim(-1.,3.5) 
plt.show()
savemyplot("Colors")
plt.close()






#TODO : Make field of view to the right size







