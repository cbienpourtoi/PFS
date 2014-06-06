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

print strftime("%Y-%m-%d %H:%M:%S", gmtime())

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
allcone = hdulist[1].data
cols = hdulist[1].columns
hdulist.close()


print "There are "+str(len(allcone))+" objects in the cone."

"""
Just playing with the files: z distribution.
"""
"""
fig = plt.figure()
plt.title("Redshift distribution for the lightcone")
plt.xlabel("Apparent Redshift (Z_APP)")
plt.ylabel("#")
plt.hist(cone['Z_APP'], bins=200)
plt.show()
savemyplot("z_dist")
plt.close()
"""


"""
Redshift and limit magnitude selection
"""

#for i in np.arange(len(selection_properties)):
for i in [1,2,3,4,5]:
	print "redshift: ~" + str(selection_properties['z'][i]) + ". Filters : " + str(selection_properties['Filter1'][i]) +" "+ str(selection_properties['Filter2'][i]) +" "+ str(selection_properties['Filter3'][i])
	mask_z = np.abs( allcone.field('Z_APP') - selection_properties['z'][i] ) < dz
	mask_mag = allcone.field(selection_properties['Filter3'][i]) < selection_properties['LimitMag'][i]
	mask = mask_mag #& mask_z
	print strftime("%Y-%m-%d %H:%M:%S", gmtime())

	# selecting all objects with z>2 takes 20 min
	cone = allcone[mask]

	print "Number of candidates with mag["+str(selection_properties['Filter3'][i])+"]>"+str(selection_properties['LimitMag'][i])+": " + str(len(cone))

	print strftime("%Y-%m-%d %H:%M:%S", gmtime())

	"""
	fig = plt.figure()
	plt.title("Redshift vs Redder Magnitude")
	plt.xlabel("Apparent Redshift (Z_APP)")
	plt.ylabel("Magnitude in redder color (?)")
	plt.hist2d(cone['Z_APP'], cone[selection_properties['Filter3'][i]], bins=1000)
	plt.show()
	savemyplot("z_vs_mag")
	plt.close()
	"""

	"""
	Color selection
	"""

	# Points for the selection between the 3 colors
	# Bottom Left point (endH, limitH)
	# Top Right point (limitV, endV)
	"""
	#z=3
	limitH=1.0 
	limitV=1.2
	endH=0.15
	endV=2.5
	"""

	#z=4
	limitH=1.0 
	limitV=1.0
	endH=0.1
	endV=2.28



	# y = m x + p
	m = (limitH-endV) / (endH-limitV)
	p = limitH - m * endH

	# color differences (axes)
	f1minusf2 = cone.field(selection_properties['Filter1'][i]) - cone.field(selection_properties['Filter2'][i])
	f2minusf3 = cone.field(selection_properties['Filter2'][i]) - cone.field(selection_properties['Filter3'][i])

	print strftime("%Y-%m-%d %H:%M:%S", gmtime())
	mask_colorX = f2minusf3 < limitV
	mask_colorY = f1minusf2 > limitH
	mask_color_mp = f1minusf2 > m * f2minusf3 + p

	mask = mask_colorX & mask_colorY & mask_color_mp
	print strftime("%Y-%m-%d %H:%M:%S", gmtime())

	fig = plt.figure()
	plt.title("Colors for z~"+str(selection_properties['z'][i]))
	plt.xlabel(selection_properties['Filter2'][i] + "-" +selection_properties['Filter3'][i])
	plt.ylabel(selection_properties['Filter1'][i] + "-" +selection_properties['Filter2'][i])
	plt.hist2d(cone[selection_properties['Filter2'][i]] - cone[selection_properties['Filter3'][i]], cone[selection_properties['Filter1'][i]] - cone[selection_properties['Filter2'][i]], bins=150, range=([-1.,2.5],[-1.,8.5]))
	#plt.plot(cone[selection_properties['Filter2'][i]] - cone[selection_properties['Filter3'][i]], cone[selection_properties['Filter1'][i]] - cone[selection_properties['Filter2'][i]], '.')

	cone = cone[mask]
	print strftime("%Y-%m-%d %H:%M:%S", gmtime())

	plt.scatter(cone[selection_properties['Filter2'][i]] - cone[selection_properties['Filter3'][i]], cone[selection_properties['Filter1'][i]] - cone[selection_properties['Filter2'][i]],c=cone['Z_APP'],vmin=selection_properties['z'][i]-1,vmax=selection_properties['z'][i]+1,cmap=plt.cm.spectral)
	cbar = plt.colorbar()

############ label cbar ???

	cbar.set_label('z') 
	plt.plot([-1., endH], [limitH, limitH], '-b')
	plt.plot([limitV, limitV], [endV, 8.5], '-b')
	plt.plot([endH, limitV], [limitH, endV], '-b')

	plt.xlim(-1.,2.5) 
	plt.ylim(-1.,8.5) 
	plt.show()
	savemyplot("Colors_z_"+str(selection_properties['z'][i]))
	plt.close()


	print "Number of galaxies selected by color : "+str(len(cone))



	fig = plt.figure()
	plt.title("Redshift distribution for the selection @ z~"+str(selection_properties['z'][i]))
	plt.xlabel("Apparent Redshift (Z_APP)")
	plt.ylabel("#")
	plt.hist(cone['Z_APP'], bins=20)
	plt.show()
	savemyplot("z_dist_selection_z_"+str(selection_properties['z'][i]))
	plt.close()



#TODO : Make field of view to the right size







