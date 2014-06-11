"""
File: lightcone_explorer.py
Copyright: Loic Le Tiran, 2014
Contact: loic.le-tiran@obspm.fr
Licence: GNU GPL v3

Description:
Tool for reading the lightcones from 
http://galformod.mpa-garching.mpg.de/qa/mrobs/pages/surveys/PFS.jsp
"""


#TODO : Make field of view to the right size



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
from numpy.random import random
#import math
from time import gmtime, strftime


plot_extension = ".png"
plot_directory = "./plots/"
def savemyplot(name):
	fig.savefig(plot_directory+name+plot_extension)
	return

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))






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










print "##################################################"
print "#######      Selection method 2:                 #" 
print "#######  just selecting statistically in z bins  #"
print "##################################################"

for i in np.arange(len(selection_properties)):

	print "Redshift: ~" + str(selection_properties['z'][i]) + "; Limit Magnitude: " + str(selection_properties['LimitMag'][i])
	mask_z = np.abs( allcone.field('Z_APP') - selection_properties['z'][i] ) < dz
	mask_mag = allcone.field(selection_properties['Filter3'][i]) < selection_properties['LimitMag'][i]
	mask = mask_mag & mask_z
	print strftime("%Y-%m-%d %H:%M:%S", gmtime())

	# selecting all objects with z>2 takes 20 min
	cone = allcone[mask]
	
	# I want to select a gaussian distribution.
	# I will take for each object a random number "choicemaker" 
	# between 0 and 1 (using random which I checked to be uniform)
	# Each z corresponds to a probability Pz (1 for z=z_center, eg z=2, 
	# and decreases acccording to a gaussian around).
	# If choicemaker<Pz, I keep the object
	#choicemaker = random(len(cone))
	distribution_parameters = [1., selection_properties['z'][i], dz/2.5]
	
	nbins = len(cone)/50
	hist, bin_edges = np.histogram(cone['Z_APP'], bins = nbins)
	#hist_densities, bin_edges = np.histogram(cone['Z_APP'], bins = 20, density = True)

	
	Pz = gauss(bin_edges[:-1], *distribution_parameters) / hist # AJOUTER NORMALISATION A 1 SUR z=2
	Pz = Pz / np.mean(Pz[0.45*len(Pz):0.55*len(Pz)])

	mask_gaussian = np.array([], dtype=bool)
	for objects in cone:
		proba = Pz[(np.where(objects.field('Z_APP') <= bin_edges))[0][0]-1]
		mask_gaussian = np.append(mask_gaussian, [proba>random()])
	
	cone_gaussian = cone[mask_gaussian]
	
	print "Number of elements in the redshift bin : " + str(len(cone))
	print "Number of elements after gaussian selection : " + str(len(cone_gaussian))
	
	fig = plt.figure()
	plt.title("Redshift distribution for the fake selection @ z~"+str(selection_properties['z'][i]))
	plt.xlabel("Apparent Redshift (Z_APP)")
	plt.ylabel("#")
	plt.hist(cone['Z_APP'], bins=nbins)
	plt.hist(cone_gaussian['Z_APP'], bins=nbins)
	plt.plot(bin_edges[:-1], Pz)
	plt.plot(bin_edges, gauss(bin_edges, *[38., distribution_parameters[1], distribution_parameters[2]]))
	#plt.plot(bin_edges[:-1], hist_densities)
	#plt.plot(bin_edges[:-1], gauss(bin_edges[:-1], *distribution_parameters))
	savemyplot("z_dist_fake_selection_z_"+str(selection_properties['z'][i]))
	plt.show()
	plt.close()

	
sys.exit()
	
























print "##################################################"
print "#######      Selection method 1:                 #" 
print "#######    real 3colors selction                 #"
print "##################################################"


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
	plt.title("Redshift distribution for the 3 color selection @ z~"+str(selection_properties['z'][i]))
	plt.xlabel("Apparent Redshift (Z_APP)")
	plt.ylabel("#")
	plt.hist(cone['Z_APP'], bins=20)
	plt.show()
	savemyplot("z_dist_color_selection_z_"+str(selection_properties['z'][i]))
	plt.close()




	


