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
from numpy.random import random
import math
from time import gmtime, strftime


plot_extension = ".png"
plot_directory = "./plots/"
def savemyplot(fig, name):
	fig.savefig(plot_directory+name+plot_extension)
	return

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


#########################
#### Fields of View #####
#########################
def info_FoV():
	
	global Ratio_FoV
	
	# PFS internal diameter of the circle inscrit in the hexagon:
	PFS_diameter = 1.3
	# PFS Field of View:
	PFS_FoV = 3./4. * PFS_diameter * PFS_diameter * math.cos(30.*np.pi/180.)

	# Lightcones FoV:
	Lightcones_side = 1.4
	Lightcones_FoV = Lightcones_side * Lightcones_side

	# FoV Ratio between PFS and the lightcones:
	Ratio_FoV = PFS_FoV / Lightcones_FoV


def open_lightcone(file_number):

	global allcone, cols

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
	savemyplot(fig, "z_dist")
	plt.close()
	"""
	
def creates_tables():

	global selection_properties, gaussian_selection, dz
	
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


	# Prepares a result table for the gaussian selection
	Nobj_zbin =      [       0,        0,        0,        0,        0,      0]
	Nobj_gauss =     [       0,        0,        0,        0,        0,      0]
	Nobj_PFS =       [       0,        0,        0,        0,        0,      0]
	Nobj_expected =  [    2700,     2000,      830,      190,       14,      4]
	gaussian_selection = Table([z_bins, mag_limit, Nobj_zbin, Nobj_gauss, Nobj_PFS, Nobj_expected], names=('z', 'LimitMag', '# objects in z bin', '# objects gaussian', '# objects PFS', '# expected objects'), meta={'name': 'table of gaussian selected objects'})


def main():

	info_FoV()

	creates_tables()

	file_number = 1
	open_lightcone(file_number)

	selec_gauss()
	#selec_3colors()



def selec_gauss():

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
		
		# I want to select a gaussian distribution from a random distribution, selecting all the objects at the peak of the gaussian.

		# Paremeters of the gaussian distribution:
		distribution_parameters = [1., selection_properties['z'][i], dz/2.5]
		
		nbins = 10
		
		# Current distribution:
		hist, bin_edges = np.histogram(cone['Z_APP'], bins = nbins, range =(selection_properties['z'][i]-dz,selection_properties['z'][i]+dz))

		# Probability of selection in order to inverse the distribution to a gaussian distribution:
		Pz = gauss(bin_edges[:-1], *distribution_parameters) / hist # This can lead to 0/0, but that's a minor point, considering it will result in an expected 0.
		Pz = Pz / np.min(Pz[0.4*len(Pz):0.6*len(Pz)]) # The min does a light over selection.

		# Objects are selected randomly depending on the probability Pz of the bin they belong to.
		mask_gaussian = np.array([], dtype=bool)
		for objects in cone:
			proba = Pz[(np.where(objects.field('Z_APP') <= bin_edges))[0][0]-1]
			mask_gaussian = np.append(mask_gaussian, [proba>random()])
		
		# gaussian distributed selection:
		cone_gaussian = cone[mask_gaussian]

		gaussian_selection['# objects in z bin'][i] = len(cone)
		gaussian_selection['# objects gaussian'][i] = len(cone_gaussian)
		gaussian_selection['# objects PFS'][i] = int(round(len(cone_gaussian)*Ratio_FoV))
		
		print "Total number of elements in the redshift bin : " + str(gaussian_selection['# objects in z bin'][i])
		print "Number of elements after gaussian selection : " + str(gaussian_selection['# objects gaussian'][i])
		print "Number of elements after FoV correction : " + str(gaussian_selection['# objects PFS'][i])
		
		fig = plt.figure()
		plt.title("Redshift distribution for the gaussian selection @ z~"+str(selection_properties['z'][i])+"\n Objects selected only inside PFS FoV: "+str(gaussian_selection['# objects PFS'][i]))
		plt.xlabel("Apparent Redshift (Z_APP)")
		plt.ylabel("#")
		plt.hist(cone['Z_APP'], bins=nbins,label="Initial distribution: "+str(gaussian_selection['# objects in z bin'][i]), range =(selection_properties['z'][i]-dz,selection_properties['z'][i]+dz))
		plt.hist(cone_gaussian['Z_APP'], bins=nbins,label="Gaussian selection: "+str(gaussian_selection['# objects gaussian'][i]), range =(selection_properties['z'][i]-dz,selection_properties['z'][i]+dz))
		#plt.plot(bin_edges[:-1], Pz)
		plt.plot(bin_edges, gauss(bin_edges, *[max(np.histogram(cone_gaussian['Z_APP'], bins=nbins, range =(selection_properties['z'][i]-dz,selection_properties['z'][i]+dz))[0]), distribution_parameters[1], distribution_parameters[2]]),label ="just a gaussian")
		#plt.plot(bin_edges[:-1], hist_densities)
		#plt.plot(bin_edges[:-1], gauss(bin_edges[:-1], *distribution_parameters))
		plt.legend()
		plt.xlim(selection_properties['z'][i]-dz, selection_properties['z'][i]+dz)
		savemyplot(fig, "z_dist_gaussian_selection_z_"+str(selection_properties['z'][i]))
	#	plt.show()
		plt.close()

	print gaussian_selection

	"""
	fig = plt.figure()
	plt.title("Comparison objects selected vs objects expected in the PFS FoV")
	plt.xlabel("Redshift bins")
	plt.ylabel("#")
	#plt.hist(cone['Z_APP'], bins=nbins,label="Initial distribution: "+str(gaussian_selection['# objects in z bin'][i]), range =(selection_properties['z'][i]-dz,selection_properties['z'][i]+dz))
	#plt.hist(cone_gaussian['Z_APP'], bins=nbins,label="Gaussian selection: "+str(gaussian_selection['# objects gaussian'][i]), range =(selection_properties['z'][i]-dz,selection_properties['z'][i]+dz))
	plt.plot(gaussian_selection['z'], gaussian_selection['# objects PFS'], label="# of selected objects")
	plt.plot(gaussian_selection['z'], gaussian_selection['# expected objects'], label="Expected # of selected objects")
	plt.legend()
	plt.yscale('log')
	savemyplot(fig, "Object_Counts")
	plt.show()
	plt.close()
	"""


def selec_3colors():

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

		# selecting all objects with z>2 takes 20 min
		cone = allcone[mask]

		print "Number of candidates with mag["+str(selection_properties['Filter3'][i])+"]>"+str(selection_properties['LimitMag'][i])+": " + str(len(cone))


		"""
		fig = plt.figure()
		plt.title("Redshift vs Redder Magnitude")
		plt.xlabel("Apparent Redshift (Z_APP)")
		plt.ylabel("Magnitude in redder color (?)")
		plt.hist2d(cone['Z_APP'], cone[selection_properties['Filter3'][i]], bins=1000)
		plt.show()
		savemyplot(fig, "z_vs_mag")
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

		mask_colorX = f2minusf3 < limitV
		mask_colorY = f1minusf2 > limitH
		mask_color_mp = f1minusf2 > m * f2minusf3 + p

		mask = mask_colorX & mask_colorY & mask_color_mp

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
		savemyplot(fig, "Colors_z_"+str(selection_properties['z'][i]))
		plt.close()


		print "Number of galaxies selected by color : "+str(len(cone))



		fig = plt.figure()
		plt.title("Redshift distribution for the 3 color selection @ z~"+str(selection_properties['z'][i]))
		plt.xlabel("Apparent Redshift (Z_APP)")
		plt.ylabel("#")
		plt.hist(cone['Z_APP'], bins=20)
		plt.show()
		savemyplot(fig, "z_dist_color_selection_z_"+str(selection_properties['z'][i]))
		plt.close()



if __name__ == '__main__':
  main()

