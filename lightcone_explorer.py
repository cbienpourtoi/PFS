"""
File: lightcone_explorer.py
Copyright: Loic Le Tiran, 2014
Contact: loic.le-tiran@obspm.fr
Licence: GNU GPL v3

Description:
Tool for reading the lightcones from 
http://galformod.mpa-garching.mpg.de/qa/mrobs/pages/surveys/PFS.jsp

Requirements:
imagemagic problem ?

TODO:
Patch NP that happears to be negative sometimes when reading the fits files. Pass int16 as uint16.

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
from mpl_toolkits.basemap import Basemap
import matplotlib.animation as animation
from numpy.lib.recfunctions import append_fields


plot_extension = ".png"
plot_directory = "./plots/"
def savemyplot(fig, name):
	fig.savefig(plot_directory+name+plot_extension)
	return

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))



#########################
####      Main      #####
#########################
def main():

	info_FoV()

	creates_tables()

	file_number = 1
	open_lightcone(file_number)

	#selec_gauss()
	selec_3colors()

	plot_sky()






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


##########################
#### Opens Lightcones ####
##########################
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

	global selection_properties, gaussian_selection, color_selection, dz
	
	"""
	Table of redshift bins, properties, limit magnitudes and selection filters
	"""

	dz = 0.5 # Delta in redshift (ie selection between z-dz and z+dz)


	# Make a table of the redshifts, selection filters and limit magnitude associated (limit mag always on the redder)
	z_bins =         [      2.,       3.,       4.,       5.,       6.,     7.]
	mag_limit =      [     24.,     24.3,     24.5,     24.9,     24.9,   25.3]
	selec_filter_1 = [      '', 'SDSS_U', 'SDSS_G', 'SDSS_R', 'SDSS_I',    'Z']
	selec_filter_2 = [      '', 'SDSS_G', 'SDSS_R', 'SDSS_I', 'SDSS_Z',    'Y']
	selec_filter_3 = ['SDSS_G', 'SDSS_R', 'SDSS_I', 'SDSS_Z',      'Y',    'J']
	# Points for the selection between the 3 colors
	# Bottom Left point (endH, limitH)
	# Top Right point (limitV, endV)
	limitH =         [      1.,       1.,     1.12,      1.1,       1.,     1.]
	limitV =         [      1.,      1.2,       1.,     0.82,       1.,     1.]
	endH =           [     0.1,      .15,     0.04,     0.24,      0.1,   0.05]
	endV =           [    2.28,      2.5,     2.45,     1.58,     2.28,   2.28]
	selection_properties = Table([z_bins, mag_limit, selec_filter_1, selec_filter_2, selec_filter_3, limitH, limitV, endH, endV], names=('z', 'LimitMag', 'Filter1', 'Filter2', 'Filter3', 'selec: limitH', 'selec: limitV', 'selec: endH', 'selec: endV'), meta={'name': 'table of the selection properties'})


	# Prepares a result table for the gaussian selection
	Nobj_zbin =      [       0,        0,        0,        0,        0,      0]
	Nobj_gauss =     [       0,        0,        0,        0,        0,      0]
	Nobj_PFS =       [       0,        0,        0,        0,        0,      0]
	Nobj_expected =  [    2700,     2000,      830,      190,       14,      4]
	gaussian_selection = Table([z_bins, mag_limit, selec_filter_3, Nobj_zbin, Nobj_gauss, Nobj_PFS, Nobj_expected], names=('z', 'LimitMag', 'Filter3', '# objects in z bin', '# objects gaussian', '# objects PFS', '# expected objects'), meta={'name': 'table of gaussian selected objects'})

	# Prepares a result table for the 3 colors selection
	Nobj_maglim =    [       0,        0,        0,        0,        0,      0]
	Nobj_3colors =   [       0,        0,        0,        0,        0,      0]
	Nobj_PFS =       [       0,        0,        0,        0,        0,      0]
	Nobj_expected =  [    2700,     2000,      830,      190,       14,      4]
	color_selection = Table([z_bins, mag_limit, selec_filter_1, selec_filter_2, selec_filter_3, Nobj_maglim, Nobj_3colors, Nobj_PFS, Nobj_expected], names=('z', 'LimitMag', 'Filter1', 'Filter2', 'Filter3', '# objects under mag lim', '# objects color selected', '# objects PFS', '# expected objects'), meta={'name': 'table of 3 colors selected objects'})



	
def plot_sky():
	
	global zi, dz_plot
	global lllon, lllat, urlon, urlat

	# Infos for selecting redshift slices
	dz_plot = 0.015
	zmin = 3.5
	zmax = 8.
	zi = np.arange(zmin, zmax, dz_plot*2)

	# Infos for the positions of the corners of the basemap
	lllon=min(allcone.field('RA'))
	lllat=min(allcone.field('Dec'))
	urlon=max(allcone.field('RA'))
	urlat=max(allcone.field('Dec'))
	
	fig = plt.figure(figsize=(10,10))  
	anim = animation.FuncAnimation(fig, animate, frames=25)
	anim.save('animation.gif', writer='imagemagick', fps = 2);
	#plt.show()
	"""
	
	fig = plt.figure()
	#m = Basemap(projection='merc',lon_0=0, lat_0=0, celestial=True)
	ani = animation.FuncAnimation(fig, animate, frames = len(zi), interval=50, blit=True)
	#ani.save('animation.gif', writer='imagemagick', fps = 4);
	plt.show()
	"""


def animate(nframe):
	print str(nframe)+'/'+str(len(zi))
	
	# Selects the data in the redshift slice
	mask = np.where(np.abs(allcone.field('Z_APP') - zi[nframe]) < dz_plot)
	conedz = allcone[mask]
	lats = conedz.field('Dec')
	lons = conedz.field('RA')

	# Selects the data selected in the previous selection 
	lats_3colors = np.array([])
	lons_3colors = np.array([])
	global common_GALID
	common_GALID = set(conedz.field('GALID')) & set(list_GALID)
	for ids in common_GALID:
		lats_3colors = np.append(lats_3colors, conedz[np.where(conedz.field('GALID') == ids)].field('Dec') )
		lons_3colors = np.append(lons_3colors, conedz[np.where(conedz.field('GALID') == ids)].field('RA') )
	

	plt.cla()
	m = Basemap(projection='merc',lon_0=0, lat_0=0, llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat, celestial=True)

	# Lattitudes and longtitudes
	poslines = [-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]
	m.drawparallels(poslines,labels=[1,0,0,0])
	m.drawmeridians(poslines,labels=[0,0,0,1])

	# draw points
	x, y = m(lons,lats)
	m.scatter(x,y,0.03,marker='o',color='b')
	x_3colors, y_3colors = m(lons_3colors,lats_3colors)
	m.scatter(x_3colors, y_3colors, 10,marker='o',color='r')

	# Adds a title
	plt.title('z='+str(zi[nframe]))

	
"""
def animate(nframe):

	print nframe
	
	plt.cla()
	
	mask = np.where(np.abs(allcone.field('Z_APP') - zi[nframe]) < dz_plot)
	conedz = allcone[mask]

	#conedz_3colors = conedz[np.where(allcone_selected_3colors == True)]

	lats = conedz.field('Dec')
	lons = conedz.field('RA')
	
	lats_3colors = np.array([])
	lons_3colors = np.array([])
	global common_GALID
	common_GALID = set(conedz.field('GALID')) & set(list_GALID)
	for ids in common_GALID:
		lats_3colors = np.append(lats_3colors, conedz[np.where(conedz.field('GALID') == ids)].field('Dec') )
		lons_3colors = np.append(lons_3colors, conedz[np.where(conedz.field('GALID') == ids)].field('RA') )
		
	
	#for ids in list_GALID:
	#	print np.where(conedz.field('GALID') == ids)
	#lats_3colors = conedz_3colors.field('Dec')
	#lons_3colors = conedz_3colors.field('RA')
	
	# Classical coordinate system:
	#lats = random(10) * 180. - 90.
	#lons = random(10) * 360.

	#m = Basemap(projection='merc',lon_0=0, lat_0=0, llcrnrlon=min(allcone.field('RA')), llcrnrlat=min(allcone.field('Dec')), urcrnrlon=max(allcone.field('RA')), urcrnrlat=max(allcone.field('Dec')), celestial=True)
	
	# draw map with markers for float locations
	#x, y = m(lons,lats)
	#x_3colors, y_3colors = m(lons_3colors,lats_3colors)
	#xz, yz = m(0.65, 0.65)
	#poslines = [-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]
	#m.drawparallels(poslines,labels=[1,0,0,0])
	#m.drawmeridians(poslines,labels=[0,0,0,1])
	
	#plt.title('Map')
	plt.title('z='+str(zi[nframe]))
	
	#plt.title('z='+str(zi[nframe]))
	#m.scatter(x,y,0.03,marker='o',color='b')
	#m.scatter(x_3colors, y_3colors, 10,marker='o',color='r')
	#t = plt.text(xz, yz, 'z='+str(zi[nframe]))
	
	
	#print x_3colors, y_3colors, str(zi[nframe])

	#return pts#points, points2, pts#, t2
	
"""

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

	global conelist, cone, list_GALID, allcone_selected_3colors

	print strftime("%Y-%m-%d %H:%M:%S", gmtime())
	allcone_selected_3colors = np.zeros(len(allcone), dtype=int)
	print strftime("%Y-%m-%d %H:%M:%S", gmtime())

	number_duplicates = 0

	bins = [1,2,3,4,5]
	#bins = [1]
	#bins = np.arange(len(selection_properties))

	conelist = []
	list_GALID = []

	for i in bins:
		print "redshift: ~" + str(selection_properties['z'][i]) + ". Filters : " + str(selection_properties['Filter1'][i]) +" "+ str(selection_properties['Filter2'][i]) +" "+ str(selection_properties['Filter3'][i])
		#mask_z = np.abs( allcone.field('Z_APP') - selection_properties['z'][i] ) < dz
		mask_mag = allcone.field(selection_properties['Filter3'][i]) < selection_properties['LimitMag'][i]

		#print "size mask_mag : "+ str(len(mask_mag))
		#print "number True :"+ str(np.count_nonzero(mask_mag))

		cone = allcone[mask_mag]
		color_selection['# objects under mag lim'][i] = len(cone)

		print "Number of candidates with mag["+str(selection_properties['Filter3'][i])+"]>"+str(selection_properties['LimitMag'][i])+": " + str(color_selection['# objects under mag lim'][i])


		"""
		fig = plt.figure()
		plt.title("Redshift vs Redder Magnitude, selection props for z~"+str(selection_properties['z'][i]))
		plt.xlabel("Apparent Redshift (Z_APP)")
		plt.ylabel("Magnitude in redder color: "+str(selection_properties['Filter3'][i]) )
		plt.hist2d(cone['Z_APP'], cone[selection_properties['Filter3'][i]], bins=1000)
		#plt.show()
		savemyplot(fig, "z_vs_mag_for_z_"+str(selection_properties['z'][i]))
		plt.close()
		"""
		
		fig = plt.figure()
		plt.title("ZOOM: Number of particles for "+str(selection_properties['Filter3'][i])+" < "+str(selection_properties['LimitMag'][i]))
		plt.xlabel("NP")
		plt.ylabel("#")
		#plt.yscale('log')
		plt.hist(cone['NP'],bins = 1000, range = (0,1000))
		#plt.show()
		savemyplot(fig, "NP_ZOOM_filter3_le_"+str(selection_properties['LimitMag'][i]))
		plt.close()
		
		
		fig = plt.figure()
		plt.title("Number of particles for "+str(selection_properties['Filter3'][i])+" < "+str(selection_properties['LimitMag'][i]))
		plt.xlabel("NP")
		plt.ylabel("#")
		#plt.yscale('log')
		plt.hist(cone['NP'],bins = 1000)
		plt.show()
		savemyplot(fig, "NP_filter3_le_"+str(selection_properties['LimitMag'][i]))
		plt.close()
		
		
		"""
		fig = plt.figure()
		plt.title("NP vs Redder Magnitude, selection props for z~"+str(selection_properties['z'][i]))
		plt.xlabel("Number of particles (NP)")
		plt.ylabel("Magnitude in redder color: "+str(selection_properties['Filter3'][i]) )
		plt.hist2d(cone['NP'], cone[selection_properties['Filter3'][i]], bins=1000)
		plt.show()
		savemyplot(fig, "NP_vs_mag_for_z_"+str(selection_properties['z'][i]))
		plt.close()
		"""
		
		###########################################
		#    3 color selection is done here:      #
		###########################################
		
		limitH = selection_properties['selec: limitH'][i]
		limitV = selection_properties['selec: limitV'][i]
		endH = selection_properties['selec: endH'][i]
		endV = selection_properties['selec: endV'][i]

		# y = m x + p
		m = (limitH-endV) / (endH-limitV)
		p = limitH - m * endH

		f1minusf2 = cone.field(selection_properties['Filter1'][i]) - cone.field(selection_properties['Filter2'][i])
		f2minusf3 = cone.field(selection_properties['Filter2'][i]) - cone.field(selection_properties['Filter3'][i])
		
		mask_colorX = f2minusf3 < limitV
		mask_colorY = f1minusf2 > limitH
		mask_color_mp = f1minusf2 > m * f2minusf3 + p
		mask = mask_colorX & mask_colorY & mask_color_mp

		#print "size mask : "+ str(len(mask))

		# Color diagram : histogram + individual points for selected objects
		fig = plt.figure()
		plt.title("Colors for z~"+str(selection_properties['z'][i]))
		plt.xlabel(selection_properties['Filter2'][i] + "-" +selection_properties['Filter3'][i])
		plt.ylabel(selection_properties['Filter1'][i] + "-" +selection_properties['Filter2'][i])

		# histogram
		plt.hist2d(cone[selection_properties['Filter2'][i]] - cone[selection_properties['Filter3'][i]], cone[selection_properties['Filter1'][i]] - cone[selection_properties['Filter2'][i]], bins=150, range=([-1.,2.5],[-1.,8.5]))
		#plt.plot(cone[selection_properties['Filter2'][i]] - cone[selection_properties['Filter3'][i]], cone[selection_properties['Filter1'][i]] - cone[selection_properties['Filter2'][i]], '.')

		# selecting objects
		cone = cone[mask]
		print "Number of objects after 3 colors selection: " + str(len(cone))

		# Marks selected objects in the initial cone: 
		#allcone_selected_3colors[mask_mag] = 1
		#allcone_selected_3colors[np.where(allcone_selected_3colors == 1)] = 1
		#print len(np.where(allcone_selected_3colors == True))
		
		# Checking that some of these objects are not duplicates
		# If there are duplicates, they are deleted from this new bin (ie they stay in the precedent lower z bin).
		duplicates = set(list_GALID) & set(cone.field('GALID'))
		print "Number of duplicates: " + str(len(duplicates))
		if len(duplicates) > 0:
			number_duplicates = number_duplicates + len(duplicates)
			mask_duplicates = np.ones(len(cone), dtype=bool)
			for dupie in duplicates:
				mask_duplicates[np.where(cone.field('GALID') == dupie)] = False
			cone = cone[mask_duplicates]
			print "After deleting the duplicates in the new cone, number of objects in the cone: "+str(len(cone)) 

		color_selection['# objects color selected'][i] = len(cone)
		print "Number of galaxies selected by color : "+str(color_selection['# objects color selected'][i])

		for ids in cone.field('GALID'):
			list_GALID.append(ids)

		# Adding this cone bin to a list of all "cones" at different redshifts
		conelist.append(cone)
		
		
		# plotting individual points for selected objects
		plt.scatter(cone[selection_properties['Filter2'][i]] - cone[selection_properties['Filter3'][i]], cone[selection_properties['Filter1'][i]] - cone[selection_properties['Filter2'][i]],c=cone['Z_APP'],vmin=selection_properties['z'][i]-1,vmax=selection_properties['z'][i]+1,cmap=plt.cm.spectral)
		cbar = plt.colorbar()
		cbar.set_label('redshift') 

		# plotting the limits
		plt.plot([-1., endH], [limitH, limitH], '-b')
		plt.plot([limitV, limitV], [endV, 8.5], '-b')
		plt.plot([endH, limitV], [limitH, endV], '-b')

		if selection_properties['z'][i] == 3:
			plt.xlim(-1.,2.5) 
			plt.ylim(-1.,8.5) 
		else:
			plt.xlim(-0.5,1.5) 
			plt.ylim(-1.,4.2) 
		

		#plt.show()
		savemyplot(fig, "Colors_z_"+str(selection_properties['z'][i]))
		plt.close()


		fig = plt.figure()
		plt.title("ZOOM after 3 color selection: Number of particles @ z~"+str(selection_properties['z'][i]))
		plt.xlabel("NP")
		plt.ylabel("#")
		#plt.yscale('log')
		plt.hist(cone['NP'],bins = 1000, range = (0,1000))
		#plt.show()
		savemyplot(fig, "NP_ZOOM_after_3c_selection_filter3_z_"+str(selection_properties['z'][i]))
		plt.close()

		fig = plt.figure()
		plt.title("After 3 color selection: Number of particles @ z~"+str(selection_properties['z'][i]))
		plt.xlabel("NP")
		plt.ylabel("#")
		#plt.yscale('log')
		plt.hist(cone['NP'])
		#plt.show()
		savemyplot(fig, "NP_after_3c_selection_filter3_z_"+str(selection_properties['z'][i]))
		plt.close()

		fig = plt.figure()
		plt.title("Redshift distribution for the 3 color selection @ z~"+str(selection_properties['z'][i]))
		plt.xlabel("Apparent Redshift (Z_APP)")
		plt.ylabel("#")
		plt.hist(cone['Z_APP'], bins=40)
		#plt.show()
		savemyplot(fig, "z_dist_color_selection_z_"+str(selection_properties['z'][i]))
		plt.close()
		
		color_selection['# objects PFS'][i] = int(round(color_selection['# objects color selected'][i]*Ratio_FoV))
		
	fig = plt.figure()
	plt.title("Redshift distribution for the 3 color selection")
	plt.xlabel("Apparent Redshift (Z_APP)")
	plt.ylabel("#")
	for conei in conelist:
		plt.hist(conei['Z_APP'], bins=40, range = (0,9), histtype = 'step', label="$z_{median} \sim"+str("%0.1f"%(np.median(conei['Z_APP'])))+"$ "+str("%6.f"%len(conei))+" objects")
	plt.legend()
	#plt.show()
	savemyplot(fig, "z_dist_color_selection")
	plt.close()

	print color_selection

	if number_duplicates !=0:
		print "There was "+str(number_duplicates)+" duplicates in the selection. They have been taken care of."


	
	



if __name__ == '__main__':
  main()

