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
Patch NP that appears to be negative sometimes when reading the fits files. Pass int16 as uint16.

"""

import numpy as np
import scipy
# import pyfits
import sys
import os
import glob
import csv
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, vstack
#from astropy.cosmology.parameters import WMAP9 # depreciated since astropy 0.4
# from astropy.cosmology import comoving_distance # depreciated since astropy 0.4
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
from astropy.cosmology import z_at_value
# from string import upper,lower
# from numpy import recfromcsv
#from numpy.random import normal
from numpy.random import random
import math
from time import gmtime, strftime
from mpl_toolkits.basemap import Basemap
import matplotlib.animation as animation
from numpy.lib.recfunctions import append_fields


table_write_format = 'fixed_width'


def savemyplot(fig, name):
    fig.savefig(plot_directory + name + plot_extension)
    return


def gauss(x, *p):
    A, mu, sigma = p
    return A * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))


#########################
####      Main      #####
#########################
def main():

    global plot_extension
    plot_extension = ".png"



    creates_tables()

    #open_CFHTLS(1)
    #selec_3colors_CFHTLS()
    #sys.exit()


    file_number = 1
    open_lightcone(file_number)

    #test_z2()
    #sys.exit()

    selec_gauss()
    #selec_3colors()

    look_overdense()
    sys.exit()

    #plot_sky_animate()

    plot_sky()


#########################
#### Fields of View #####
#########################
def info_FoV(Lightcones_FoV):
    global Ratio_FoV, Ratio_FoV_PFS_CFHTLS

    # PFS internal diameter of the circle inscrit in the hexagon:
    PFS_diameter = 1.3  #deg
    # PFS Field of View:
    PFS_FoV = 3. / 4. * PFS_diameter * PFS_diameter * math.cos(30. * np.pi / 180.)


    # FoV Ratio between PFS and the lightcones:
    Ratio_FoV = PFS_FoV / Lightcones_FoV


    # CFHTLS - each FoV is 1x1 sq deg
    CFHTLS_FoV = 1. * 1.  #sq deg
    Ratio_FoV_PFS_CFHTLS = PFS_FoV / CFHTLS_FoV


##########################
#### Opens Lightcones ####
##########################
def open_lightcone(file_number):
    global allcone, cols, plot_directory, galid, conename
    global lllon, lllat, urlon, urlat

    # Lightcones path
    conepath = "./data/lightcones/"
    conename = "wmap1_bc03_" + str(file_number).zfill(3) + "_igm1.fits"
    #conename = "P1_M05_Ks28.fits"
    plot_directory = "./plots/"+conename[0:-5]+"/"
    if not os.path.exists(plot_directory) : os.mkdir(plot_directory)

    hdulist = fits.open(conepath + conename)
    allcone = hdulist[1].data
    cols = hdulist[1].columns
    # Some files have keyword GALAXYID, others have GALID: we check here what we want:
    if 'GALAXYID' in str(cols):
        galid = 'GALAXYID'
        # Lightcones FoV:
        Lightcones_radius = 1.  #deg
        Lightcones_FoV = np.pi * Lightcones_radius * Lightcones_radius
        lllon, lllat, urlon, urlat = -1., -1., 1., 1.
    elif 'GALID' in str(cols):
        galid = 'GALID'
        # Lightcones FoV:
        Lightcones_side = 1.4  #deg
        Lightcones_FoV = Lightcones_side * Lightcones_side
        lllon, lllat, urlon, urlat = -0.7, -0.7, 0.7, 0.7
    else:
        print "Cannot find either GALID or GALAXYID in the keywords: exiting savagely."
        sys.exit()

    hdulist.close()

    print "There are " + str(len(allcone)) + " objects in the cone."

    info_FoV(Lightcones_FoV)



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


##################################
#### Makes Properties Tables  ####
##################################
def creates_tables():
    global selection_properties, gaussian_selection, color_selection, color_selection_CFHTLS, dz

    """
    Table of redshift bins, properties, limit magnitudes and selection filters
    """

    dz = 0.5  # Delta in redshift (ie selection between z-dz and z+dz)


    # Make a table of the redshifts, selection filters and limit magnitude associated (limit mag always on the redder)
    z_bins = [2., 3., 4., 5., 6., 7.]
    mag_limit = [24., 24.3, 24.5, 24.9, 24.9, 25.3]

    selec_filter_1 = ['', 'SDSS_U', 'SDSS_G', 'SDSS_R', 'SDSS_I', 'Z']
    selec_filter_2 = ['', 'SDSS_G', 'SDSS_R', 'SDSS_I', 'SDSS_Z', 'Y']
    selec_filter_3 = ['SDSS_G', 'SDSS_R', 'SDSS_I', 'SDSS_Z', 'Y', 'J']

    #CFHTLS
    selec_filter_1_CFHTLS = ['', 'u', 'g', 'r', 'i', 'z']
    selec_filter_2_CFHTLS = ['', 'g', 'r', 'i', 'z', '']
    selec_filter_3_CFHTLS = ['g', 'r', 'i', 'z', '', '']

    # Points for the selection between the 3 colors
    # Bottom Left point (endH, limitH)
    # Top Right point (limitV, endV)
    limitH = [1., 1., 1.12, 1.1, 1., 1.]
    limitV = [1., 1.2, 1., 0.82, 1., 1.]
    endH = [0.1, .15, 0.04, 0.24, 0.1, 0.05]
    endV = [2.28, 2.5, 2.45, 1.58, 2.28, 2.28]
    selection_properties = Table \
        ([z_bins, mag_limit, selec_filter_1, selec_filter_2, selec_filter_3, selec_filter_1_CFHTLS,
          selec_filter_2_CFHTLS, selec_filter_3_CFHTLS, limitH, limitV, endH, endV], names=
         ('z', 'LimitMag', 'Filter1', 'Filter2', 'Filter3', 'Filter1_CFHTLS', 'Filter2_CFHTLS', 'Filter3_CFHTLS',
          'selec: limitH', 'selec: limitV', 'selec: endH', 'selec: endV'),
         meta={'name': 'table of the selection properties'})  # Prepares a result table for the gaussian selection
    Nobj_zbin = [0, 0, 0, 0, 0, 0]
    Nobj_gauss = [0, 0, 0, 0, 0, 0]
    Nobj_PFS = [0, 0, 0, 0, 0, 0]
    Nobj_expected = [2700, 2000, 830, 190, 14, 4]
    gaussian_selection = Table([z_bins, mag_limit, selec_filter_3, Nobj_zbin, Nobj_gauss, Nobj_PFS, Nobj_expected],
                               names=
                               ('z', 'LimitMag', 'Filter3', '# objects in z bin', '# objects gaussian', '# objects PFS',
                                '# expected objects'), meta={'name': 'table of gaussian selected objects'})

    # Prepares a result table for the 3 colors selection
    Nobj_maglim = [0, 0, 0, 0, 0, 0]
    Nobj_3colors = [0, 0, 0, 0, 0, 0]
    Nobj_PFS = [0, 0, 0, 0, 0, 0]
    Nobj_expected = [2700, 2000, 830, 190, 14, 4]
    color_selection = Table \
        ([z_bins, mag_limit, selec_filter_1, selec_filter_2, selec_filter_3, Nobj_maglim, Nobj_3colors, Nobj_PFS,
          Nobj_expected], names=
         ('z', 'LimitMag', 'Filter1', 'Filter2', 'Filter3', '# objects under mag lim', '# objects color selected',
          '# objects PFS', '# expected objects'), meta={'name': 'table of 3 colors selected objects'})

    # Prepares a result table for the 3 colors selection for CFHTLS
    color_selection_CFHTLS = Table \
        ([z_bins, mag_limit, selec_filter_1_CFHTLS, selec_filter_2_CFHTLS, selec_filter_3_CFHTLS, Nobj_maglim,
          Nobj_3colors, Nobj_PFS, Nobj_expected], names=
         ('z', 'LimitMag', 'Filter1_CFHTLS', 'Filter2_CFHTLS', 'Filter3_CFHTLS', '# objects under mag lim CFHTLS',
          '# objects color selected CFHTLS', '# objects PFS CFHTLS', '# expected objects'),
         meta={'name': 'CFHTLS: table of 3 colors selected objects'})



################################
#### look for overdensities ####
################################
def look_overdense():
    # Sizes of bins: spatial, redshift
    xybin = 5./60. #deg = 5arcsec

    print "XYBIN should change with z !!!!!!"

    print "CHANGE VALUE OF LBIN FOR 10, ! 1000 IS JUST A TEST VALUE!"
    lbin = 1000 * u.Mpc

    # Initial redshift
    zdown = 1.5

    # Cut a redshift slice of depth lbin (converted in z)
    zup = z_at_value(cosmo.comoving_distance, cosmo.comoving_distance(zdown) + lbin)
    mask_down = selection.field('Z_APP') > zdown
    mask_up = selection.field('Z_APP') <= zup
    mask_slice = mask_down & mask_up
    cone_slice = selection[mask_slice]

    n_objects_slice = len(cone_slice)
    print n_objects_slice





    # HAVE TO CHANGE RA BECAUSE IT MUST BE 360 ! FOR THE OTHER FILE !

    # selects only the objects in teh right RA-Dec square
    xsize = urlon - lllon
    ysize = urlat - lllat

    nbinsx = xsize/xybin
    nbinsy = ysize/xybin

    density = np.empty([nbinsx, nbinsy])

    nx = 0
    ny = 0

    for nx in np.arange(0, nbinsx-1, 1):
        for ny in np.arange(0, nbinsy-1, 1):

            mask_RA_down = cone_slice.field('RA') < urlon - nx * xybin
            mask_RA_up = cone_slice.field('RA') >= urlon - (nx+1) * xybin
            mask_RA = mask_RA_down & mask_RA_up
            mask_Dec_down = cone_slice.field('DEC') > lllat + ny * xybin
            mask_Dec_up = cone_slice.field('DEC') <= lllat + (ny+1) * xybin
            mask_Dec = mask_Dec_down & mask_Dec_up
            mask_coord = mask_RA & mask_Dec
            cone_cube = cone_slice[mask_coord]
            n_objects_cube = len(cone_cube)
            density[nx][ny] = n_objects_cube / float(n_objects_slice)

            print n_objects_cube
            print density[nx][ny]

    fig = plt.figure(figsize=(6, 3.2))

    ax = fig.add_subplot(111)
    ax.set_title('colorMap')
    plt.imshow(density)
    ax.set_aspect('equal')

    plt.colorbar(orientation='vertical')
    plt.show()

    x=2

    sys.exit()















    print cosmo.comoving_distance([1,2,3,4])

    z=np.arange(0.5, 8.5, 0.1)
    comodist = cosmo.comoving_distance(z+.002) - cosmo.comoving_distance(z-.002)


    fig = plt.figure(figsize=(10, 10))
    #plt.title(conename)
    #plt.xkcd()
    plt.plot(z, comodist, label="Comobile Distance z-.002 to z+.002")
    plt.plot(z, 1./(1.+z)*20., label="prop 1/(1+z)")
    plt.legend()
    plt.xlabel("z")
    plt.show()
    #savemyplot(fig, "comobile")
    plt.close()




###################
#### Plot Sky  ####
###################
def plot_sky():
    global zi, dz_plot

    # Infos for the positions of the corners of the basemap lllon= min(allcone.field('RA'))


    fig = plt.figure(figsize=(10, 10))


    '''
    # Selects the data in the redshift slice
    mask = np.where(np.abs(allcone.field('Z_APP') - zi[nframe]) < dz_plot)
    conedz = allcone[mask]
    lats = conedz.field('Dec')
    lons = conedz.field('RA')
    '''

    # Selects the data selected in the previous selection
    lats_selection = selection.field('DEC')
    lons_selection = selection.field('RA')
    z_selection = selection.field('Z_APP')

    lons_selection[np.where(lons_selection > 180.)] -= 360.


    #plt.cla()
    m = Basemap(projection='merc', lon_0=0, lat_0=0, llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat, celestial=True)  # Lattitudes and longtitudes
    poslines = [-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]
    m.drawparallels(poslines, labels=[1, 0, 0, 0])
    m.drawmeridians(poslines, labels=[0, 0, 0, 1])
    plt.title(conename)

    # draw points
    #x, y = m(lons,lats)
    #m.scatter(x,y,0.03,marker='o',color='b')
    #x_selection, y_selection = m(lons_selection[np.where(selection.field('Z_APP') < 2.)], lats_selection[np.where(selection.field('Z_APP') < 2.)])
    #m.scatter(x_selection, y_selection, 10, marker='o', color='r')
    x_selection, y_selection = m(lons_selection, lats_selection)

    print x_selection
    print y_selection

    m.scatter(x_selection, y_selection, 10, marker='o', c=np.round(z_selection), vmin=min(z_selection), vmax=max(z_selection),cmap=plt.cm.spectral, lw = 0)
    #m.scatter(x_selection, y_selection, 10, marker='o')
    cbar = plt.colorbar()
    cbar.set_label('redshift')
    # Adds a title
    #plt.title('z='+str(zi[nframe]))


    plt.show()
    savemyplot(fig, "sky_map")
    plt.close()
    """
    fig = plt.figure()
    #m = Basemap(projection='merc',lon_0=0, lat_0=0, celestial=True)
    ani = animation.FuncAnimation(fig, animate, frames = len(zi), interval=50, blit=True)
    #ani.save('animation.gif', writer='imagemagick', fps = 4);
    plt.show()
    """


###############################
#### Plot Sky + animation  ####
###############################
def plot_sky_animate():
    global zi, dz_plot
    global lllon, lllat, urlon, urlat

    # Infos for selecting redshift slices
    dz_plot = 0.015
    zmin = 3.5
    zmax = 8.
    zi = np.arange(zmin, zmax, dz_plot * 2)

    # Infos for the positions of the corners of the basemap


    lllon = min(allcone.field('RA'))
    lllat = min(allcone.field('Dec'))
    urlon = max(allcone.field('RA'))
    urlat = max(allcone.field('Dec'))

    fig = plt.figure(figsize=(10, 10))
    anim = animation.FuncAnimation(fig, animate, frames=25)
    anim.save(plot_directory+'animation.gif', writer='imagemagick', fps=2);
    #plt.show()

    """
    fig = plt.figure()
    #m = Basemap(projection='merc',lon_0=0, lat_0=0, celestial=True)
    ani = animation.FuncAnimation(fig, animate, frames = len(zi), interval=50, blit=True)
    #ani.save('animation.gif', writer='imagemagick', fps = 4);
    plt.show()
    """


################################
#### Animation sub-routine  ####
################################
def animate(nframe):
    print str(nframe) + '/' + str(len(zi))

    # Selects the data in the redshift slice
    mask = np.abs(allcone.field('Z_APP') - zi[nframe]) < dz_plot
    conedz = allcone[mask]

    lats = conedz.field('Dec')
    lons = conedz.field('RA')

    # Selects the data selected in the previous selection
    lats_3colors = np.array([])
    lons_3colors = np.array([])
    global common_GALID
    common_GALID = set(conedz.field(galid)) & set(list_GALID)
    for ids in common_GALID:
        lats_3colors = np.append(lats_3colors, conedz[np.where(conedz.field(galid) == ids)].field('Dec'))
        lons_3colors = np.append(lons_3colors, conedz[np.where(conedz.field(galid) == ids)].field('RA'))

    plt.cla()
    m = Basemap(projection='merc', lon_0=0, lat_0=0, llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat,
                celestial=True)

    # Lattitudes and longtitudes
    poslines = [-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]
    m.drawparallels(poslines, labels=[1, 0, 0, 0])
    m.drawmeridians(poslines, labels=[0, 0, 0, 1])

    # draw points
    x, y = m(lons, lats)
    m.scatter(x, y, 0.03, marker='o', color='b')
    x_3colors, y_3colors = m(lons_3colors, lats_3colors)
    m.scatter(x_3colors, y_3colors, 10, marker='o', color='r')

    # Adds a title
    plt.title('z=' + str(zi[nframe]))


#############################
#### Gaussian Selection  ####
#############################
def selec_gauss():
    print "##################################################"
    print "#######      Selection method 2:                 #"
    print "#######  just selecting statistically in z bins  #"
    print "##################################################"

    global list_GALID, cone_selection, selection

    list_GALID = []
    selection = [] # Table of all data for selected objects (all redshift samples)

    for i in np.arange(len(selection_properties)):

        print "Redshift: ~" + str(selection_properties['z'][i]) + "; Limit Magnitude: " + str(selection_properties['LimitMag'][i])
        mask_z = np.abs(allcone.field('Z_APP') - selection_properties['z'][i]) < dz
        mask_mag = allcone.field(selection_properties['Filter3'][i]) < selection_properties['LimitMag'][i]
        mask = mask_mag & mask_z
        print strftime("%Y-%m-%d %H:%M:%S", gmtime())

        # Redshift vs Magnitude in filter3
        cone_z = allcone[mask_z]
        fig = plt.figure()
        plt.title("Redshift vs "+selection_properties['Filter3'][i])
        plt.xlabel("Apparent Redshift (Z_APP)")
        plt.ylabel(selection_properties['Filter3'][i])
        plt.hist2d(cone_z['Z_APP'], cone_z[selection_properties['Filter3'][i]], bins=100, range=[[min(cone_z['Z_APP']),max(cone_z['Z_APP'])], [23, 35]])
        plt.plot([min(cone_z['Z_APP']), max(cone_z['Z_APP'])], [selection_properties['LimitMag'][i], selection_properties['LimitMag'][i]])
        #plt.ylim(24, 35)
        #plt.show()
        savemyplot(fig, "z_vs_" + str(selection_properties['Filter3'][i])+ "_" + str(selection_properties['z'][i]))
        plt.close()

        # Magnitude in filter3 distribution for a certain z bin
        fig = plt.figure()
        plt.title(selection_properties['Filter3'][i]+" for z~"+str(selection_properties['z'][i]))
        plt.xlabel(selection_properties['Filter3'][i])
        plt.ylabel("#")
        plt.hist(cone_z[selection_properties['Filter3'][i]], bins=100, range=[23, 35])
        plt.plot([selection_properties['LimitMag'][i], selection_properties['LimitMag'][i]], [0, 10000])
        #plt.ylim(24, 35)
        #plt.show()
        savemyplot(fig, "Filter "+str(selection_properties['Filter3'][i])+ "_for_z_" + str(selection_properties['z'][i]))
        plt.close()




        # selecting all objects with z>2 takes 20 min
        cone = allcone[mask]

        # I want to select a gaussian distribution from a random distribution, selecting all the objects at the peak of the gaussian.

        # Paremeters of the gaussian distribution:
        distribution_parameters = [1., selection_properties['z'][i], dz / 2.5]

        nbins = 10

        # Current distribution:
        hist, bin_edges = np.histogram(cone['Z_APP'], bins=nbins, range=(selection_properties['z'][i] - dz, selection_properties['z'][i] + dz))

        # Probability of selection in order to inverse the distribution to a gaussian distribution:
        Pz = gauss(bin_edges[:-1],
                   *distribution_parameters) / hist  # This can lead to 0/0, but that's a minor point, considering it will result in an expected 0.
        Pz = Pz / np.min(Pz[0.4 * len(Pz):0.6 * len(Pz)])  # The min does a light over selection.

        # Objects are selected randomly depending on the probability Pz of the bin they belong to.
        mask_gaussian = np.array([], dtype=bool)
        for objects in cone:
            proba = Pz[(np.where(objects.field('Z_APP') <= bin_edges))[0][0] - 1]
            mask_gaussian = np.append(mask_gaussian, [proba > random()])

        # gaussian distributed selection:
        cone_gaussian = cone[mask_gaussian]

        gaussian_selection['# objects in z bin'][i] = len(cone)
        gaussian_selection['# objects gaussian'][i] = len(cone_gaussian)
        gaussian_selection['# objects PFS'][i] = int(round(len(cone_gaussian) * Ratio_FoV))

        print "Total number of elements in the redshift bin : " + str(gaussian_selection['# objects in z bin'][i])
        print "Number of elements after gaussian selection : " + str(gaussian_selection['# objects gaussian'][i])
        print "Number of elements after FoV correction : " + str(gaussian_selection['# objects PFS'][i])

        for ids in cone_gaussian.field(galid):
            list_GALID.append(ids)

        selection.append(Table(cone_gaussian))

        fig = plt.figure()
        plt.title("Redshift distribution for the gaussian selection @ z~" + str(selection_properties['z'][i]) + "\n Objects selected only inside PFS FoV: " + str(gaussian_selection['# objects PFS'][i]))
        plt.xlabel("Apparent Redshift (Z_APP)")
        plt.ylabel("#")
        plt.hist(cone['Z_APP'], bins=nbins,label="Initial distribution: " + str(gaussian_selection['# objects in z bin'][i]), range=(selection_properties['z'][i] - dz, selection_properties['z'][i] + dz))
        plt.hist(cone_gaussian['Z_APP'], bins=nbins, label="Gaussian selection: " + str(gaussian_selection['# objects gaussian'][i]), range=(selection_properties['z'][i] - dz, selection_properties['z'][i] + dz))
        #plt.plot(bin_edges[:-1], Pz)
        plt.plot(bin_edges, gauss(bin_edges, *[max(np.histogram(cone_gaussian['Z_APP'], bins=nbins, range=(selection_properties['z'][i] - dz, selection_properties['z'][i] + dz))[0]), distribution_parameters[1], distribution_parameters[2]]), label="just a gaussian")
        #plt.plot(bin_edges[:-1], hist_densities)
        #plt.plot(bin_edges[:-1], gauss(bin_edges[:-1], *distribution_parameters))
        plt.legend()
        plt.xlim(selection_properties['z'][i] - dz, selection_properties['z'][i] + dz)
        savemyplot(fig, "z_dist_gaussian_selection_z_" + str(selection_properties['z'][i]))
        #	plt.show()
        plt.close()

    selection = vstack(selection)
    #print min(selection.field('Z_APP'))
    #print max(selection.field('Z_APP'))


    print gaussian_selection
    ascii.write(gaussian_selection, plot_directory+'gaussian_selection.txt', format=table_write_format)

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


#############################
#### Dropouts selection  ####
#############################
def selec_3colors():
    print "##################################################"
    print "#######      Selection method 1:                 #"
    print "#######    real 3colors selction                 #"
    print "##################################################"

    global conelist, cone, list_GALID, allcone_selected_3colors, selection

    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
    allcone_selected_3colors = np.zeros(len(allcone), dtype=int)
    print strftime("%Y-%m-%d %H:%M:%S", gmtime())

    number_duplicates = 0

    bins = [1, 2, 3, 4, 5]
    #bins = [1]
    #bins = np.arange(len(selection_properties))

    conelist = []
    list_GALID = []
    selection = [] # Table of all data for selected objects (all redshift samples)


    for i in bins:
        print "redshift: ~" + str(selection_properties['z'][i]) + ". Filters : " + str(
            selection_properties['Filter1'][i]) + " " + str(selection_properties['Filter2'][i]) + " " + str(
            selection_properties['Filter3'][i])
        #mask_z = np.abs( allcone.field('Z_APP') - selection_properties['z'][i] ) < dz
        mask_mag = allcone.field(selection_properties['Filter3'][i]) < selection_properties['LimitMag'][i]

        #print "size mask_mag : "+ str(len(mask_mag))
        #print "number True :"+ str(np.count_nonzero(mask_mag))

        cone = allcone[mask_mag]
        color_selection['# objects under mag lim'][i] = len(cone)

        print "Number of candidates with mag[" + str(selection_properties['Filter3'][i]) + "]>" + str(
            selection_properties['LimitMag'][i]) + ": " + str(color_selection['# objects under mag lim'][i])

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
        plt.title("ZOOM: Number of particles for " + str(selection_properties['Filter3'][i]) + " < " + str(
            selection_properties['LimitMag'][i]))
        plt.xlabel("NP")
        plt.ylabel("#")
        #plt.yscale('log')
        plt.hist(cone['NP'], bins=1000, range=(0, 1000))
        #plt.show()
        savemyplot(fig, "NP_ZOOM_filter3_le_" + str(selection_properties['LimitMag'][i]))
        plt.close()

        fig = plt.figure()
        plt.title("Number of particles for " + str(selection_properties['Filter3'][i]) + " < " + str(
            selection_properties['LimitMag'][i]))
        plt.xlabel("NP")
        plt.ylabel("#")
        #plt.yscale('log')
        plt.hist(cone['NP'], bins=1000)
        #plt.show()
        savemyplot(fig, "NP_filter3_le_" + str(selection_properties['LimitMag'][i]))
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
        m = (limitH - endV) / (endH - limitV)
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
        plt.title("Colors for z~" + str(selection_properties['z'][i]))
        plt.xlabel(selection_properties['Filter2'][i] + "-" + selection_properties['Filter3'][i])
        plt.ylabel(selection_properties['Filter1'][i] + "-" + selection_properties['Filter2'][i])

        # histogram
        plt.hist2d(cone[selection_properties['Filter2'][i]] - cone[selection_properties['Filter3'][i]],
                   cone[selection_properties['Filter1'][i]] - cone[selection_properties['Filter2'][i]], bins=150,
                   range=([-1., 2.5], [-1., 8.5]))
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
        duplicates = set(list_GALID) & set(cone.field(galid))
        print "Number of duplicates: " + str(len(duplicates))
        if len(duplicates) > 0:
            number_duplicates = number_duplicates + len(duplicates)
            mask_duplicates = np.ones(len(cone), dtype=bool)
            for dupie in duplicates:
                mask_duplicates[np.where(cone.field(galid) == dupie)] = False
            cone = cone[mask_duplicates]
            print "After deleting the duplicates in the new cone, number of objects in the cone: " + str(len(cone))

        color_selection['# objects color selected'][i] = len(cone)
        print "Number of galaxies selected by color : " + str(color_selection['# objects color selected'][i])

        for ids in cone.field(galid):
            list_GALID.append(ids)

        selection.append(Table(cone))

        # Adding this cone bin to a list of all "cones" at different redshifts
        conelist.append(cone)


        # plotting individual points for selected objects
        plt.scatter(cone[selection_properties['Filter2'][i]] - cone[selection_properties['Filter3'][i]],
                    cone[selection_properties['Filter1'][i]] - cone[selection_properties['Filter2'][i]],
                    c=cone['Z_APP'], vmin=selection_properties['z'][i] - 1, vmax=selection_properties['z'][i] + 1,
                    cmap=plt.cm.spectral)
        cbar = plt.colorbar()
        cbar.set_label('redshift')

        # plotting the limits
        plt.plot([-1., endH], [limitH, limitH], '-b')
        plt.plot([limitV, limitV], [endV, 8.5], '-b')
        plt.plot([endH, limitV], [limitH, endV], '-b')

        if selection_properties['z'][i] == 3:
            plt.xlim(-1., 2.5)
            plt.ylim(-1., 8.5)
        else:
            plt.xlim(-0.5, 1.5)
            plt.ylim(-1., 4.2)


        #plt.show()
        savemyplot(fig, "Colors_z_" + str(selection_properties['z'][i]))
        plt.close()

        fig = plt.figure()
        plt.title("ZOOM after 3 color selection: Number of particles @ z~" + str(selection_properties['z'][i]))
        plt.xlabel("NP")
        plt.ylabel("#")
        #plt.yscale('log')
        plt.hist(cone['NP'], bins=1000, range=(0, 1000))
        #plt.show()
        savemyplot(fig, "NP_ZOOM_after_3c_selection_filter3_z_" + str(selection_properties['z'][i]))
        plt.close()

        fig = plt.figure()
        plt.title("After 3 color selection: Number of particles @ z~" + str(selection_properties['z'][i]))
        plt.xlabel("NP")
        plt.ylabel("#")
        #plt.yscale('log')
        plt.hist(cone['NP'])
        #plt.show()
        savemyplot(fig, "NP_after_3c_selection_filter3_z_" + str(selection_properties['z'][i]))
        plt.close()

        fig = plt.figure()
        plt.title("Redshift distribution for the 3 color selection @ z~" + str(selection_properties['z'][i]))
        plt.xlabel("Apparent Redshift (Z_APP)")
        plt.ylabel("#")
        plt.hist(cone['Z_APP'], bins=40)
        #plt.show()
        savemyplot(fig, "z_dist_color_selection_z_" + str(selection_properties['z'][i]))
        plt.close()

        color_selection['# objects PFS'][i] = int(round(color_selection['# objects color selected'][i] * Ratio_FoV))

    fig = plt.figure()
    plt.title("Redshift distribution for the 3 color selection")
    plt.xlabel("Apparent Redshift (Z_APP)")
    plt.ylabel("#")
    for conei in conelist:
        plt.hist(conei['Z_APP'], bins=40, range=(0, 9), histtype='step',
                 label="$z_{median} \sim" + str("%0.1f" % (np.median(conei['Z_APP']))) + "$ " + str(
                     "%6.f" % len(conei)) + " objects")
    plt.legend()
    #plt.show()
    savemyplot(fig, "z_dist_color_selection")
    plt.close()

    print color_selection
    ascii.write(color_selection, plot_directory+'color_selection.txt', format=table_write_format)

    if number_duplicates != 0:
        print "There was " + str(number_duplicates) + " duplicates in the selection. They have been taken care of."

    selection = vstack(selection)


#############################
####  A simple test for  ####
####   numeration and    ####
####   compare with R    ####
#############################
def test_z2():
    mask_z = np.abs(allcone.field('Z_APP') - 2.) < .5
    mask_mag = allcone.field('SDSS_I') < 24.
    mask = mask_mag & mask_z
    # Redshift vs Magnitude in filter3
    cone = allcone[mask]
    print len(cone)




###########################
#### Opens CFHTLS cats ####
###########################
def open_CFHTLS(field_number):
    global CFHTLS

    # Table of the CFHTLS fields
    field = [1, 2, 3, 4]
    field_name = ["D1", "D2", "D3", "D4"]
    pointing = ["022559-042940", "100028+021230", "141927+524056", "221531-174356"]
    CFHTLS_fields = Table([field, field_name, pointing], names=('field', 'field name', 'pointing'),
                          meta={'name': 'table of the CFHTLS fields'})

    # Catalog path
    catpath = "./data/CFHTLS/"
    catname = "CFHTLS_D-85_ugriyz_" + CFHTLS_fields['pointing'][field_number] + "_T0007_SIGWEI_MAGAUTO.cat"

    CFHTLS = Table.read(catpath + catname, format='ascii', header_start=15)
    #CFHTLS = Table.read(catpath+catname, format='ascii', data_end=10000, header_start=15)

    print "There are " + str(len(CFHTLS)) + " objects in the CFHTLS catalog: " + catname

    print CFHTLS


########################################
#### Dropouts selection  FOR CFHTLS ####
########################################
def selec_3colors_CFHTLS():
    print "##################################################"
    print "#######    CFHTLS dropout selction               #"
    print "##################################################"

    bins = [1, 2, 3]
    #bins = [1]
    #bins = np.arange(len(selection_properties))

    #conelist = []
    #list_GALID = []

    for i in bins:
        print "redshift: ~" + str(selection_properties['z'][i]) + ". Filters : " + str(
            selection_properties['Filter1_CFHTLS'][i]) + " " + str(
            selection_properties['Filter2_CFHTLS'][i]) + " " + str(selection_properties['Filter3_CFHTLS'][i])
        #mask_z = np.abs( allcone.field('Z_APP') - selection_properties['z'][i] ) < dz
        print selection_properties['Filter3_CFHTLS'][i]
        print selection_properties['LimitMag'][i]
        print CFHTLS[selection_properties['Filter3_CFHTLS'][i]]

        mask_mag = CFHTLS[selection_properties['Filter3_CFHTLS'][i]] < selection_properties['LimitMag'][i]

        #print "size mask_mag : "+ str(len(mask_mag))
        #print "number True :"+ str(np.count_nonzero(mask_mag))

        sub_CFHTLS = CFHTLS[mask_mag]
        color_selection_CFHTLS['# objects under mag lim CFHTLS'][i] = len(sub_CFHTLS)

        print "Number of candidates with mag[" + str(selection_properties['Filter3_CFHTLS'][i]) + "]>" + str(
            selection_properties['LimitMag'][i]) + ": " + str(len(sub_CFHTLS))

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


        ###########################################
        #    3 color selection is done here:      #
        ###########################################

        limitH = selection_properties['selec: limitH'][i]
        limitV = selection_properties['selec: limitV'][i]
        endH = selection_properties['selec: endH'][i]
        endV = selection_properties['selec: endV'][i]

        # y = m x + p
        m = (limitH - endV) / (endH - limitV)
        p = limitH - m * endH

        f1minusf2 = sub_CFHTLS.field(selection_properties['Filter1_CFHTLS'][i]) - sub_CFHTLS.field(
            selection_properties['Filter2_CFHTLS'][i])
        f2minusf3 = sub_CFHTLS.field(selection_properties['Filter2_CFHTLS'][i]) - sub_CFHTLS.field(
            selection_properties['Filter3_CFHTLS'][i])

        mask_colorX = f2minusf3 < limitV
        mask_colorY = f1minusf2 > limitH
        mask_color_mp = f1minusf2 > m * f2minusf3 + p
        mask = mask_colorX & mask_colorY & mask_color_mp

        #print "size mask : "+ str(len(mask))

        # Color diagram : histogram + individual points for selected objects
        fig = plt.figure()
        plt.title("Colors for z~" + str(selection_properties['z'][i]))
        plt.xlabel(selection_properties['Filter2_CFHTLS'][i] + "-" + selection_properties['Filter3_CFHTLS'][i])
        plt.ylabel(selection_properties['Filter1_CFHTLS'][i] + "-" + selection_properties['Filter2_CFHTLS'][i])

        # histogram
        plt.hist2d(sub_CFHTLS[selection_properties['Filter2_CFHTLS'][i]] - sub_CFHTLS[
            selection_properties['Filter3_CFHTLS'][i]],
                   sub_CFHTLS[selection_properties['Filter1_CFHTLS'][i]] - sub_CFHTLS[
                       selection_properties['Filter2_CFHTLS'][i]], bins=150, range=([-1., 2.5], [-1., 8.5]),
                   norm=LogNorm())
        #plt.plot(sub_CFHTLS[selection_properties['Filter2_CFHTLS'][i]] - sub_CFHTLS[selection_properties['Filter3_CFHTLS'][i]], sub_CFHTLS[selection_properties['Filter1_CFHTLS'][i]] - sub_CFHTLS[selection_properties['Filter2_CFHTLS'][i]], '.')

        # selecting objects
        sub_CFHTLS = sub_CFHTLS[mask]
        print "Number of objects after 3 colors selection: " + str(len(sub_CFHTLS))

        # Marks selected objects in the initial cone:
        #allcone_selected_3colors[mask_mag] = 1
        #allcone_selected_3colors[np.where(allcone_selected_3colors == 1)] = 1
        #print len(np.where(allcone_selected_3colors == True))


        color_selection_CFHTLS['# objects color selected CFHTLS'][i] = len(sub_CFHTLS)
        print "Number of galaxies selected by color : " + str(
            color_selection_CFHTLS['# objects color selected CFHTLS'][i])



        # plotting individual points for selected objects
        #plt.scatter(sub_CFHTLS[selection_properties['Filter2_CFHTLS'][i]] - sub_CFHTLS[selection_properties['Filter3_CFHTLS'][i]], sub_CFHTLS[selection_properties['Filter1_CFHTLS'][i]] - sub_CFHTLS[selection_properties['Filter2_CFHTLS'][i]])

        # plotting the limits
        plt.plot([-1., endH], [limitH, limitH], '-b')
        plt.plot([limitV, limitV], [endV, 8.5], '-b')
        plt.plot([endH, limitV], [limitH, endV], '-b')

        if selection_properties['z'][i] == 3:
            plt.xlim(-1., 2.5)
            plt.ylim(-1., 8.5)
        else:
            plt.xlim(-0.5, 1.5)
            plt.ylim(-1., 4.2)

        plt.show()
        savemyplot(fig, "CFHTLS_" + str(selection_properties['z'][i]))
        plt.close()

        color_selection_CFHTLS['# objects PFS CFHTLS'][i] = int(
            round(color_selection_CFHTLS['# objects color selected CFHTLS'][i] * Ratio_FoV_PFS_CFHTLS))

    print color_selection_CFHTLS
    ascii.write(color_selection_CFHTLS, plot_directory+'color_selection_CFHTLS.txt', format=table_write_format)




if __name__ == '__main__':
    main()

