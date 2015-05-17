__author__ = 'loic'

""" Plots the figure needed for the ssp draft
with Jenny and Roderik. """

import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import matplotlib.cm as cm
from astropy.io import fits
import sys
from matplotlib import rcParams
import os

selection = "COLSEL3"
#selection = "gaussian"
dir = "plots/planck1_m05_002_igm1_nov7/"+selection+"/"
catname = selection+"_selection_with_densities"
dirplus = dir+"proposal/"

if not os.path.exists(dirplus):
    os.makedirs(dirplus)


"""
# Creates a numpy file from the table in ascii
# (something I should have done in the original code)
# Actually just saves it in fits for Jenny
create_numpy_file = True
if create_numpy_file:
    cat = Table.read(dir+catname+".txt", format="ascii", delimiter='|')
    #print cat
    #np.save(dir+catname+".npy", cat)


#cat = np.load(dir+catname+".npy")
#print cat
"""


cat = Table.read(dir+catname+".txt", format="ascii", delimiter='|')
cat.remove_column('col0')
cat.remove_column('col54')

"""
cat.write(dir+catname+'.fits', overwrite=True)

hdulist = fits.open(dir+catname+'.fits')
print hdulist[1].columns

#print (fits.open("data/lightcones/planck1_m05_002_igm1_nov7.fits"))[1].data

sys.exit()
"""





fig = plt.figure()
plt.title("CENTRALMVIR vs 9th NN distance")
plt.hist2d(np.log10(cat["CENTRALMVIR"]), np.log10(cat["Dist_nearest_9_in_Mpch"]), bins=100)
plt.xlabel("CENTRALMVIR")
plt.ylabel("Dist_nearest_9_in_Mpch")
#plt.show()
plt.savefig(dirplus+"CENTRALMVIRvs9.png")
plt.close()

fig = plt.figure()
plt.title("MVIR vs 9th NN distance")
plt.hist2d(np.log10(cat["MVIR"]), np.log10(cat["Dist_nearest_9_in_Mpch"]), bins=100)
plt.xlabel("MVIR")
plt.ylabel("Dist_nearest_9_in_Mpch")
#plt.show()
plt.savefig(dirplus+"MVIRvs9.png")
plt.close()


fig = plt.figure()
plt.title("redshift distribution")
plt.hist(cat["Z_APP"])
plt.show()
plt.close()


cat = cat[cat["Z_APP"] >= 0.5]
cat = cat[cat["Z_APP"] <= 2.]



y = "CENTRALMVIR"
z = "SFR"
x = "Dist_nearest_9_in_Mpch"
logx = np.log10(cat[x])
nbins = 20.
minx = np.min(logx)
maxx = np.max(logx)

medx = np.array([])
medy = np.array([])
medz = np.array([])

bins = np.arange(minx, maxx, (maxx-minx)/nbins)
for i in np.arange(nbins-1.):
    mask = (logx >= bins[i]) & (logx < bins[i+1])
    subcat = cat[mask]
    medx = np.append(medx, np.median(subcat[x]))
    medy = np.append(medy, np.median(subcat[y]))
    medz = np.append(medz, np.median(subcat[z]))



fig = plt.figure()
plt.title("Whole sample")
plt.plot(np.log10(medx), np.log10(medy), '-')
plt.scatter(np.log10(medx), np.log10(medy), c=np.log10(medz), s=80) #, cmap=cm.colormap_name
plt.xlabel("log "+x)
plt.ylabel("log "+y)
cb = plt.colorbar()
cb.set_label("log "+z)
#plt.show()
#plt.savefig(dirplus+"nameithere.png")
plt.close()


###########
# Version with redshift separation:


y = "CENTRALMVIR"
z = "SFR"
x = "Dist_nearest_9_in_Mpch"

xtit = "9th Nearest Neighbor Distance [Mpc/h]"
ytit = "Central Mvir [M$_{\odot}$]"
ztit = "SFR [M$_{\odot}$yr$^{-1}$]"

logx = np.log10(cat[x])
nbins = 20.
minx = np.min(logx)
maxx = np.max(logx)


fig = plt.figure(figsize=(10,10))

rcParams.update({'font.size': 20})




plt.title("")
#cnorm = Normalize(vmin=np.min(), vmax=np.max())

#plt.hist2d(np.log10(cat[x]), np.log10(cat[y]), bins=100, cmap="Greys")

znbins = 3
minz = np.min(cat["Z_APP"])
maxz = np.max(cat["Z_APP"])
zbins = np.arange(minz, maxz+(maxz-minz)/znbins, (maxz-minz)/znbins)

colorSFRmin = 0.
colorSFRmax = 0.

for j in np.arange(znbins):
    zmask = (cat["Z_APP"] >= zbins[j]) & (cat["Z_APP"] < zbins[j+1])
    zcat = cat[zmask]

    zlabel = "z = "+str("{0:.1f}".format(np.median(cat["Z_APP"][zmask])))

    medx = np.array([])
    medy = np.array([])
    medz = np.array([])

    bins = np.arange(minx, maxx+(maxx-minx)/nbins, (maxx-minx)/nbins)
    for i in np.arange(nbins):
        logx = np.log10(zcat[x])
        mask = (logx >= bins[i]) & (logx < bins[i+1])
        subcat = zcat[mask]
        medx = np.append(medx, np.median(subcat[x]))
        medy = np.append(medy, np.median(subcat[y]))
        medz = np.append(medz, np.median(subcat[z]))

    thiscolorSFRmin = np.log10(np.min(medz[medz !=0]))

    print colorSFRmin
    colorSFRmin = np.min([thiscolorSFRmin, colorSFRmin])
    print colorSFRmin
    colorSFRmax = np.max([np.max(np.log10(medz)), colorSFRmax])
    #print colorSFRmax


    plt.plot(np.log10(medx), np.log10(medy), '-', label=zlabel, linewidth=2.0, zorder=1)
    plt.scatter(np.log10(medx), np.log10(medy), c=np.log10(medz), s=80, vmin=-0.7, vmax=colorSFRmax, zorder=2) #, cmap=cm.colormap_name

plt.legend()
#plt.xlabel("log "+x, fontsize = 20)
plt.xlabel("log$_{10}$ "+xtit)
plt.ylabel("log$_{10}$ "+ytit)
cb = plt.colorbar()
cb.set_label("log$_{10}$ "+ztit)
#plt.show()
plt.savefig(dirplus+x+"_vs_"+y+"_vs_"+z+"_3z.png")
plt.close()


