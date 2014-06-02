import numpy
import astropy
from astropy.io import ascii
import pyfits
import math
from time import gmtime, strftime
import sys
from sys import exit,argv
import csv
import os
from string import upper,lower
from numpy import recfromcsv
import matplotlib.pyplot as plt
from numpy.random import normal

def main():

#  Opening Millennium Run lightcone catalog in fits format 
  path = './data/'
  conename = path + 'wmap1_bc03_001_igm1.fits'
  hdulist = pyfits.open(conename)
  cone = hdulist[1].data
  cols = hdulist[1].columns
  hdulist.close()

  #some examples using it
  
  #print the number of objects in the cone
  print len(cone)

  #see the column names:
  print cols

  #print the values of column named 'Z_APP' (redshift):
  print cone['Z_APP']
  
  
  

  #print the highest star formation rate in the catalog:
  print max(cone['SFR'])

  #now make some plots with matplotlib!
  
if __name__ == '__main__':
  main()
