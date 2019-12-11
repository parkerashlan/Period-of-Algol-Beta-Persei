# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 16:08:02 2019

@author: parke
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as spis
import scipy.interpolate as spin
import os

path = './Algolfits'

algolfluxtxt = open('Algolfits.txt','w+')
algolfluxtxt.write('AlgolFlux'+' '+'JulianDate'+'\n')
algolfluxtxt.close()

algolfluxdata = []
JD = []
Nr = []
df = []

for filename in os.listdir(path):
    varstar = fits.open(path+'/'+filename)
    
    """Loads image data into image variable and loads julian date into date
    variable"""
    image = varstar[0].data 
    date = varstar[0].header['MJD-OBS'] 
    Nr.append(varstar[0].header['rdnoise'])
    name = varstar[0].header['OBJECT']
    varstar.close()

    """Finds star by looking at range of values which would only be within the object of interest"""
    starcenter = np.where(image <= 3500) and np.where(image >= 2000)
    
    """Uses the middle of the array to use as the center of the star because the array is sorted and 
    the star is near the middle of the image"""
    middlex = int(len(starcenter[0])/2)
    middley = int(len(starcenter[1])/2)

    """Creates a box around the star 30 pixels wide and 30 pixels high"""
    starx = [starcenter[0][middlex]-15,starcenter[0][middlex]+15]
    stary = [starcenter[1][middley]-15,starcenter[1][middley]+15]

    """Slicing the image around the star"""
    starbox = image[starx[0]:starx[1],stary[0]:stary[1]]
    
    """Finds the indexes of the star pixels by looking at all values above 1000
    the same with the background pixels"""
    starloc = np.where(starbox > 1200)
    backgroundloc = np.where(starbox <= 1200)
    
    """Uses the indexes for the star location to find the pixels which correlate to the star"""
    algolstar = starbox[starloc[0],starloc[1]]
    background = starbox[backgroundloc[0],backgroundloc[1]]

    """Takes average of the star pixels and the background pixels and subtracts them
    to find the flux of the star"""
    algolg = np.average(algolstar)
    backgroundg = np.average(background)
    algolflux = algolg - backgroundg
    print(algolflux)
    
    """Find the peaks of brightness along a line through the middle of the
    starbox then calculates the FWHM which is the inner aperture radius"""
    aperture = starbox[int(len(starbox)/2),:]
    peaks,_ = spis.find_peaks(aperture)
    half = spis.peak_widths(aperture, peaks, 0.5)
    fwhm = np.max(half[0])

    """Error Analysis"""
    Nia = algolstar.size
    Noa = background.size
    dF = np.sqrt(np.sum((- np.square(14.5)+algolstar))-((np.square((Nia/Noa)))*(np.sum(np.square(14.5)+background))))
    df.append(dF)
    print('error',dF)

    """Graph various data to check it"""
    plt.figure()
#    x = np.linspace(0,30,30)
#    plt.plot(x,aperture)
    plt.imshow(image,cmap='gray')
#    plt.plot(peaks,aperture[peaks],'x')
    plt.show()
#    print(name)
    
    """Appends data to lists to be plotted later"""
    algolfluxdata.append(algolflux)
    JD.append(date)
    
    """Writes the data to a text file"""
    algolfluxtxt = open('Algolfits.txt','a')
    algolfluxtxt.write(str(algolflux)+' '+str(date)+'\n')
    algolfluxtxt.close()

"""Plots data points along with error bars"""
plt.figure()
plt.xlabel('Time (MJD)')
plt.ylabel('Flux (e-)')
plt.title('Algol Flux')
plt.plot(JD,algolfluxdata,'ro')
plt.errorbar(JD,algolfluxdata,yerr=df,fmt='none',ecolor='r',capsize=3)

"""Fits polynomial line to the data"""
plt.figure()
plt.title('Extrapolated Light Curve of Algol')
x_new = np.linspace(JD[0],JD[-1])
f = spin.interp1d(JD,algolfluxdata,kind='cubic', fill_value='extrapolate')
y_new = f(x_new)
plt.errorbar(x_new,y_new)
plt.show()
