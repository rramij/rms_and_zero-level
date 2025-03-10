# To calculate local RMS and zero level around the phase center
# Written by Ramij Raja
# Date: 20 June 2023 # Last Update: 15 Nov 2023
# RUN in Python3
# conda3
#########################################
import sys
import numpy as np
from astropy.io import fits
from astropy.io.fits import getdata
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
#############################
# Main Code 
#############################
fits_image=sys.argv[1]
data, hdr = getdata(fits_image, header=True)

pc = int((hdr['CRPIX1']+hdr['CRPIX2'])/2) - 2 #Phase center pixel
box_size=1500 #arcsec
print('Calculating image stats around the Phase Center within 1500 arcsec box')
cell = abs(hdr['CDELT1']*3600) #arcsec
box_sizeP=int(round(box_size/cell)) #pixels
# Region limits
box_lim = {'xl1':int(round(pc-box_sizeP/2)), 'yl1':int(round(pc-box_sizeP/2)), 'xl2':int(round(pc+box_sizeP/2)), 'yl2':int(round(pc+box_sizeP/2))}

# Local data
class RadioError(Exception):
    """Base class for exceptions in this module."""
    pass

nax = hdr['NAXIS']
if nax == 4:
  a = data[0,0,box_lim['yl1']:box_lim['yl2'],box_lim['xl1']:box_lim['xl2']]
elif nax == 2:
  a = data[box_lim['yl1']:box_lim['yl2'],box_lim['xl1']:box_lim['xl2']]
else:
  raise RadioError('Too many or too few axes to proceed (%i)' % nax)

# Stats before sigma clipping
img_min = np.nanmin(a)
img_max = np.nanmax(a)
img_mean = np.nanmean(a)
img_median = np.nanmedian(a)
img_std = np.nanstd(a)
print('Image stats before sigma clipping')
print('  Min	    Max	    Mean     Median   std')
print('{:.3g}'.format(img_min),'{:.3g}'.format(img_max),'{:.3g}'.format(img_mean),'{:.3g}'.format(img_median),'{:.3g}'.format(img_std))
'''
sigclip=SigmaClip(sigma=2,maxiters=None)
clipD = sigclip(a)
imin = np.nanmin(clipD)
imax = np.nanmax(clipD)
imean = np.nanmean(clipD)
imedian = np.nanmedian(clipD)
istd = np.std(clipD)
print('Image stats after sigma clipping')
print('Min		Max		Mean		Median		std')
print(imin,imax,imean,imedian,istd)
'''
image_mean, image_median, image_stddev = sigma_clipped_stats(a,sigma=2.7,maxiters=None)
print("------------------------------")
print('Image stats after sigma clipping at 2.7 sigma around the median')
print('  Min	    Max	    Mean     Median   std')
print('{:.3g}'.format(image_median-(3*image_stddev)),'{:.3g}'.format(image_median+(3*image_stddev)),'{:.3g}'.format(image_mean),'{:.3g}'.format(image_median),'{:.3g}'.format(image_stddev))
print('-------------------------------')
print('Local RMS: ', '{:.3g}'.format(image_stddev), ' Jy/beam')
print('-------------------------------')
print('Zero level:', '{:.3g}'.format((image_mean+image_median)/2), ' Jy/beam')
print('-------------------------------')
