## To calculate RMS noise of an FITS image, run : python rms_calculation.py FITSfile
## The python should be python2
############################################
import sys

import glob

import astropy

from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt

from astropy.stats import SigmaClip

from astropy.stats import sigma_clipped_stats

from astropy.wcs import WCS

from astropy.io import fits

from astropy.utils.data import get_pkg_data_filename
####################

file1=sys.argv[1]

image = fits.open(file1)

prihdr = image[0].header

data = image[0].data

image_data = data#[0,0,:,:]

image_data=np.nan_to_num(image_data)



img_min, img_max = np.min(image_data),np.max(image_data)

img_mean=np.mean(image_data);img_median=np.median(image_data);img_std=np.std(image_data)

print img_min,img_max,img_mean,img_median,img_std



sigclip = SigmaClip(sigma = 3.0, iters=5)

image_data_sigclip = sigclip(image_data)

image_min_clip, image_max_clip = np.min(image_data_sigclip),np.max(image_data_sigclip)

image_mean_clip=np.mean(image_data_sigclip);image_median_clip=np.median(image_data_sigclip);image_std_clip=np.std(image_data_sigclip)

print image_min_clip,image_max_clip,image_mean_clip,image_median_clip,image_std_clip





image_mean, image_median, image_stddev = sigma_clipped_stats(image_data, sigma = 3.0, iters = 5)

print image_mean,image_median

print "The calculated rms using sigma_clipped_stats(image_data, sigma = 3.0, iters = 5) is ="

print image_stddev
