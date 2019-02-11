"""
Title:          Fractal Dimension Calculation on Polaris FITS data
Author:         James Beattie
Creation Date:  12th Feb, 2019
"""

# Imports
###########################################################################

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from skimage import measure, filters, morphology   # for drawing contours and Gaussian filters
from matplotlib import rc                          # nicer text in matplotlib
from skimage.morphology import convex_hull_image

###########################################################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


# Directories
###########################################################################

dataAdd     = "./FITSfiles/"
dataWrite   = "./PolarisProcessing"

# Function Definitions
###########################################################################

def FITSreader(file,dataAddress):
    """
    This reads in FITS files and immediatly extracts the data from themself.

    INPUTS:
    file            - the file name
    dataAddress     - the address of the file
    """

    data    = fits.open(dataAddress+file)
    data    = data[0].data

    return data


def MaxPixelWindow(image):
    max_pix_x   = []
    max_pix_y   = []
    start       = 0
    dy          = 5;
    top         = image.shape[0]/dy
    BrickTop    = 70
    BrickBottom = 580

    for end in xrange(1,top):
        if dy*end > BrickTop and dy*end < BrickBottom:
            y,x = np.where(image[start:dy*end,:] == image[start:dy*end,:].max())
            max_pix_y.append(y+start)
            max_pix_x.append(x)
            start += dy
        else:
            start += dy
            continue

    return np.array([max_pix_x, max_pix_y])

# Coordinates
###########################################################################

# # some info on coordinates:
# coord_x_pix     = Brick[0].header["CRPIX1"]
# coord_y_pix     = Brick[0].header["CRPIX2"]
# coord_x_offset  = Brick[0].header["CRVAL1"]
# coord_y_offset  = Brick[0].header["CRVAL2"]
#
# # compute the size of a pixel
# pix = Brick[0].header["CDELT2"]*3600.0                                      # in arcsec
# print Brick[0].header["CDELT1"]*3600.0, Brick[0].header["CDELT2"]*3600.0    # in arcsec
# print 'pixel size in arcsec = ', pix
#
# distance    = np.array([8.3,8.0,8.6])   # kpc
# arcsec2rad  = 206265.0                  # arcseconds to radians
# pix_pc = distance * 1000.0 * pix/arcsec2rad
# print 'pix_pc = ', pix_pc


# Column Density
###########################################################################

# sensitivity     = 25e-3*1.9e23
# outline1        = 5e22
# outline2        = 1e23


# Polaris_image   = np.flipud(Brickcf[0].data)   # Federrath 2016 preprocessed
# Brick_image     = Brick[0].data[0,0,:,:]        # additional indexing needed for the first fits file
#
# # Column Density Indexes (all non-zero values)
# nz_index    = np.where(Brickcf_image.ravel() >= 2*sensitivity)  # indexes of densites
# sigma       = Brickcf_image.ravel()[nz_index];      # column density array cm^-2
# sigma_mean  = np.mean(sigma)                        # mean column density
# sigma_std   = np.std(sigma)                         # standard deviation of column density

# Contours
###########################################################################

#plt.hist(np.log10(sigma),bins=100)
#plt.axvline(np.log10(sigma_mean), color='r')
#plt.show()

# countours_cf = measure.find_contours( Brickcf_image, outline1)
# countours_JB = measure.find_contours( Brickcf_image, outline2)

# Image Processing
###########################################################################

# mask    = morphology.binary_closing(Brickcf_image >= outline1)
# mask    = morphology.binary_closing(mask)



# max_pix = MaxPixelWindow(image=Brickcf_image)

# Working Script
###########################################################################

Pol_qui     = FITSreader(,dataAdd)
Pol_sax     = FITSreader(,dataAdd)



# Plots
###########################################################################

# pix_1pc     = 1/pix_pc[0]
# dx          = 10
# dy          = 600
# scale_bar   = np.array([[dx,dx+pix_1pc],[dy,dy]])
# fs          = 18

# fig, ax = plt.subplots()
# for contour in countours_cf:
#     ax.plot(contour[:, 1], contour[:, 0], linewidth=2, color='black')
# for contour in countours_JB:
#     plt.plot(contour[:, 1], contour[:, 0], linewidth=2, color='black')
# plt.plot(max_pix[0],max_pix[1],'ko',linewidth=2)
# plt.imshow( np.log10( Brickcf_image*(Brickcf_image >= outline1) ), vmin=22.5,vmax=23.5, cmap=plt.cm.seismic)
# plt.annotate(r'ALMA + $Herschel$',xy=(10, 50),fontsize=fs,color='black')
# plt.annotate(r'G0.253+0.016',xy=(10, 65),fontsize=fs,color='black')
# plt.annotate(r'1pc',xy=(45, 595),fontsize=fs-2,color='black')
# plt.plot(scale_bar[0],scale_bar[1],color='black',linewidth=2)
# ax.set_xticks([])
# ax.set_yticks([])
# cbar = plt.colorbar()
# cbar.set_label(r"$\log_{10} \Sigma$ [cm$^{-2}$]",fontsize=fs)
# plt.show()
