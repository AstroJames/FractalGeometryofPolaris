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
from matplotlib.gridspec import GridSpec

###########################################################################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


# Directories
###########################################################################

dataAdd     = "./FITSfiles/"
dataWrite   = "./PolarisProcessing"

# Function Definitions
###########################################################################

def ColDensityFITSextract(file,dataAddress):
    """
    This reads in FITS files and immediatly extracts the column density data from themself.

    INPUTS:
    file            - the file name
    dataAddress     - the address of the file
    """

    # Read in the FITS files and extract column density.
    data    = fits.open(dataAddress+file)
    data    = data[0].data[0,0,:,:]

    # Check if the data is actually just the column denisty
    if data.ndim != 2:
        print("WARNING: The dimensions of the FITS read is not 2.")

    # The first two moments
    mean    = data.mean()
    var     = data.var()

    return data, mean, var


def MaxPixelWindow(image,windowSize,windowTop,windowBottom):
    max_pix_x   = []
    max_pix_y   = []
    start       = 0
    dy          = 5;
    top         = image.shape[0]/dy

    for end in xrange(1,top):
        if dy*end > windowTop and dy*end < windowBottom:
            y,x = np.where(image[start:dy*end,:] == image[start:dy*end,:].max())
            max_pix_y.append(y+start)
            max_pix_x.append(x)
            start += dy
        else:
            start += dy
            continue

    return np.array([max_pix_x, max_pix_y])


def PlottingFunc(ax,data,label,fs):
    ax.imshow( np.log10( data ), cmap=plt.cm.plasma,vmin=20.5,vmax=21.5)
    ax.annotate(r'{}'.format(label),xy=(50-3, 150-3),fontsize=fs,color='black',xycoords='data')
    ax.annotate(r'{}'.format(label),xy=(50, 150),fontsize=fs,color='white',xycoords='data')
    #ax[0].plot(scale_bar[0],scale_bar[1],color='black',linewidth=2)
    ax.set_xticks([])
    ax.set_yticks([])

    Ymax, Xmax = np.where(data == data.max())

    # Black shadow
    ax.scatter(Xmax[0]-3, Ymax[0]-3, marker ='o', color='black',s=2)
    ax.annotate(r'$\Sigma_{max}$',xy=(Xmax[0] - 23, Ymax[0] - 23),fontsize=fs,
                color='black',xycoords='data')
    # White
    ax.scatter(Xmax[0], Ymax[0], marker='o', color='white',s=2)
    ax.annotate(r'$\Sigma_{max}$',xy=(Xmax[0] - 20, Ymax[0] - 20),fontsize=fs,
                color='white',xycoords='data')

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

Pol_qui, mean_qui, var_qui     = ColDensityFITSextract("herschel_polaris_coldens_18_quiet.fits",dataAdd)
Pol_sax, mean_sax, var_sax     = ColDensityFITSextract("herschel_polaris_coldens_18_saxophone.fits",dataAdd)
Pol_ful, mean_ful, var_ful     = ColDensityFITSextract("herschel_polaris_coldens_36.fits",dataAdd)


# Plots
###########################################################################

# pix_1pc     = 1/pix_pc[0]
# dx          = 10
# dy          = 600
# scale_bar   = np.array([[dx,dx+pix_1pc],[dy,dy]])

f, ax = plt.subplots(1,2, figsize=(9, 3), dpi=250, facecolor='w')
f.subplots_adjust(hspace=0.01)
fs = 12;

PlottingFunc(ax[0],Pol_sax,'Saxophone Region',fs)
PlottingFunc(ax[1],Pol_qui,'Quiet Region',fs)





#cbar1 = f.colorbar(Pol_sax,ax = ax[0],pad=0.01)
#cbar1.set_label(r"$\log_{10} \Sigma$ [cm$^{-2}$]",fontsize=fs)

#cbar2 = f.colorbar(Pol_qui,ax = ax[1],pad=0.01)
#cbar2.set_label(r"$\log_{10} \Sigma$ [cm$^{-2}$]",fontsize=fs)
plt.tight_layout()


plt.show()
