"""
Title:          Fractal Dimension Calculation on Polaris FITS data
Author:         James Beattie
Creation Date:  12th Feb, 2019
"""

# Imports
###########################################################################

import numpy as np
from astropy.coordinates import SkyCoord
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
    """


    """

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


def PlottingFunc(ax,data,label,fs,colorbar=None):
    """

    """

    img = ax.imshow( np.log10( data ), cmap=plt.cm.plasma,vmin=20.5,vmax=21.5,interpolation='none')
    ax.annotate(r'{}'.format(label),xy=(50-3, 150-3),fontsize=fs,color='black',xycoords='data')
    ax.annotate(r'{}'.format(label),xy=(50, 150),fontsize=fs,color='white',xycoords='data')
    #ax[0].plot(scale_bar[0],scale_bar[1],color='black',linewidth=2)
    ax.set_xticks([])
    ax.set_yticks([])

    Ymax, Xmax = np.where(data == data.max())

    # Black shadow
    ax.scatter(Xmax[0]-3, Ymax[0]-3, marker ='o', color='black',s=3)
    ax.annotate(r'$\Sigma_{max}$',xy=(Xmax[0] - 23, Ymax[0] - 23),fontsize=fs+2,
                color='black',xycoords='data')
    # White
    ax.scatter(Xmax[0], Ymax[0], marker='o', color='white',s=3)
    ax.annotate(r'$\Sigma_{max}$',xy=(Xmax[0] - 20, Ymax[0] - 20),fontsize=fs+2,
                color='white',xycoords='data')


    if colorbar is not None:
        cbar = plt.colorbar(img, ax=ax,pad=0.01)
        cbar.set_label(r"$\log_{10} \Sigma$ [cm$^{-2}$]",fontsize=fs,labelpad=20,rotation=270)


def CoordFITSExtract(file,dataAddress):
    """

    """
    data    = fits.open(dataAddress+file)

    # some info on coordinates:
    coord_x_pix     = data[0].header["CRPIX1"]
    coord_y_pix     = data[0].header["CRPIX2"]
    coord_x_offset  = data[0].header["CRVAL1"]
    coord_y_offset  = data[0].header["CRVAL2"]

    # compute the size of a pixel
    pix             = data[0].header["CDELT2"]*3600.0       # in arcsec
    distance        = 0.15                                  # kpc
    arcsec2rad      = 206265.0                              # arcseconds to radians
    pix_pc          = distance * 1000.0 * pix/arcsec2rad    # pixels to parsecs
    print 'pixel size in parsecs = ', pix_pc

    return (coord_x_pix, coord_y_pix), (coord_x_offset, coord_y_offset), pix_pc


# Coordinates
###########################################################################



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

Pol_qui, mean_qui, var_qui                  = ColDensityFITSextract("herschel_polaris_coldens_18_quiet.fits",dataAdd)
coord_pix_qui, coord_offset_qui, pix_pc_qui = CoordFITSExtract("herschel_polaris_coldens_18_quiet.fits",dataAdd)

Pol_sax, mean_sax, var_sax                  = ColDensityFITSextract("herschel_polaris_coldens_18_saxophone.fits",dataAdd)
coord_pix_sax, coord_offset_sax, pix_pc_sax = CoordFITSExtract("herschel_polaris_coldens_18_saxophone.fits",dataAdd)

# Plots
###########################################################################

pix_1pc     = 1/pix_pc_sax  # the number of pixels required for a parsec
dx          = 150           # coordinates for the scale bar
dy          = 1800          # coordinates for the scale bar
scale_bar   = np.array([[dx,dx+pix_1pc],[dy,dy]])

f, ax = plt.subplots(1,2, figsize=(4, 4), dpi=250, facecolor='w')
f.subplots_adjust(hspace=0.01,wspace=0)
fs = 12;

PlottingFunc(ax[0],Pol_sax,'Saxophone Subregion',fs,None)
ax[0].annotate('1pc',xy=(315,1750),color='black',fontsize=fs-2)
ax[0].plot(scale_bar[0],scale_bar[1],color='black',linewidth=2)

PlottingFunc(ax[1],np.flipud(np.transpose(np.flipud(Pol_qui))),'Quiet Subregion',fs,True)
#cbar1.set_label(r"$\log_{10} \Sigma$ [cm$^{-2}$]",fontsize=fs)

#cbar2 = f.colorbar(Pol_qui,ax = ax[1],pad=0.01)
#cbar2.set_label(r"$\log_{10} \Sigma$ [cm$^{-2}$]",fontsize=fs)
#plt.tight_layout()


plt.show()
