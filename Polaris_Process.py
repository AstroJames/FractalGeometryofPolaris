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
    This function creates a sliding window along the vertical axis of the image,
    finding the coordinates for the maximum value for each windowself.

    INPUTS:
    image           - the 2D image data, should be a np.array.
    windowSize      - the size of the window in pixels (along the vertical axis).
    windowTop       - where the window should start in pixels coordinates (for clouds that don't begin )
                        at the top boundary.
    windowBottom    - where the window should terminate in pixel cooridnates.

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


def PlottingFunc(ax,data,label,fs,axisExtent,colorbar=None,refPixel=None):
    """
    This function plots the saxophone and quite subregions of polaris.

    INPUTS:
    ax          - the axis from matplotlib.
    data        - the 2D column density data.
    label       - the label for the subregiong.
    fs          - the fontsize for the labels.
    colobar     - the colorbar (set to true if you want one).
    refPixel    - the reference pixel (set to true if one exists).

    """

    # Axis labels
    numOfLabels = 10
    x           = axisExtent[0]
    y           = axisExtent[1]
    xStep       = int(data.shape[0]/ (numOfLabels-1))
    yStep       = int(data.shape[1]/ (numOfLabels-1))
    xPositions  = np.arange(0,data.shape[0],xStep)
    yPositions  = np.arange(0,data.shape[1],yStep)
    xLabels     = x[::xStep]
    yLabels     = y[::yStep]


    # The 2D column density and the labels for the regions.
    img = ax.imshow( np.log10( data ), cmap=plt.cm.plasma,vmin=20.5,vmax=21.5,interpolation='none')
    ax.annotate(r'{}'.format(label),xy=(50-3, 150-3),fontsize=fs,color='black',xycoords='data')
    ax.annotate(r'{}'.format(label),xy=(50, 150),fontsize=fs,color='white',xycoords='data')

    # Turn of axis ticks for now.
    ax.set_xticks([])
    ax.set_yticks([])

    # The global maximum value of the column density.
    Ymax, Xmax = np.where(data == data.max())

    # The label for the maximum column density and the pixel coordinate for the density.
    ax.scatter(Xmax[0]-3, Ymax[0]-3, marker ='o', color='black',s=3)
    ax.annotate(r'$\Sigma_{max}$',xy=(Xmax[0] - 23, Ymax[0] - 23),fontsize=fs+2,
                color='black',xycoords='data')
    ax.scatter(Xmax[0], Ymax[0], marker='o', color='white',s=3)
    ax.annotate(r'$\Sigma_{max}$',xy=(Xmax[0] - 20, Ymax[0] - 20),fontsize=fs+2,
                color='white',xycoords='data')

    # Add a colourbar if the user wants one.
    if colorbar is not None:
        cbar = plt.colorbar(img, ax=ax,pad=0.01)
        cbar.set_label(r"$\log_{10} \Sigma$ [cm$^{-2}$]",fontsize=fs,labelpad=20,rotation=270)

    # Add a reference pixel if the user wants one.
    if refPixel is not None:
        ax.scatter(refPixel[0],refPixel[1],color='white',marker='*',s=3)


def CoordFITSExtract(file,dataAddress,density):
    """
    This function extracts the coordiante information from the FITS file.

    file        - the file name.
    dataAddress - the directory of the FITS files.
    density     - the 2D column denisty.

    """
    data    = fits.open(dataAddress+file)

    # some info on coordinates:
    x_axis          = data[0].header["NAXIS1"]
    y_axis          = data[0].header["NAXIS2"]
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

    # Create an array for the y axis in offset coordinates
    if (density.shape[0] == y_axis):
        print("The y-axis is correct.")

        y_axis = ( coord_y_pix  - np.array(range(0,density.shape[0]-1)) ) * data[0].header["CDELT2"]

    # Create an array for the x axis in the offset coordinates
    if (density.shape[1] == x_axis):
        print("The x-axis is correct.")

        x_axis = ( coord_x_pix  - np.array(range(0,density.shape[1]-1)) ) * data[0].header["CDELT1"]


    return (coord_x_offset, coord_y_offset), (x_axis, y_axis), pix_pc


# Coordinates
###########################################################################



# Column Density
###########################################################################
# Saxaphone
# NAXIS1  = 1567
# NAXIS2  = 1948

# Quiet
# NAXIS1  = 2000
# NAXIS2  = 1366


# Working Script
###########################################################################
quietFile   = "herschel_polaris_coldens_18_quiet.fits"
saxFile     = "herschel_polaris_coldens_18_saxophone.fits"
fullFile    = "herschel_polaris_coldens_36.fits"

Pol_qui, mean_qui, var_qui              = ColDensityFITSextract(quietFile,dataAdd)
coord_offset_qui, axis_qui, pix_pc_qui  = CoordFITSExtract(quietFile,dataAdd,Pol_qui)

Pol_sax, mean_sax, var_sax              = ColDensityFITSextract(saxFile,dataAdd)
coord_offset_sax, axis_sax, pix_pc_sax  = CoordFITSExtract(saxFile,dataAdd,Pol_sax)

# Plots
###########################################################################

pix_1pc     = 1/pix_pc_sax  # the number of pixels required for a parsec
dx          = 150           # coordinates for the scale bar
dy          = 1800          # coordinates for the scale bar
scale_bar   = np.array([[dx,dx+pix_1pc],[dy,dy]])

f, ax = plt.subplots(1,2, figsize=(4, 4), dpi=250, facecolor='w')
f.subplots_adjust(hspace=0.01,wspace=0)
fs = 12;

PlottingFunc(ax[0],Pol_sax,'Saxophone Subregion',fs,axis_sax)
ax[0].annotate('1pc',xy=(315,1750),color='white',fontsize=fs-2)
ax[0].plot(scale_bar[0],scale_bar[1],color='white',linewidth=2)

PlottingFunc(ax[1],np.flipud(np.transpose(np.flipud(Pol_qui))),'Quiet Subregion',fs,axis_qui,colorbar=True)
#cbar1.set_label(r"$\log_{10} \Sigma$ [cm$^{-2}$]",fontsize=fs)

#cbar2 = f.colorbar(Pol_qui,ax = ax[1],pad=0.01)
#cbar2.set_label(r"$\log_{10} \Sigma$ [cm$^{-2}$]",fontsize=fs)
#plt.tight_layout()


plt.show()
