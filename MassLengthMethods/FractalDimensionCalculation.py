"""
Title:      Mass-Length Fractal Dimension Curve Constructor
Author:     James Beattie
Comments:   Implementation of Beattie et al. 2019 on Polaris FITS data.

"""

##############################################################################################################################
# Dependencies
##############################################################################################################################

import numpy as np                             # mathematics
import argparse                                # command line arguments
from skimage import measure
import imp
import sys
from astropy.io import fits
import matplotlib.patches as patches
from matplotlib import pyplot as plt
from matplotlib import rc
import FractalDimensionUtilities
from FractalDimensionUtilities import *         # FractalDimensionUtilities module written by me, for me.
imp.reload(FractalDimensionUtilities)

##############################################################################################################################
# Comments
##############################################################################################################################

# last edited:
# 1 May, 2018: changed the sampling to be linear
# 5 June, 2018: May have found a bug in the boundary conditions.
# 6 January, 2019: editting for application on the Brick

# Example Command line:
""" run FractalDimensionCalculation -lmin 3 -boundaries 'Terminating' -FITSfile 'FITS.fits' """

##############################################################################################################################
# Any functions for processing: brick.dust_column_density_cf.fits
##############################################################################################################################

def MaxPixelWindow(image,windowsize,method,windowBottom,windowTop):
    """
    This function is used to get all of the maximum values of the column density,
    either in a vertically sliding window, or a single, global maximum


    INPUTS:
    image           - the 2D column density data.
    windowsize      - the sie of the window in pixels.
    method          - the method, either: "SlidingWindow" or "MaxPixel".
    windowBottom    - the bottom boundary of the sliding window.
    windowTop       - the top boundary of the sliding window.
    """

    # Error Handling
    method_options = ["SlidingWindow", "MaxPixel"]
    if method not in method_options:
        raise ValueError("You have not input a correct method option.")
        sys.exit(1)


    max_pix_x   = []                            # max pixel coordinates for x
    max_pix_y   = []                            # max pixel coordinates for y
    start       = 0                             # bottom y coordinate for the sliding window
    top         = image.shape[0]/windowsize     # top y coordinate for the sliding window

    # These are the coordinates for the brick
    #windowTop    = 70
    #windowBottom = 620

    # Sliding window method
    if method == "SlidingWindow":
        # Slide the window down the brick image, searching for the maximum pixel in each window
        for end in xrange(1,top):
            if windowsize*end > windowTop and windowsize*end < windowBottom:
                y,x = np.where(image[start:windowsize*end,:] == image[start:windowsize*end,:].max())
                max_pix_y.append(y[0]+start)
                max_pix_x.append(x[0])
                start += windowsize
            else:
                start += windowsize
                continue

    # Single Max pixel method
    elif method == "MaxPixel":
        y, x = np.where(image[windowTop:windowBottom,:] == image[windowTop:windowBottom,:].max())
        max_pix_y.append(y[0]+start)
        max_pix_x.append(x[0])

    return max_pix_x, max_pix_y

def FITSRead(filename,RBthreshold,CFthreshold):
    image           = fits.open(filename)             # open fits file
    image           = np.flipud(image[0].data)        # Federrath 2016 preprocessed
    RB_mask         = image >= RBthreshold            # sensitivity treshold in JR2014
    CF_sigma        = image*(image >= CFthreshold)    # indexes of the contour in CF2016

    return CF_sigma, RB_mask

def ColDensityFITSextract(file):
    """
    This reads in FITS files and immediatly extracts the column density data from themself.

    INPUTS:
    file            - the file name
    """

    dataAddress = "../FITSfiles/"

    # Read in the FITS files and extract column density.
    data    = fits.open(dataAddress+file)
    data    = data[0].data[0,0,:,:]

    # Check if the data is actually just the column denisty
    if data.ndim != 2:
        print("WARNING: The dimensions of the FITS read is not 2.")

    return data

##############################################################################################################################
# Command Line Arguments
##############################################################################################################################

ap 			= argparse.ArgumentParser(description = 'Just a bunch of input arguments');
ap.add_argument('-FITSfile', '--FITSfile',required=True, help = 'the FITS file', type=str);
ap.add_argument('-boundaries', '--boundaries',required=True, help = 'the boundary type: "Terminating" or "Periodic" ', type=str);
ap.add_argument('-lmin', '--lmin',required=True, help = 'the intial length of the box', type=int);
args 		= vars(ap.parse_args());

##############################################################################################################################
# Error handling
##############################################################################################################################

readtype_options = ["RegionCentroid", "RegionMaxpixel", "Maxpixel"]
boundary_options = ["Terminating", "Periodic"]

if args['boundaries'] not in boundary_options:
    raise ValueError("You have not input a correct boundary option.")
    sys.exit(1)

# Make sure that dxdy is odd.
if np.mod(args['lmin'],2) == 0:
     args['lmin'] += 1;
     print "I have increased lmin by 1, because it needs to be odd."

##############################################################################################################################
# Mass-Length Algorithm
##############################################################################################################################

#############################
# Initialisations
#############################

dxdy                = args['lmin'];         # the length of the square
lmin                = args['lmin']          # store the minimum l for sampling
dxdy_arr            = [];                   # store the dxdy values
log_dxdy_arr        = [];                   # store log values for plotting
log_mass_arr        = [];                   # store the mass.
mass_arr            = [];                   # store the masses
std_arr             = [];                   # standard deviation array.
box_count           = [];                   # store the amount of lxl regions
fractal_dim         = [];                   # store the fractal dimensions
iter_count          = 0;                    # just a count of the iterations in the while loop
Box_state           = True;                 # this will be used to break the while loop if there are no longer any boxes left
Box_state_count     = 0;                    # box state counter required for the iteration
terminator_count    = 0;                    # count how many boxes get terminated, using terminating BCs
Region_label        = []                    # label the region that terminates
JRthreshold         = 2.0*(25e-3*1.9e23)    # see Ratheborne et al. (2014, page 2, 3 left column)
CFthreshold         = 5e22                  # cm^{-2} # CF2016, contour

image       = ColDensityFITSextract(args['FITSfile'])
lmax        = image.shape[1];   # the maximum size of the squares for the counting method

#############################
# Calculate Centroids
#############################

centroid_x, centroid_y = MaxPixelWindow(image,50,"MaxPixel",image.shape[1],0)

#############################
# Initialise Figure
#############################

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

#pix_pc      = 0.01408382    # pixels per parsec
#pix_1pc     = 1/pix_pc   # 1 parsec in pixels
#dx          = 10
#dy          = 600
#scale_bar   = np.array([[dx,dx+pix_1pc],[dy,dy]])
fs          = 20

fig, ax = plt.subplots()
plt.scatter(centroid_x,centroid_y,marker='o',color='white')
plt.imshow( np.log10( image ),cmap=plt.cm.plasma,vmax=22,vmin=20.5)
#plt.annotate(r'ALMA + $Herschel$',xy=(10, 50),fontsize=fs,color='black')
#plt.annotate(r'G0.253+0.016',xy=(10, 65),fontsize=fs,color='black')
#plt.annotate(r'1pc',xy=(45, 595),fontsize=fs-2,color='black')
#plt.plot(scale_bar[0],scale_bar[1],color='black',linewidth=2)
ax.set_xticks([])
ax.set_yticks([])
cbar = plt.colorbar()
cbar.set_label(r"$\log_{10} \Sigma$ [cm$^{-2}$]",fontsize=fs,labelpad=20)

#############################
# Region Expansion
#############################

while dxdy < lmax and Box_state is True:
    x                   = 0                   # initialise the x coordinate
    y                   = 0                   # initialise the y coordinate
    num_boxes           = len(centroid_x)     # total number of boxes equals total amount of centroids
    box_counter         = 0                   # initialise the box counter
    Mass_region         = [];


    # for each set of centroid coordinates
    for i in xrange(0,num_boxes):

        # calcuate where the median of the odd dxdy is to center the box on the centroid.
        x = int(centroid_x[i] - ((dxdy-1)/2 + 1))
        y = int(centroid_y[i] - ((dxdy-1)/2 + 1))

        ####################################################
        # Conditions for terminating boundary conditions
        ####################################################

        if (args['boundaries'] == "Terminating"):

            # Skip on regions that have already terminated
            if any(np.array(Region_label) == i):
                continue

            #if the entire region is outside of the mask or if the expansion goes past the dimensions

            if y < 0 or x < 0 or y > image.shape[0] or x > image.shape[1] or x + dxdy > image.shape[1]  or y + dxdy > image.shape[0]:
                # move onto the next box if anything fails
                terminator_count += 1   # store the termination
                Region_label.append(i)  # store the box label

                # Calculate old x and y
                dxdy_old    = dxdy - 2
                x_old       = int(centroid_x[i] - ((dxdy_old-1)/2 + 1))
                y_old       = int(centroid_y[i] - ((dxdy_old-1)/2 + 1))

                rect        = patches.Rectangle((x_old,y_old),dxdy_old,dxdy_old,linewidth=1,edgecolor='white',fill=False);
                ax.add_patch(rect);
                continue
            else:
                # else just sum the column densities in the region
                Mass_sum    = sum(sum(image[y:y+dxdy,x:x+dxdy]));
                Box_state_count += 1
                box_counter += 1;
                Mass_region.append(Mass_sum);


                if np.mod(dxdy,100) == 1:
                    rect        = patches.Rectangle((x,y),dxdy,dxdy,linewidth=1,edgecolor='white',fill=False);
                    ax.add_patch(rect);

                #if dxdy == 225:
                #    rect        = patches.Rectangle((x,y),dxdy,dxdy,linewidth=1,edgecolor='r',facecolor='r',alpha=0.1);
                #    ax.add_patch(rect);

        ####################################################
        # Conditions on x for periodic boundary conditions
        ####################################################
        # if the new length from the coordinate is greater than the image width, i.e. the RHS
        # We have to split the sum over the pixel values up across the two rectangular regions.

        elif (args['boundaries'] == "Periodic"):

            if x + dxdy > image.shape[0] and y >= 0 and x >= 0 and y + dxdy <= image.shape[1]:
                #print('x + dxdy > img.shape[0]')
                Mass_Psum_1    = sum(sum(image[y:y+dxdy,x:image.shape[0]]));

                # the partial sum over the second rectangular region
                Mass_Psum_2    = sum(sum(image[y:y+dxdy,0:(dxdy-(image.shape[0]-x))]));

                Mass_sum       = Mass_Psum_1 + Mass_Psum_2;


            # if the coordinate passes through the LHS
            elif x < 0 and y >= 0 and y + dxdy <= image.shape[1] and x + dxdy <= image.shape[0]:
                #print('x < 0')
                Mass_Psum_1    = sum(sum(image[y:y+dxdy,np.mod(x,image.shape[0]):image.shape[0]]));

                delx           = (dxdy-(image.shape[0]-np.mod(x,image.shape[0])))
                Mass_Psum_2    = sum(sum(image[y:y+dxdy,0:delx]));

                Mass_sum       = Mass_Psum_1 + Mass_Psum_2;


            # if the y or the y + dxdy is outside of the image, ignore it.
            elif y + dxdy > image.shape[1] and x >= 0 and y >= 0 and x + dxdy <= image.shape[0]:
                #print('y + dxdy > img.shape[1]')
                # the partial sum over the first rectangular region
                Mass_Psum_1    = sum(sum(image[y:image.shape[1],x:x+dxdy]));

                # the partial sum over the second rectangular region
                Mass_Psum_2    = sum(sum(image[0:(dxdy-(image.shape[1]-y)),x:x+dxdy]));

                Mass_sum       = Mass_Psum_1 + Mass_Psum_2;


            elif y < 0 and x >= 0 and y + dxdy <= image.shape[1] and x + dxdy <= image.shape[0]:
                #print('y < 0')
                # the partial sum over the first rectangular region
                Mass_Psum_1    = sum(sum(image[0:(dxdy-(image.shape[1]-np.mod(y,image.shape[1]))),x:x+dxdy]));

                # the partial sum over the second rectangular region
                Mass_Psum_2    = sum(sum(image[(image.shape[1]+y):image.shape[1],x:x+dxdy]));

                Mass_sum       = Mass_Psum_1 + Mass_Psum_2;

            elif y < 0 and x < 0:
                #print('y < 0 and x < 0')

                yP1   = 0;
                dyP1  = dxdy - (image.shape[1] - np.mod(y,image.shape[1]));

                xP1   = 0;
                dxP1  = dxdy - (image.shape[0] - np.mod(x,image.shape[0]));

                # the partial sum over the first rectangular region
                Mass_Psum_1    = sum(sum(image[yP1:dyP1,xP1:dxP1]));


                yP2   = 0;
                dyP2  = dxdy - (image.shape[1] - np.mod(y,image.shape[1]));

                xP2   = np.mod(x,image.shape[0]);
                dxP2  = image.shape[0];
                # the partial sum over the second rectangular region
                Mass_Psum_2    = sum(sum(image[yP2:dyP2,xP2:dxP2]));

                # the partial sum over the second rectangular region
                Mass_Psum_3    = sum(sum(image[np.mod(y,image.shape[1]):image.shape[1],0:dxP1]));

                # the partial sum over the second rectangular region
                Mass_Psum_4    = sum(sum(image[np.mod(y,image.shape[1]):image.shape[1],np.mod(x,image.shape[0]):image.shape[0]]));

                Mass_sum       = Mass_Psum_1 + Mass_Psum_2 + Mass_Psum_3 + Mass_Psum_4;


            elif y + dxdy > image.shape[1] and x + dxdy > image.shape[0]:
                #print('y + dxdy > img.shape[1] and x + dxdy > img.shape[0]')
                # the partial sum over the first rectangular region
                Mass_Psum_1    = sum(sum(image[y:image.shape[1],x:image.shape[0]]));

                # the partial sum over the second rectangular region
                Mass_Psum_2    = sum(sum(image[y:image.shape[1],0:(dxdy-(image.shape[0]-x))]));

                # the partial sum over the second rectangular region
                Mass_Psum_3    = sum(sum(image[0:(dxdy-(image.shape[1]-y)),0:(dxdy-(image.shape[0]-x))]));

                # the partial sum over the second rectangular region
                Mass_Psum_4    = sum(sum(image[0:(dxdy-(image.shape[1]-y)),x:image.shape[0]]));

                Mass_sum       = Mass_Psum_1 + Mass_Psum_2 + Mass_Psum_3 + Mass_Psum_4;


            elif y+dxdy > image.shape[1] and x < 0:
                #print('y + dxdy > img.shape[1] and x < 0')
                Mass_Psum_1    = sum(sum(image[y:image.shape[1],np.mod(x,image.shape[0]):image.shape[0]]));

                delx           = (dxdy-(image.shape[0]-np.mod(x,image.shape[0])))
                Mass_Psum_2    = sum(sum(image[y:y+image.shape[1],0:delx]));
                # the partial sum over the second rectangular region
                Mass_Psum_3    = sum(sum(image[0:dxdy-(image.shape[1]-y),0:delx]));
                # the partial sum over the second rectangular region
                Mass_Psum_4    = sum(sum(image[0:dxdy-(image.shape[1]-y),np.mod(x,image.shape[0]):image.shape[0]]));

                Mass_sum       = Mass_Psum_1 + Mass_Psum_2 + Mass_Psum_3 + Mass_Psum_4;

            elif x+dxdy > image.shape[0] and y < 0:
                #print('x + dxdy > img.shape[0] and y < 0')

                Mass_Psum_1    = sum(sum(image[0:(dxdy - (image.shape[1] - np.mod(y,image.shape[1]))),x:image.shape[0]]));

                Mass_Psum_2    = sum(sum(image[0:(dxdy - (image.shape[1] - np.mod(y,image.shape[1]))),0:(dxdy-(image.shape[0]-x))]));

                Mass_Psum_3    = sum(sum(image[np.mod(y,image.shape[1]):image.shape[1],x:image.shape[0]]));

                Mass_Psum_4    = sum(sum(image[np.mod(y,image.shape[1]):image.shape[1],0:(dxdy-(image.shape[0]-x))]));

                Mass_sum       = Mass_Psum_1 + Mass_Psum_2 + Mass_Psum_3 + Mass_Psum_4;

            else:
                Mass_sum    = sum(sum(image[y:y+dxdy,x:x+dxdy]));

            box_counter += 1;
            Mass_region.append(Mass_sum);

    #############################
    # Exit Statement
    #############################

    # Control the boxcounting state. Terminate the while loop is there were no
    # boxes counted in the for loop iteration
    if Box_state_count == 0:
        Box_state = False      # set the box state to false
        print('length:')
        print(dxdy_arr)
        print('\n')
        print('box_count:')
        print(box_count)
        print('\n')
        print('mass:')
        print(mass_arr)
        print('\n')
        print('std:')
        print(std_arr)
        print('\n')
        print('fractal_dim:')
        print(fractal_dim)

        plt.show()

        continue
    else:
        Box_state_count = 0

    print("{} boxes have been terminated.".format(terminator_count))
    print("There were a total of {} boxes in the l = {} iteration ({}) \n".format(box_counter,dxdy,iter_count))

    #############################
    # Construct Fractal Dimension
    #############################

    # Calculate the mean over each of the squares
    Mass_mean  = np.mean(Mass_region);
    Mass_std   = np.std(Mass_region);

    # Store the amount of boxes and the size of the boxes

    log_mass_arr.append(np.log10(Mass_mean));
    log_dxdy_arr.append(np.log10(dxdy)); # in units parsecs


    # After the first iteration calculate the linear regression and fractal dimension
    if iter_count >= 1:
        #print('Starting regression on iteration {}'.format(iter_count))

        dxdy_arr.append(dxdy); # in units parsecs
        mass_arr.append(Mass_mean);
        std_arr.append(Mass_std);
        box_count.append(len(Mass_region));

        slope, intercept, domain, abline_values, log_log_predict, X, Y = linear_regression(log_mass_arr,log_dxdy_arr,iter_count);

        # Storing all of the predictions
        if iter_count == 1:
            predict_storage = np.array(abline_values);
        else:
            predict_storage = np.vstack([predict_storage,abline_values]);

        # Store the fractal dimension
        fractal_dim.append(slope);

    iter_count          += 1    # add one to the iteration counter
    dxdy                += 2    # sample uniformly.

    # Make sure that dxdy is always odd for symmetry about the centroid
    if np.mod(dxdy,2) == 0:
        dxdy += 1

maxdxdy = np.double(max( np.array(dxdy_arr) ))

plt.show()
plt.plot(np.log10(dxdy_arr/maxdxdy ),fractal_dim)
#plt.axhline(min(fractal_dim),color='red')
plt.xlabel(r"$\log_{10}\left(\ell/\ell_{max}\right)$",fontsize=16)
plt.ylabel(r"$\mathcal{D}(\ell/\ell_{max})$",fontsize=16)
plt.show()
