import numpy as np;
import heapq;

import imageio;                                 # reading in mp4 data
import h5py;                                    # importing in hdf5 files

# image processing modules
import skimage;                                 # import image data
from skimage import filters, measure;           # some more image processing tools for getting area and perimeter data

# Statistics modules
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score


# Suppress warnings
import warnings
warnings.filterwarnings(action="ignore", module="scipy",message="^internal gelsd")
warnings.filterwarnings(action="ignore", module="sklearn", message="^internal gelsd")

# Classes
##############################################################################################################################
class h5py_data:
    density    = [];
    time       = [];
    direction  = [];
    minmax_xyz = [];
    velocity   = [];

# Reading functions
##############################################################################################################################
def h5py_reader(data,readtype):

    if readtype == 'all':
        h5_output   = os.system('h5ls ' + data);
        print '\n'
        f           = h5py.File(data, 'r+');
        h5py_data.density     = f['dens_proj_xy'];
        h5py_data.time        = f['time'];
        h5py_data.direction   = f['direction'];
        h5py_data.minmax_xyz  = f['minmax_xyz'];
        h5py_data.velocity    = (f['velx_proj_xy'],f['vely_proj_xy'],f['velz_proj_xy']);
        f.close;
    elif readtype == 'projection':
        f           = h5py.File(data, 'r+');
        h5py_data.density     = f['dens_proj_xy'];
        h5py_data.time        = f['time'];
        h5py_data.direction   = f['direction'];
        h5py_data.minmax_xyz  = f['minmax_xyz'];
        f.close;
    elif readtype == 'slice':
        f           = h5py.File(data, 'r+');
        h5py_data.density     = f['dens_slice_xy'];
        h5py_data.time        = f['time'];
        h5py_data.direction   = f['direction'];
        h5py_data.minmax_xyz  = f['minmax_xyz'];
    elif readtype == 'slice_highres':
        f           = h5py.File(data, 'r+');
        h5py_data.density     = f['slice_dens'];
        h5py_data.time        = f['time'];
        h5py_data.direction   = f['direction'];
        h5py_data.minmax_xyz  = f['minmax_xyz'];
    else:
        print('Invalid readtype argument.')
    return h5py_data


def reader_functions(video,hdf5_file):
    vid         = [];
    hdf5_list   = [];
    # Video reader
    if video is not None:
        vid = read_video()
    # HDF5 reader
    elif hdf5_file is not None:

        if hdf5_file == 'Data.txt':

            f = open(hdf5_file,'r');

            for line in f:
                hdf5_list.append(line)
    # No reader
    else:
        print 'No reading present.'

    return vid, hdf5_list


def read_and_threshold(video_args,hdf5_args,readtype_args,threshold_args):

    vid, hdf5_list  = reader_functions(video_args,hdf5_args);

    if video_args is not None:# Process image, scale image, mean, sd of image
        image           = vid.get_data(num);
        image_thres     = image[:,:,2];
    else:
        # If I pass it a text file it will iterate through the text file.
        if hdf5_args == 'Data.txt':
            h5py_data       = h5py_reader(hdf5_list[num][0:22],readtype_args);    # read the first line of Data.txt
            # create the density PDF
            #image           = np.log(h5py_data.density[:,:,0]/np.mean(h5py_data.density[:,:,0].ravel()));  # read in the density image.
            image           = h5py_data.density[:,:,0];  # read in the density image.
            image_thres     = image;
        # If I pass it a single file just read the single file
        else:
            h5py_data       = h5py_reader(hdf5_args,readtype_args);    # read the first line of Data.txt
            # create the density PDF
            #image           = np.log10(h5py_data.density[:,:]/h5py_data.density[:,:].mean());  # read in the density image.
            image           = h5py_data.density[:,:];  # read in the density image.
            image_thres     = image;
            time            = h5py_data.time[:];

    mask_3d, mask   = thresholding_preprocess(image_thres,threshold_args);

    return mask_3d, mask, image, image_thres, time

# Preprocessing functions
##############################################################################################################################

def thresholding_preprocess(image,threshold):
    #print('Thresholding the 2D data')

    if threshold == 'triangle':
        val         = filters.threshold_triangle(image);
    elif threshold == 'mean':
        val         = filters.threshold_mean(image);
    elif threshold == 'MassLength': # mass perimenter method
        val         = 3*image.max()/4
    else:
        print('You need to select either triangle, mean or MassPerimeter.')

    mask        = image > val;
    mask_ones   = np.ones([mask.shape[0],mask.shape[1]]);
    mask_inv    = (-1)*(mask*mask_ones-1);
    mask_inv3d  = np.repeat(mask_inv[:, :, np.newaxis], 3, axis=2)
    return mask_inv3d, mask


def contour_finder(image):
    # Create the image threshold
    image_float                     = image.astype(np.float)
    image_guass                     = filters.gaussian(image_float,sigma=10);
    contours                        = measure.find_contours(image_guass, 0.2);
    return contours


def label_creator(image):
    #print('Using segmentation on the image')
    #image_guass                     = filters.gaussian(image,sigma=10);
    all_labels                      = measure.label(image)
    return all_labels


def read_video(vid_file):
    vid     = imageio.get_reader(vid_file,  'ffmpeg');
    return vid


def linear_regression(Y,X,iteration):
    #print('Least squares fitting on iteration {}'.format(iteration+1))

    Y      = np.transpose([Y]);
    X      = np.transpose([X]);


    # Perform the linear regression
    regr            = linear_model.LinearRegression();
    regr.fit(X,Y);
    log_log_predict = regr.predict(X);

    # Just make a bunch of values over the entire x domain for the line
    slope           = regr.coef_[0][0];
    intercept       = regr.intercept_[0];
    domain          = np.arange(0,20);
    abline_values   = [slope * i + intercept for i in domain];

    return slope, intercept, domain, abline_values, log_log_predict, X, Y


def centroid_coordinates(mask,image,numofpeaks,type):
    #print('Calculating centroids')

    centroid_x              = [];       # the x array of centroid coordinates
    centroid_y              = [];       # the y array of centroid coordinates
    max_pixel_x             = [];       # max pixel x coordinates
    max_pixel_y             = [];       # max pixel y coordinates
    image_pixels            = [];       # storing image pixels from region coordinates.
    centroid_coord          = [];       # initialising the centroid coordinate
    area                    = [];

    if type == 'Maxpixel':
        copy = image.copy()

        for i in xrange(numofpeaks):


            # find the maximum pixel
            y,x = np.where(copy==copy.max());
            centroid_x.append(x[0]);   # add it to centroid x
            centroid_y.append(y[0]);   # add it to centroid y
            copy[y[0],x[0]] = 0 # set it to 0 to find the next global maximum

        centroid_coord = centroid_x # just create a dummy max peak vector


    else:

        # Generate the centroids from segmented areas in the mask
        centroid_coord          = measure.regionprops(measure.label(mask));

        # Down-sample to a uniform sample from the array
        if len(centroid_coord) > numofpeaks: # make sure I only use the number of peaks specified
            #centroid_downsample.append(np.random.choice(centroid_coord,numofpeaks,replace=False));

            # Take the largest area regions that are all above the threshold
            for i in xrange(len(centroid_coord)):
                area.append(centroid_coord[i].area)

            # Take the indexes on those regions and fill the downsampled centroid vector
            indexes                 = map(area.index,heapq.nlargest(numofpeaks, area))
            centroid_coord          = [centroid_coord[x] for x in indexes]


        if type == 'RegionCentroid':

            for i in xrange(0,len(centroid_coord)):
                y0,x0  = centroid_coord[i].centroid
                centroid_y.append(y0)
                centroid_x.append(x0)

        elif type == 'RegionMaxpixel':

            for i in xrange(0,len(centroid_coord)):

                y0 = centroid_coord[i].coords[xrange(0,len(centroid_coord[i].coords))][:,1]
                x0 = centroid_coord[i].coords[xrange(0,len(centroid_coord[i].coords))][:,0]

                # store all image pixels for the coordinates of that region
                image_pixels = list(image[y0,x0]);

                # calculate the maximum pixel
                maxpixel   = np.max(image_pixels);
                arraycoord = image_pixels.index(maxpixel);

                # Add the coordinates of the maximum pixel to the centroid_x and centroid_y arrayself.
                # these aren't actual geometric centroids in this context.
                centroid_x.append(x0[arraycoord]);
                centroid_y.append(y0[arraycoord]);

                # Reset the arrays
                max_pixel_x             = [];       # max pixel x coordinates
                max_pixel_y             = [];       # max pixel y coordinates
                image_pixels            = [];       # storing image pixels from region coordinates.

    return centroid_x, centroid_y, len(centroid_coord)
