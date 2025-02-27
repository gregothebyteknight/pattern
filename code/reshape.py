
# IMPORT LIBRARIES
import numpy as np
import pandas as pd
import os
import skimage
import warnings
warnings.filterwarnings("ignore")

from skimage.io import imread
from pathlib import Path

# DOWNLOADING DATA AND LIST VARIABLES
masks = imread("../data/init/full_mask_final_segmentation_hwatershed_500.00_90%.tif")
print("Shape of the masks before:", masks.shape)
masks = np.swapaxes(np.swapaxes(masks, 0, 1), 1, 2) # swapping x and y axes, then z and y axes
print("Shape of the masks after:", masks.shape)

fold_list = [f.path for f in os.scandir("../data/init") if f.is_dir()]
n_chan = len(fold_list) # number of channels == n of cell markers
chan_list = [] # list of channel names

x_dim, y_dim, z_dim = 1655, 1604, 152 # set according to masks.shape
image = np.zeros((x_dim, y_dim, z_dim, n_chan), dtype = np.uint16) # with channels dim added

def reshape_image(fold_list):
    """
    Reshape the image data into a 4D numpy array (with channels added)
    @fold_list: list of folders containing the original image data
    """
    for k, folder in enumerate(fold_list):
        directory = Path(folder)
        file_list = sorted([f for f in directory.iterdir() if f.is_file()])
        image_temp = np.zeros((x_dim, y_dim, z_dim), dtype = np.uint16)
        
        for i, file in enumerate(file_list):
            image_temp[:, :, i] = imread(file, plugin = 'tifffile')
        image[:, :, :, k] = image_temp
    return image


def channel_order(fold_list, chan_list = chan_list):
    """
    Save the order of the ion channels in a csv file
    @fold_list: list of folders containing the original image data
    """
    for k in range(len(fold_list)):
        chan_list.append(fold_list[k][-5:])

# CELL INTENSITIES, COORDINATES AND CHANNEL ORDER
image = reshape_image(fold_list)
print("Shape of the image:", image.shape)
cell_int = skimage.measure.regionprops_table(masks, intensity_image = image,
                                                   properties = ['mean_intensity'])
cell_props = skimage.measure.regionprops_table(masks, 
                                               properties = ['area', 'centroid'])
cell_props['centroid-2'] = np.array(cell_props['centroid-2']) * 2 # width of slices (z axis) is 2 micrometers

channel_order(fold_list)
pd.DataFrame(cell_props).to_csv("../data/cell_coordinates.csv", header = ["area", "x", "y", "z"])
pd.DataFrame(cell_int).to_csv("../data/cell_intensities.csv", header = chan_list)