
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
masks = imread("../data/full_mask_final_segmentation_hwatershed_bg500_90%.tif")
masks = np.swapaxes(np.swapaxes(masks, 0, 1), 1, 2) # swapping x and y axes, then z and y axes
print("Shape of the masks:", masks.shape)

fold_list = [f.path for f in os.scandir("../data/SIMILARITY10_In115") if f.is_dir()]
n_chan = len(fold_list) # number of channels == n of cell markers
chan_list = [] # list of channel names

x_dim, y_dim, z_dim = 1079, 1095, 33
image = np.zeros((x_dim, y_dim, z_dim, n_chan)) # with channels dim added

def reshape_image(fold_list):
    """
    Reshape the image data into a 4D numpy array (with channels added)

    @fold_list: list of folders containing the original image data
    """
    for k in range(len(fold_list)):
        directory = Path(fold_list[k])
        file_list = sorted([f for f in directory.iterdir() if f.is_file()])

        image_temp = np.zeros((x_dim, y_dim, z_dim))
        for i in range(len(file_list)):
            image_temp[:, :, i] = imread(file_list[i], plugin = 'tifffile')
        image[:, :, :, k] = image_temp
    return image

def channel_order(fold_list, chan_list = chan_list):
    """
    Save the order of the ion channels in a csv file
    
    @fold_list: list of folders containing the original image data
    """
    for k in range(len(fold_list)):
        chan_list.append(fold_list[k][-5:])

    pd.DataFrame(chan_list).to_csv("../data/order_channels.csv", header = ["channel"])

# CELL INTENSITIES, COORDINATES AND CHANNEL ORDER
image = reshape_image(fold_list)
print("Shape of the image:", image.shape)
cell_int = skimage.measure.regionprops_table(np.squeeze(masks), intensity_image = image,
                                                   properties = ['mean_intensity'])
cell_props = skimage.measure.regionprops_table(np.squeeze(masks), 
                                               properties = ['area', 'centroid'])
channel_order(fold_list)
pd.DataFrame(cell_props).to_csv("../data/cell_coordinates.csv", header = ["area","x","y","z"])
pd.DataFrame(cell_int).to_csv("../data/cell_intensities.csv", header = chan_list)