# IMPORT LIBRARIES
import os
import torch
import warnings
warnings.filterwarnings("ignore")

import numpy as np

from STalign import STalign
from scipy.stats import trim_mean 
from contextlib import redirect_stdout # to supress output

# INITIALIZE VARIABLES
map_slice = {0: 2, 1: 7, 2: 14, 3: 20, 4: 25, 5: 29, 6: 34, 7: 39, 8: 44, 
             9: 49, 10: 50, 11: 51, 12: 52, 13: 54, 14: 59, 15: 64, 16: 69, 
             17: 74, 18: 78, 19: 84, 20: 86, 21: 91, 22: 97, 23: 102, 24: 106}

def affine(deg, df_s, df_t = None):
    """
    Performs affine transformation on the df_s
    Applies both rotation and translational matrices
    Additionally - scaling in respect of df_t
    @def(int): degree of rotation (clockwise)
    @df_s(pd.DataFrame): source slice dataset
    @df_t(pd.DataFrame): target slice dataset
    """
    rad = np.deg2rad(-deg)
    x_s, y_s = df_s['x'].to_numpy(), df_s['y'].to_numpy()

    # rotation matrix
    l_mat = np.array([[np.cos(rad), -np.sin(rad)],
                  [np.sin(rad), np.cos(rad)]])
    x_rot, y_rot = l_mat @ np.vstack((x_s, y_s))

    if df_t is not None:
        x_t, y_t = df_t['x'].to_numpy(), df_t['y'].to_numpy()

        # applying scaling
        x_scl = np.std(x_t) / np.std(x_rot) * x_rot
        y_scl = np.std(y_t) / np.std(y_rot) * y_rot
        
        # translation matrix
        t_mat = np.array([
            trim_mean(x_t, proportiontocut = 0.1) - trim_mean(x_scl, proportiontocut = 0.1),
            trim_mean(y_t, proportiontocut = 0.1) - trim_mean(y_scl, proportiontocut = 0.1)
    ])

    else: 
        t_mat = [0, 0]
        x_scl, y_scl = x_rot, y_rot

    df_s['x'] = x_scl + t_mat[0] 
    df_s['y'] = y_scl + t_mat[1] 

    return df_s

def coord_align(df_s, df_t):
    """
    Shift coordinates of sample dataset towards target
    using STalign approach
    @df_s(pd.DataFrame): source dataframe
    @df_t(pd.DataFrame): target dataframe
    """
    pts_s = df_s[['x', 'y']].to_numpy()
    pts_t = df_t[['x', 'y']].to_numpy()

    # RASTERIZE
    with open(os.devnull, 'w') as fnull, redirect_stdout(fnull):
        x_s_rast, y_s_rast, s_tensor, _ = STalign.rasterize(pts_s[:, 0], pts_s[:, 1],
                                                            dx = 30, blur = 1.5)
        x_t_rast, y_t_rast, t_tensor, _ = STalign.rasterize(pts_t[:, 0], pts_t[:, 1],
                                                            dx = 30, blur = 1.5)

    # SOLVE MAPPING
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    torch.set_default_device(device) # for point mapping

    params = {'niter': 10000, 'diffeo_start': 10000, 'device': device,
              'epL': 2e-08, 'epV': 0}

    out = STalign.LDDMM([y_s_rast, x_s_rast], s_tensor, [y_t_rast, x_t_rast], t_tensor, **params)

    a = out['A']
    v = out['v']
    xv = out['xv']

    # APPLY MAPPING
    points = STalign.transform_points_source_to_target(xv, v, a, pts_s[:, [1, 0]])

    # switch tensor from cuda to cpu for plotting with numpy
    if points.is_cuda:
        points = points.cpu()

    df_s['x'] = points[:, 1]
    df_s['y'] = points[:, 0]

    return df_s

def scatter(df, ax, title):
    """
    Maps plot to the panel
    @df(pd.DataFrame): dataframe
    @ax(matplotlib): matplotlib axes object
    @title(str): coordinates of axis on the panel 
    """
    ax.scatter(x = df['x'], y = df['y'], c = df['area'], s = 0.5)
    ax.set_title(title, fontsize = 8)