'''

TIME DIFFERENCE OF ARRIVAL (TDOA) LOCALISATION
USING MAXIMUM LIKELIHOOD ESTIMATION

AUTHOR: ABIJITH J KAMATH
abijithj@iisc.ac.in

'''

# %% IMPORT PACKAGES
import os
import numpy as np

from skimage.io import imread
from scipy.io import loadmat

from matplotlib import pyplot as plt
from matplotlib import style
from matplotlib import rcParams

import utils

# %% PLOT SETTINGS

plt.style.use(['science','ieee'])

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["cm"],
    "mathtext.fontset": "cm",
    "font.size": 11})

# %% LOAD DATA

img = imread('mapimage.jpeg')

data = loadmat('TDOA_data.mat')

anchor_location = data['anchor_location'].astype(np.float)
target_location = data['target_location'].astype(np.float)
true_distances = data['true_distances']

noisy_distances = data['noisy_distances']
var_vec = data['sigma2']

# %% ESTIMATE USING MLE

init_pos = np.array([0.0,0.0])[:,np.newaxis]
mle_pos = utils.mle_tdoa(noisy_distances[:,0,0], anchor_location, init_pos,
    var=0.001, step_size=.0005, max_iter=100)

# %% PLOT LOCATIONS

os.makedirs('./results', exist_ok=True)
path = './results/'

utils.plot_sensors(target=mle_pos, anchors=anchor_location, image=img,
    xaxis_label=r'$x$', yaxis_label=r'$y$', title_text=r'TDOA Localisation using MLE',
    marker='+', markersize=11, target_colour='red', legend_label=r'MLE',
    show=True, save=path+'mle')
# %%
