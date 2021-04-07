'''

MONTE-CARLO ANALYSIS OF TDOA LOCALISATION
COMPARING MLE, BLUE WITH CRLB

AUTHOR: ABIJITH J KAMATH
abijithj@iisc.ac.in

'''

# %% IMPORT PACKAGES
import os
import numpy as np

from tqdm import tqdm
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

data = loadmat('TDOA_data.mat')

anchor_location = data['anchor_location'].astype(np.float)
target_location = data['target_location'].astype(np.float)
true_distances = data['true_distances']

noisy_distances = data['noisy_distances']
var_vec = data['sigma2']


# %%

_, num_iter, var_len = noisy_distances.shape
init_pos = np.array([0.0,0.0])[:,np.newaxis]

min_x_var = np.zeros(var_len)
min_y_var = np.zeros(var_len)

mle_x_mse = np.zeros(var_len)
mle_y_mse = np.zeros(var_len)

blue_x_mse = np.zeros(var_len)
blue_y_mse = np.zeros(var_len)

for noise_iter in tqdm(range(var_len)):

    noise = var_vec[noise_iter]

    # Compute CRLB
    fim = utils.crlb_tdoa(target_location, anchor_location, noise)
    cov_est = np.linalg.inv(fim)
    min_x_var[noise_iter] = cov_est[0,0]
    min_y_var[noise_iter] = cov_est[1,1]

    # Compute error in MLE and BLUE estimates
    for avg_iter in range(num_iter):
        mle_est = utils.mle_tdoa(noisy_distances[:,avg_iter,noise_iter],
            anchor_location, init_pos, var=noise, step_size=noise/2.)

        blue_est = utils.blue_tdoa(noisy_distances[:,avg_iter,noise_iter],
            anchor_location, var=noise)

        mle_x_mse[noise_iter] += (mle_est[0]-target_location[0])**2
        mle_y_mse[noise_iter] += (mle_est[1]-target_location[1])**2

        blue_x_mse[noise_iter] += (blue_est[0]-target_location[0])**2
        blue_y_mse[noise_iter] += (blue_est[1]-target_location[1])**2

    mle_x_mse[noise_iter] /= num_iter
    mle_y_mse[noise_iter] /= num_iter

    blue_x_mse[noise_iter] /= num_iter
    blue_y_mse[noise_iter] /= num_iter

# %% PLOT MSE

os.makedirs('./results', exist_ok=True)
path = './results/'

# Plotting for x-coordinate
plt.figure(figsize=(12,6))
ax = plt.gca()

utils.plot_curve(var_vec, mle_x_mse, ax=ax, plot_colour='red',
    legend_label=r'MLE', line_width=4)
utils.plot_curve(var_vec, blue_x_mse, ax=ax, plot_colour='blue',
    legend_label=r'BLUE', line_width=4)
utils.plot_curve(var_vec, min_x_var, ax=ax, plot_colour='green',
    legend_label=r'CRLB', line_width=4, xaxis_label=r'$\sigma^2$',
    title_text=r'MSE in $x$-Coordinate', xlimits=[0,10], ylimits=[0,10],
    show=True, save=path+'crlb_x')

# Plotting for y-coordinate
plt.figure(figsize=(12,6))
ax = plt.gca()

utils.plot_curve(var_vec, mle_y_mse, ax=ax, plot_colour='red',
    legend_label=r'MLE', line_width=4)
utils.plot_curve(var_vec, blue_y_mse, ax=ax, plot_colour='blue',
    legend_label=r'BLUE', line_width=4)
utils.plot_curve(var_vec, min_y_var, ax=ax, plot_colour='green',
    legend_label=r'CRLB', line_width=4, xaxis_label=r'$\sigma^2$',
    title_text=r'MSE in $y$-Coordinate', xlimits=[0,10], ylimits=[0,10],
    show=True, save=path+'crlb_y')