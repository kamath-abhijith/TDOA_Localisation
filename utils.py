'''

UTILITY FUNCTIONS FOR TDOA LOACALISATION

AUTHOR: ABIJITH J KAMATH
abijithj@iisc.ac.in

'''

# %% IMPORT LIBRARIES
import numpy as np

from matplotlib import pyplot as plt

# %% PLOTTING

def plot_sensors(target, anchors, image=None, ax=None, xaxis_label=None,
    yaxis_label=None, title_text=None, marker='+', markersize=4,
    target_colour='red', legend_label=None, legend_show=True,
    legend_loc='upper left', show=True, save=None):
    ''' Plots target and anchors on image '''
    if ax is None:
        fig = plt.figure(figsize=(12,6))
        ax = plt.gca()

    plt.plot(target[0,0], target[1,0], marker, color=target_colour,
        marker=marker, markersize=markersize, label=legend_label)
    plt.plot(anchors[0,:], anchors[1,:], marker, color='black',
        marker=marker, markersize=markersize)
    
    if image is not None:
        plt.imshow(image)

    if legend_label and legend_show:
        plt.legend(loc=legend_loc, frameon=True, framealpha=0.8, facecolor='white')

    if show:
        plt.xlabel(xaxis_label)
        plt.ylabel(yaxis_label)
        plt.title(title_text)
        plt.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
        if save:
            plt.savefig(save + '.pdf', format='pdf')
        plt.show()

    return

def plot_curve(x, y, ax=None, plot_colour='blue', xaxis_label=None,
    yaxis_label=None, title_text=None, legend_label=None, legend_show=True,
    legend_loc='upper left', line_style='-', line_width=None,
    show=False, xlimits=[0,1], ylimits=[-1,1], save=None):
    '''
    Plots signal with abscissa in x and ordinates in y 

    '''
    if ax is None:
        fig = plt.figure(figsize=(12,6))
        ax = plt.gca()

    plt.plot(x, y, linestyle=line_style, linewidth=line_width, color=plot_colour, label=legend_label)
    if legend_label and legend_show:
        plt.legend(loc=legend_loc, frameon=True, framealpha=0.8, facecolor='white')
    plt.xlabel(xaxis_label)
    plt.ylabel(yaxis_label)

    if show:
        plt.xlim(xlimits)
        plt.ylim(ylimits)
        plt.title(title_text)
        if save:
            plt.savefig(save + '.pdf', format='pdf')
        plt.show()

    return

# %% OPERATORS

def covariance_mtx(var):
    ''' Returns covariance matrix of range differences '''
    return var * np.array([[2,-1,0], [-1,2,-1], [0,-1,2]])

def target_anchor_distance(x, anchors):
    ''' Returns distance between target and anchors '''
    return np.linalg.norm(x-anchors, axis=0)

def dist_diff(distances):
    ''' Returns distance differences '''
    return -1.0 * np.diff(distances)

def dist_derivative(x, anchors, distances):
    ''' Returns derivative of distance differences '''
    num_anchors = anchors.shape[1]
    return (((x - anchors)/distances)[:,:num_anchors-1] - 
        ((x - anchors)/distances)[:,1:]).T

# %% ACCURACY ANALYSIS

def crlb_tdoa(target, anchors, var):
    '''
    Returns the Fisher information matrix of 2D TDOA

    :param target: 2D position of the target
    :param anchors: 2D position of the anchors
    :param var: Variance of the noise in range measurements

    :return fim: Fisher information matrix

    '''

    num_anchors = anchors.shape[1]

    # Noise statistics
    covW = covariance_mtx(var)
    presW = np.linalg.inv(covW)

    # Distance vector
    target_anchor_dist = target_anchor_distance(target, anchors)

    del_dist_vec = dist_derivative(target, anchors, target_anchor_dist)
    
    # Fisher information matrix
    fim = (del_dist_vec.T.dot(presW)).dot(del_dist_vec)
        
    return fim

# %% SOLVERS

def mle_tdoa(distances, anchors, init_pos, var=1., step_size=0.1, max_iter=100):
    '''
    Returns maximum-likelihood estimate (MLE) for TDOA localisation.
    Gradient descent with fixed step-size is used as the solver.

    :param distances: distance measurements
    :param anchors: 2D position of anchors
    :param init_pos: Initialisation for gradient descent
    :param step_size: Step size of gradient descent optimiser
    :param max_iter: Maximum iterations before exit
    :param var: variance of noise in distance measurements

    :returns: mle estimate for position

    '''

    covN = covariance_mtx(var)
    presN = np.linalg.inv(covN)
    
    range_diff = dist_diff(distances)
    mle_estimate = init_pos

    for _ in range(max_iter):
        dx = target_anchor_distance(mle_estimate, anchors)
        dx_vec = dist_diff(dx)
        del_dx = dist_derivative(mle_estimate, anchors, dx)
        mle_estimate += step_size * ((del_dx.T).dot(presN)).dot(range_diff -
            dx_vec)[:,np.newaxis]

    return mle_estimate

def blue_tdoa(distances, anchors, var=1.):
    '''
    Returns best linear unbiased estimate (BLUE) for TDOA localisation.

    :param distances: distance measurements
    :param anchors: 2D position of anchors
    :param var: variance of noise in distance measurements

    :returns: blue estimate for position

    '''

    num_anchors = anchors.shape[1]

    cov_mtx = 2.0*var**2 * np.array([[4, 1, 1, 1, 1, 0],
                                    [1, 4, 1, 1, 0, 1],
                                    [1, 1, 4, 0, 1, 1],
                                    [1, 1, 0, 4, 1, 1],
                                    [1, 0, 1, 1, 4, 1],
                                    [0, 1, 1, 1, 1, 4]], dtype=np.float)
    pres_mtx = np.linalg.inv(cov_mtx)
    r = np.zeros(6)
    gamma_vec = np.zeros(6)
    H_mtx = np.zeros([6, 5])

    idx = 0
    for i in range(num_anchors):
        for j in range(i+1, num_anchors):
            r[idx] = distances[i]-distances[j]
            gamma_vec[idx] = (r[idx])**2 - np.linalg.norm(anchors[:,i])**2 + np.linalg.norm(anchors[:,j])**2 - 2*var
            idx = idx+1

    idx = 0
    for anchor_idx in range(num_anchors - 1):
        lst = list(range((anchor_idx), num_anchors - 1))
        for pos in lst:
            H_mtx[idx,0] = -2*(anchors[0,anchor_idx] - anchors[0,(pos+1)])
            H_mtx[idx,1] = -2*(anchors[1,anchor_idx] - anchors[1,(pos+1)])

            H_mtx[idx,(pos+2)] = -2*r[idx]
            idx = idx+1

    den = np.linalg.inv((((H_mtx.T).dot(pres_mtx)).dot(H_mtx)))
    num = ((H_mtx.T).dot(pres_mtx)).dot(gamma_vec-2*var)
    blue_est = den.dot(num)

    return (blue_est[:2])[:,np.newaxis]