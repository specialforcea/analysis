# -*- coding: utf-8 -*-
"""
Created on Fri Dec 09 10:37:47 2016

@author: rubidium
"""

from __future__ import division
import numpy as np

def fluctuations_from_mean(measurement_array):
    measurement_array = np.array(measurement_array)
    mean = np.average(measurement_array, axis=0)
    return mean, np.array([(shot-mean) for shot in measurement_array])

def covariance(measurement_array, norm=False):
    measurement_array = np.array(measurement_array)
    y_bar, delta_y = fluctuations_from_mean(measurement_array)
    N, M = measurement_array.shape[0], measurement_array.shape[1]
    xlim = np.int(np.amax(measurement_array.shape)/2)
    x_ax = np.linspace(-xlim, xlim, np.amax(measurement_array.shape))
    sigma_j = np.array([np.std(delta_y[j]) for j in range(N)])
    if norm:
        product = np.array([delta_y[j, n]*delta_y[j, m]/sigma_j[j]**2 for j in range(N) for n in range(M) for m in range(M)])
    else:
        product = np.array([delta_y[j, n]*delta_y[j, m] for j in range(N) for n in range(M) for m in range(M)])
    product = product.reshape((N, M, M))
    print np.mean(product)
    return x_ax, np.average(product, axis=0) 
    
# Test
if False:    
    x_test = [[0, 1, 4, 2, 6, 5], [0.2, 1, 0.5, 0, 6, 5], [0.4, 1.6, 5, 4.1, 6.1, 5.2], [0, 3, 5, 2, 6.4, 5], [0.11, 12, 2, 4, 6, -5]]
    x_mean, delta_x = fluctuations_from_mean(x_test)
    cov_axis, cov_xdata = covariance(x_test, norm=True)
    
    from matplotlib import pyplot as plt
    fig = plt.figure(frameon=False)
    plt.pcolormesh(cov_axis, cov_axis, cov_xdata, cmap='viridis')
    plt.show()