from __future__ import division
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [18, 12]
import numpy as np

import pandoramoon as pandora

import pickle
import math
import random
import os
import csv
import sys

from statsmodels.stats import stattools

#MoonPy's cofiam code ---------------------------------------------------
import numpy as np
try:
	from numba.core.decorators import jit
except:
	from numba import jit 


def max_order(times, duration, baseline=0, kmaximum=30):
	if baseline == 0:
		baseline = np.nanmax(times) - np.nanmin(times)
	assert duration > 0
	kmax = int((2*baseline) / (12*duration))
	if kmax > kmaximum:
		kmax = kmaximum
	if kmax == 0:
		kmax = 1
	return kmax 


def DurbinWatson(residuals):
	residual_terms = []
	for nres, res in enumerate(residuals):
		try:
			residual_terms.append(residuals[nres+1] - residuals[nres])
		except:
			pass
	residual_terms = np.array(residual_terms)
	numerator = np.nansum(residual_terms**2)
	denominator = np.nansum(residuals**2)
	assert denominator != 0.
	return numerator / denominator

@jit(debug=True, fastmath=True, nopython=True, cache=True)
def cofiam_matrix_gen(times, degree):
	baseline = np.nanmax(times) - np.nanmin(times)
	assert baseline > 0
	rows = len(times)
	cols = 2 * (degree+1)
	X_matrix = np.ones(shape=(rows,cols))
	for x in range(rows):
		for y in range(1, int(cols/2)):
			sinarg = (2*np.pi*times[x] * y) / baseline
			X_matrix[x,y*2] = np.sin(sinarg)
			X_matrix[x,y*2 + 1] = np.cos(sinarg)
		X_matrix[x,1] = times[x]
	return X_matrix 


def cofiam_matrix_coeffs(times, fluxes, degree, solve=True):
	assert len(times) > 0
	Xmat = cofiam_matrix_gen(times, degree)
	beta_coefs = np.linalg.lstsq(Xmat, fluxes, rcond=None)[0]
	return Xmat, beta_coefs

def cofiam_function(times, fluxes, degree, solve=True, cofiam_coefficients=None):
	input_times = times.astype('f8')
	input_fluxes = fluxes.astype('f8')
	
	if solve == True:
		cofiam_matrix, cofiam_coefficients = cofiam_matrix_coeffs(input_times, input_fluxes, degree, solve=True)
	
	elif solve == False:
		assert type(cofiam_coefficients) != type(None)
		cofiam_matrix = cofiam_matrix_gen(input_times, degree)

	model = np.matmul(cofiam_matrix, cofiam_coefficients)
	return model, cofiam_coefficients


def cofiam_iterative(times, fluxes, max_degree=30, min_degree=1):
    vals_to_min = []
    degs_to_try = np.arange(min_degree,max_degree+1,1)
    DWstats = []

    for deg in degs_to_try:
        output_model, cofiam_coefficients = cofiam_function(times=times, fluxes=fluxes, degree=deg, solve=True)

        residuals = fluxes - output_model

        DWstat = DurbinWatson(residuals)
        DWstats.append(DWstat)

        val_to_minimize = (DWstat - 2)**2
        vals_to_min.append(val_to_minimize)

    best_degree = degs_to_try[np.argmin(np.array(vals_to_min))]
    best_DW = DWstats[np.argmin(np.array(vals_to_min))]

    best_model, best_coefficients = cofiam_function(times=times, fluxes=fluxes, degree=best_degree, solve=True)

    return best_model, best_degree, best_coefficients, best_DW, max_degree 




def cofiam_detrend(times, fluxes, errors, mask_idxs=None, max_degree=30):
	if type(mask_idxs) != type(None):
		# print('len(mask_idxs) = ', len(mask_idxs))
		if len(mask_idxs) > 0:
			oot_times, oot_fluxes, oot_errors = np.delete(times, mask_idxs), np.delete(fluxes, mask_idxs), np.delete(errors, mask_idxs)
		elif len(mask_idxs) == 0:
			oot_times, oot_fluxes, oot_errors = times, fluxes, errors
	
	# print('len(times), len(oot_times) = ', len(times), len(oot_times))

	__, best_degree, best_coefficients, best_DW, __ = cofiam_iterative(times=np.array(oot_times, dtype=np.float64), fluxes=np.array(oot_fluxes, dtype=np.float64), max_degree=int(max_degree))
	# print('best_degree: ', best_degree)
	# print('Durbin-Watson statistic: ', best_DW)
	# print(' ')
	# print(' ')

	
	best_model = cofiam_function(times=times, fluxes=fluxes, degree=best_degree, solve=False, cofiam_coefficients=best_coefficients)[0]

	fluxes_detrend = fluxes / best_model
	errors_detrend = errors / fluxes
	return best_model, fluxes_detrend, errors_detrend #### best model is the COFIAM

#-----------------------------------------------------------------------

def find_nearest(array, value):
    idx = np.argmin(np.abs(array - value))
    return idx





catalogue = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'r').read().split('\n')

for n in range(len(catalogue)-3):
    systemno = n+3
    i= systemno -1
    line = catalogue[i].split('\t')

    num_moons = int(line[2])

    kepler_lightcurve_file = line[1]
    f = open("/home/garvit/Downloads/Exomoon_project/good_kepler_data/{}".format(kepler_lightcurve_file), 'r')
    kepler_data = f.readlines()
    f.close()
    kepler_starttime = float(kepler_data[1].split('\t')[0])
    kepler_endtime = float(kepler_data[-1].split('\t')[0])
    data_time_interval = kepler_endtime - kepler_starttime

    N = len(kepler_data) -1
    kepler_flux = np.ndarray((N))
    timestamps = np.ndarray((N))
    flux_error = np.ndarray((N))
    quarter_arr = np.ndarray((N))

    for j in range(N):
        timestamps[j] = float(kepler_data[j+1].split('\t')[0])
        kepler_flux[j] = float(kepler_data[j+1].split('\t')[1])
        flux_error[j] = float(kepler_data[j+1].split('\t')[2])
        quarter_arr[j] = int(kepler_data[j+1].split('\t')[4])

    quarter_list = list(dict.fromkeys(quarter_arr))
    num_quarters = len(quarter_list)
    quarters = []
    for j in quarter_list:
        idxs = np.where(quarter_arr == j)[0]
        if(idxs.shape[0]==1): 
            num_quarters -= 1
            continue
        quarters.append(idxs)

    system_file = line[0].split('_')
    print(line[0])
    date, res_nores, sim_no = system_file[0], system_file[1], system_file[2]
    
    best_offset = float(line[16])
    

    
    transit_duration = float(line[17])
    per_bary = float(line[3])
    
    trends = [0 for j in range(num_quarters)]
    detrended_fluxes = [0 for j in range(num_quarters)]
    for j in range(num_quarters):
        mask_idxs = []
        trends[j], detrended_fluxes[j], errors_detrend = cofiam_detrend(times=timestamps[quarters[j]], fluxes=kepler_flux[quarters[j]], errors=flux_error[quarters[j]], mask_idxs=mask_idxs, max_degree=30)

    trend = np.concatenate(trends)

    # plt.scatter(timestamps, kepler_flux, s=5)
    # plt.scatter(timestamps, trend, s=5)
    # plt.show()

    f = open('/media/Datas/Exomoon_project/lightcurve_trends/lightcurve_trend_{}_{}_{}.csv'.format(date, res_nores, sim_no), 'w')
    f.write("Kepler light curve file: {}\n".format(kepler_lightcurve_file))
    writer1 = csv.writer(f, delimiter=',')
    writer1.writerows(zip(timestamps, trend))
    f.close()