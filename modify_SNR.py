
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


def inject():

    rebinned_flux = []
    rebinned_time = []
    
    for j in range(1,timestamps.shape[0]-1):
        bin_edges = [np.mean([timestamps[j-1], timestamps[j]]), np.mean([timestamps[j], timestamps[j+1]])]

        if(timestamps[j+1] - timestamps[j] > 1/48):
            bin_edges[1] = timestamps[j] + 1/96

        if(timestamps[j] - timestamps[j-1] > 1/48):
            bin_edges[0] = timestamps[j] - 1/96

        flux_bin = flux_total[np.where((time_kepler > bin_edges[0]) & (time_kepler < bin_edges[1]))]
        if(flux_bin.shape[0]>0):
            rebinned_flux.append(np.mean(flux_bin))
            rebinned_time.append(timestamps[j])

    rebinned_flux = np.array(rebinned_flux)
    rebinned_time = np.array(rebinned_time)

    injected_lc = np.copy(mod_keplerflux)

    for i in range(rebinned_time.shape[0]):
        injected_lc[np.where(timestamps == rebinned_time[i])] *= rebinned_flux[i]

    return injected_lc, rebinned_time


def mask_indices(time_stamps):

    mask_size = 3.0*transit_duration
    if(mask_size/per_bary > 0.5): mask_size = 0.5 * per_bary

    masks = []
    for i in range(num_epochs):
        masks.append(np.where((time_stamps > TT_kepler[i] - mask_size/2) & (time_stamps < TT_kepler[i] + mask_size/2))[0])

    mask_idxs = np.concatenate(masks)

    return mask_idxs


def modify_keplerflux():
    global trend

    f = open('/media/Datas/Exomoon_project/lightcurve_trends/lightcurve_trend_{}.csv'.format(system), 'r')
    trend = np.genfromtxt(f, delimiter=',', skip_header=1)[:,1]
    # plt.scatter(timestamps, kepler_flux, s=5)
    # plt.scatter(timestamps, trend, s=5)
    # plt.show()

    residual = kepler_flux - trend
    actual_noise = np.std(kepler_flux/trend)
    kepler_noise = np.median(flux_error/kepler_flux)

    avg_snr = np.mean([float(snr) for snr in line[5:5+num_moons]])

    req_factor =  (kepler_noise/actual_noise) * (avg_snr/target_SNR)

    mod_keplerflux = trend + residual*req_factor

    return mod_keplerflux




catalogue = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'r').read().split('\n')
target_SNR = 100.0
exclregion_file = open('/home/garvit/Downloads/Exomoon_project/excl_regions.txt', 'r').read().split('\n')[1:-1]


for n in range(len(catalogue)-3):
    systemno = n+3
    index= systemno -1
    line = catalogue[index].split('\t')
    system = line[0]
    print(n, system)

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

    mod_keplerflux = modify_keplerflux()

    system_file = line[0].split('_')
    date, res_nores, sim_no = system_file[0], system_file[1], system_file[2]
    f = open("/media/Datas/Exomoon_project/final_model_LCs/{}_{}_{}.txt".format(sim_no, res_nores, date), 'r')
    model_data = f.readlines()
    f.close()

    flux_planet_allmoons = np.ndarray((len(model_data) - 1))
    time = np.ndarray((len(model_data)-1))
    flux_moons = np.ndarray((num_moons, len(model_data)-1))

    for j in range(len(model_data)-1):
        time[j] = float(model_data[j+1].split('\t')[0])
        flux_planet_allmoons[j] = float(model_data[j+1].split('\t')[1])
        for k in range(num_moons):
            flux_moons[k][j] = float(model_data[j+1].split('\t')[2+k])

    best_offset = float(line[16])
    time_kepler = time - time[0] + kepler_starttime + best_offset
    model_time_interval = time[-1] - time[0]

    flux_planet_allmoons_temp = flux_planet_allmoons.copy()
    flux_moons_temp = flux_moons.copy()
    time_temp = time.copy()
    if(data_time_interval < model_time_interval + best_offset):
        flux_planet_allmoons_temp = flux_planet_allmoons[time_kepler < kepler_endtime]
        flux_moons_temp = flux_moons[:, time_kepler < kepler_endtime]
        time_temp = time[time_kepler < kepler_endtime]
        time_kepler = time_kepler[time_kepler < kepler_endtime]

    flux_total = np.copy(flux_planet_allmoons_temp) 
    for j in range(num_moons): flux_total+= flux_moons_temp[j]
    flux_total-=(num_moons)
    #plt.scatter(time_kepler, flux_total, s=0.5)
    # plt.show()

    TTs = model_data[0].split('\t')[:-1]
    num_epochs = len(TTs)
    TT_planet_allmoons = np.ndarray((num_epochs))
    for j in range(num_epochs):
        TT_planet_allmoons[j] = TTs[j]

    TT_kepler = TT_planet_allmoons - time[0] + kepler_starttime + best_offset
    TT_kepler = TT_kepler[TT_kepler < kepler_endtime]
    num_epochs = TT_kepler.shape[0]
    # print(TT_kepler)

    injected_flux, rebinned_time = inject()
    # injected_flux = kepler_flux.copy()
    print("injection complete\n")

    transit_duration = float(line[17])
    per_bary = float(line[3])
    
    trends = [0 for j in range(num_quarters)]
    detrended_fluxes = [0 for j in range(num_quarters)]
    for j in range(num_quarters):
        mask_idxs = mask_indices(timestamps[quarters[j]])
        trends[j], detrended_fluxes[j], errors_detrend = cofiam_detrend(times=timestamps[quarters[j]], fluxes=injected_flux[quarters[j]], errors=flux_error[quarters[j]], mask_idxs=mask_idxs, max_degree=30)

    trend = np.concatenate(trends)
    detrended_flux = np.concatenate(detrended_fluxes)
    print("detrending complete\n")
    # detrended_flux = injected_flux/trend
    # plt.scatter(timestamps, injected_flux, s=5)
    # plt.scatter(timestamps, trend, s=5)
    # plt.show()

    line2= exclregion_file[index-2].split(',')
    num_regions = int((len(line2)-3)/2)
    if(num_regions>0):
        excl_regions=[]
        for i in range(num_regions):
            excl_regions.append([float(line2[3+2*i]),float(line2[3+2*i+1])])

    timestamps_excl = timestamps.copy()
    injected_flux_excl = injected_flux.copy()
    trend_excl = trend.copy()
    detrended_flux_excl = detrended_flux.copy()
    rebinned_time_excl = rebinned_time.copy()

    for i in range(num_regions):
        start = excl_regions[i][0]
        end = excl_regions[i][1]
        indices = np.where((timestamps_excl < start) | (timestamps_excl > end))
        injected_flux_excl = injected_flux_excl[indices]
        trend_excl = trend_excl[indices]
        detrended_flux_excl = detrended_flux_excl[indices]
        timestamps_excl = timestamps_excl[indices]
        rebinned_time_excl = rebinned_time_excl[np.where((rebinned_time_excl < start) | (rebinned_time_excl > end))]

    final_lc = np.ndarray((rebinned_time_excl.shape[0]))
    for j in range(rebinned_time_excl.shape[0]):
        final_lc[j] = detrended_flux_excl[np.where(timestamps_excl == rebinned_time_excl[j])]
    plt.scatter(rebinned_time_excl, final_lc, s=5, color='black')
    plt.show()

    rebinned_time = rebinned_time_excl

    
    print("Writing the injected light curve for this system in a file\n-----------------------------------------------------\n")
    f = open('/home/garvit/Downloads/Exomoon_project/final_lightcurves_SNR/final_lightcurve_{}_{}_{}.csv'.format(date, res_nores, sim_no), 'w')
    f.write("System architecture file: {} {} {}\n".format(date, res_nores, sim_no))
    f.write("Kepler light curve file: {}\n".format(kepler_lightcurve_file))
    writer1 = csv.writer(f, delimiter=',')
    writer1.writerows(zip(rebinned_time, final_lc))
    f.close()


    


    
