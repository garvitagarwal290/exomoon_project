from __future__ import division
import matplotlib.pyplot as plt
#plt.rcParams['figure.figsize'] = [6, 4]
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


def get_modelLCs(sim_params, data_time_interval):

    params = pandora.model_params()

    #Star parameters
    params.R_star = sim_params['Planet']['Rstar'] # [m]	 
    params.u1 = 2*sim_params['Planet']['q2']*np.sqrt(sim_params['Planet']['q1'])
    params.u2 = -2*(sim_params['Planet']['q2'] - 0.5) * np.sqrt(sim_params['Planet']['q1'])

    #Planet parameters
    params.per_bary = sim_params['Planet']['Pp']  # [days]
    params.a_bary = sim_params['Planet']['a']/params.R_star  # [R_star]
    params.r_planet = sim_params['Planet']['RpRstar'] # [R_star]
    params.b_bary = sim_params['Planet']['impact']   # [0..1]
    params.t0_bary = 10  # [days]
    params.t0_bary_offset = 0.0 #0.01  # [days]
    params.M_planet =  sim_params['Planet']['m'] # [kg]

    w_bary = sim_params['Planet']['pomega']
    params.w_bary = math.degrees(w_bary) if w_bary!=None else random.random()*2*np.pi  # [deg]
    ecc_bary= sim_params['Planet']['e']
    params.ecc_bary = ecc_bary if ecc_bary!=None else 0  # [0..1]  

    transit_duration = params.per_bary/np.pi * np.arcsin(np.sqrt((params.R_star + params.r_planet*params.R_star)**2 - (params.b_bary*params.R_star)**2)/(params.a_bary*params.R_star))
    v_planet =   2*np.pi * sim_params['Planet']['a']/ params.per_bary  #(2*params.R_star * (1 - params.b_bary)) / transit_duration
    t_hill = sim_params['Planet']['RHill'] / v_planet

    #Other parameters
    epoch_duration = int(np.ceil(t_hill * 2 + transit_duration))
    params.epoch_duration = epoch_duration # [days]

    params.cadences_per_day = 30000 # [int]         #Using short cadence to generate model light curves and caluclating TTVs
    params.epoch_distance = params.per_bary   # [days]
    params.supersampling_factor = 1  # [int]
    params.occult_small_threshold = 0.01  # [0..1]
    params.hill_sphere_threshold = 1.1

    num_moons = len(sim_params.keys()) - 1

    moon_transit_depths = np.ndarray((num_moons))
    # for j in range(num_moons):
    #     moon_key = list(sim_params.keys())[j+1]
    #     moon_transit_depths[j] = (sim_params[moon_key]['r']/params.R_star)**2

    flux_planet_allmoons = np.ndarray((len(model_data) - 1))
    time = np.ndarray((len(model_data)-1))
    flux_moons = np.ndarray((num_moons, len(model_data)-1))

    for j in range(len(model_data)-1):
        time[j] = float(model_data[j+1].split('\t')[0])
        flux_planet_allmoons[j] = float(model_data[j+1].split('\t')[1])
        for i in range(num_moons):
            flux_moons[i][j] = float(model_data[j+1].split('\t')[2+i])

    TTs = model_data[0].split('\t')[:-1]
    params.epochs = len(TTs)
    TT_planet_allmoons = np.ndarray((params.epochs))
    for j in range(params.epochs):
        TT_planet_allmoons[j] = TTs[j]

    return time, flux_planet_allmoons, flux_moons, params, num_moons, TT_planet_allmoons, transit_duration, moon_transit_depths


def moon_transits_details():

    noise_arr = flux_error/kepler_flux
    noise = np.mean(noise_arr)
    
    N_transit = int(np.round(transit_duration * 48))

    SNR = np.ndarray((num_moons))
    abs_best_SNR = np.ndarray((num_moons))
    moon_transits10 = np.ndarray((num_moons, params.epochs))  
    N = int(params.epoch_duration * params.cadences_per_day)
    for i in range(num_moons):

        best_transit_depth = 0.0
        for j in range(params.epochs):

            if((flux_moons_temp[i][j*N : (j+1)*N]==1.0).all() or planet_transits10[j]==0): moon_transits10[i][j]=0
            else: 
                moon_transits10[i][j]=1
                transit_depth = 1 - np.min(flux_moons_temp[i][j*N: (j+1)*N])
                if(transit_depth > best_transit_depth): best_transit_depth = transit_depth

        #SNR[i] = moon_transit_depths[i]* np.sqrt(N_transit) / noise
        SNR[i] = best_transit_depth * np.sqrt(N_transit) / noise
        #print(transit_depth*1e6, end=' ')

        abs_best_transit_depth = 1.0 - np.min(flux_moons[i])
        abs_best_SNR[i] = abs_best_transit_depth * np.sqrt(N_transit) / noise

    return SNR, moon_transits10, abs_best_SNR


def best_planet_transits():

    tau0_offsets = np.linspace(0, params.per_bary - TT_planet_allmoons[0], 100)
    ngap_transits_list=np.zeros((100))

    if(data_gaps.shape[0]!=0):
        for i in range(100):
            ntt = TT_planet_allmoons - time[0] + kepler_starttime + tau0_offsets[i]
            ngap_transits = 0

            for tt in ntt:
                if (np.any((tt > data_gaps.T[0]) * (tt < data_gaps.T[1])) or tt > kepler_endtime):
                    ngap_transits += 1 

            ngap_transits_list[i] = ngap_transits

        best_offset = tau0_offsets[np.argmin(ngap_transits_list)]

    else:
        best_offset = 0.0

    TT_planet_allmoons_kepler = TT_planet_allmoons - time[0] + kepler_starttime 
    params.epochs = (TT_planet_allmoons_kepler + best_offset < kepler_endtime).sum()
    planet_transits10 = np.ones(params.epochs)

    if(data_gaps.shape[0]!=0):
        for j in range(params.epochs):
            if np.any((TT_planet_allmoons_kepler[j] + best_offset > data_gaps.T[0]) & (TT_planet_allmoons_kepler[j] + best_offset < data_gaps.T[1])): planet_transits10[j] = 0
            else: planet_transits10[j] = 1

    return planet_transits10, best_offset



def kepler_data_gaps(timestamps):
    
    data_gaps=[]
    min_gap_size = .25 * transit_duration #days
    for i in range(timestamps.shape[0]-1):
        if(timestamps[i+1] - timestamps[i] > min_gap_size):
            data_gaps.append((timestamps[i], timestamps[i+1]))
    return data_gaps


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

    injected_lc = np.copy(kepler_flux)

    for i in range(rebinned_time.shape[0]):
        injected_lc[np.where(timestamps == rebinned_time[i])] *= rebinned_flux[i]

    return injected_lc, rebinned_time


def mask_indices(time_stamps):

    mask_size = 8*transit_duration
    if(mask_size/params.per_bary > 0.5): mask_size = 0.5 * params.per_bary

    masks = []
    for i in range(params.epochs):
        masks.append(np.where((time_stamps > TT_kepler[i] - mask_size/2) & (time_stamps < TT_kepler[i] + mask_size/2))[0])

    mask_idxs = np.concatenate(masks)

    return mask_idxs



h = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'w')
h.write('system_file\tkepler_file\tnum_moons\tPlanet_period(days)\tP_Transits_present/absent\tSNR1\tSNR2\tSNR3\tSNR4\tSNR5\tM_Transits_present/absent1\tM_Transits_present/absent2\tM_Transits_present/absent3\tM_Transits_present/absent4\tM_Transits_present/absent5\tDW_statistic\tBest_offset\tPlanet_Transit_duration\tt0_bary\n\n')

sim_nos = np.ndarray((100))
dates = []
res_nores_list = []
f = open('/home/garvit/Downloads/Exomoon_project/most_detectable_systems_final.txt', 'r')
lines = f.read().split('\n')
f.close()

for i in range(100): 
    sim_nos[i] = lines[i].split(' ')[0]
    dates.append(lines[i].split(' ')[1])
    res_nores_list.append(lines[i].split(' ')[2])

f = open('/home/garvit/Downloads/Exomoon_project/good_kepler_data.txt', 'r')
lightcurves_files = f.read().split('\n')[:-1]
f.close()

success_count = 0
used_keplerLC = []
g = open('/home/garvit/Downloads/Exomoon_project/leftover_systems2.txt', 'w')

for n in [1]:#range(100):
    old_success_count = success_count

    sim_no = int(sim_nos[n])
    res_nores = res_nores_list[n]
    date = dates[n]

    f = open("/media/Datas/Exomoon_project/final_model_LCs/{}_{}_{}.txt".format(sim_no, res_nores, date), 'r')
    model_data = f.readlines()
    f.close()

    f = open('/home/garvit/Downloads/Exomoon_project/sim_data/{}_{}/TTVsim{}_system_dictionary.pkl'.format(date, res_nores, sim_no), 'rb')
    sim_params= pickle.load(f)
    f.close()

    time, flux_planet_allmoons, flux_moons, params, num_moons, TT_planet_allmoons, transit_duration, moon_transit_depths = get_modelLCs(sim_params, 1471)
    model_time_interval = time[-1] - time[0]
    
    for i in range(1,len(lightcurves_files)):
        if (used_keplerLC.count(lightcurves_files[i]) != 0): continue
        
        f = open("/home/garvit/Downloads/Exomoon_project/good_kepler_data/{}".format(lightcurves_files[i]), 'r')
        kepler_data = f.readlines()
        f.close()
        kepler_starttime = float(kepler_data[1].split('\t')[0])
        kepler_endtime = float(kepler_data[-1].split('\t')[0])
        data_time_interval = kepler_endtime - kepler_starttime

        N = len(kepler_data) -1
        # if(data_time_interval < 950.0 or N<46000): continue

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

        data_gaps = kepler_data_gaps(timestamps)
        data_gaps = np.array(data_gaps)

        planet_transits10, best_offset = best_planet_transits()

        flux_planet_allmoons_temp = flux_planet_allmoons.copy()
        flux_moons_temp = flux_moons.copy()
        time_temp = time.copy()
        time_kepler = time - time[0] + kepler_starttime + best_offset
        if(data_time_interval < model_time_interval + best_offset):
            flux_planet_allmoons_temp = flux_planet_allmoons[time_kepler < kepler_endtime]
            flux_moons_temp = flux_moons[:, time_kepler  < kepler_endtime]
            time_temp = time[time_kepler < kepler_endtime]
            time_kepler = time_kepler[time_kepler < kepler_endtime]

        flux_total = np.copy(flux_planet_allmoons_temp) 
        for j in range(num_moons): flux_total+= flux_moons_temp[j]
        flux_total-=(num_moons)
  
        TT_kepler = TT_planet_allmoons - time[0] + kepler_starttime + best_offset
        TT_kepler = TT_kepler[TT_kepler < kepler_endtime]

        snr, moon_transits10, abs_best_snr = moon_transits_details()
        if((snr < 2.0).any()):
            print("SNR failed..................................................", end='\r')
            if((abs_best_snr < 2.0).any()): 
                break 
            continue

        actual_planet_transits = planet_transits10.sum()
        if(actual_planet_transits < 3): 
            print("number of actual transits is less than 3........", end='\r')
            continue
        
        moon_transit_fraction = np.ndarray((num_moons))
        for j in range(num_moons):
            moon_transit_fraction[j] = moon_transits10[j].sum()/actual_planet_transits

        if((moon_transit_fraction < 0.5).any()): 
            print("moon transit fraction failed.................", end='\r')
            continue
    
        trends = [0 for j in range(num_quarters)]
        detrended_fluxes = [0 for j in range(num_quarters)]
        for j in range(num_quarters):
            trends[j], detrended_fluxes[j], errors_detrend = cofiam_detrend(times=timestamps[quarters[j]], fluxes=kepler_flux[quarters[j]], errors=flux_error[quarters[j]], mask_idxs=[])

        trend = np.concatenate(trends)
        detrended_flux = np.concatenate(detrended_fluxes)
        timestamps = timestamps[np.concatenate(quarters)]

        residual = kepler_flux[np.concatenate(quarters)] - trend
        mask_idxs = mask_indices(timestamps)
        masked_residual = np.delete(residual, mask_idxs)
        dw_statistic = stattools.durbin_watson(masked_residual[np.logical_not(np.isnan(masked_residual))])
        if(abs( dw_statistic- 2.0) > 0.5 ): 
            print("DW failed..................................................", end='\r')
            continue

        print(snr , num_moons, dw_statistic, moon_transit_fraction)
        success_count+=1
    
        string = '{}_{}_{}'.format(date, res_nores, sim_no) + '\t' + lightcurves_files[i] + '\t' + str(num_moons) + '\t' + str(params.per_bary) + '\t' + ''.join(str(int(x)) for x in planet_transits10)
        for j in range(num_moons):
            string= string + '\t' + str(snr[j])
        for j in range(5 - num_moons):
            string= string + '\tNA'
        for j in range(num_moons):
            string= string + '\t' + ''.join(str(int(x)) for x in moon_transits10[j])
        for j in range(5 - num_moons):
            string= string + '\tNA'

        string = string + '\t' + str(dw_statistic) + '\t' + str(best_offset) + '\t' + str(transit_duration) + '\t' + str(TT_kepler[0])     
        h.write(string + '\n')

        used_keplerLC.append(lightcurves_files[i])

        break

    if(old_success_count == success_count):
        print('\nsystem failed. ', sim_no, date, res_nores, snr, num_moons, transit_duration) 
        g.write('{} {} {} {} {} {}\n'.format(date, res_nores, sim_no,snr, num_moons, transit_duration))
        
    

# moon_labels=['I', 'II', 'III', 'IV', 'V']
# colors=['red', 'green', 'brown','purple','cyan']

# plt.scatter(time, flux_planet_allmoons, color="blue", label='planet', s=0.5)

# for i in range(num_moons):
#     plt.scatter(time, flux_moons[i], color=colors[i], label= 'moon {}'.format(moon_labels[i]), s=0.5)

# #plt.plot(time, flux_total, color="black", label='total flux')

# #plt.ylim(.9996, 1.0001)
# #plt.xlim(min(time)+1*params.epoch_distance, min(time)+1*params.epoch_distance+ params.epoch_duration)
# plt.xlabel("Time (days)")
# plt.ylabel("Relative flux")
# #plt.legend()
# plt.show()
