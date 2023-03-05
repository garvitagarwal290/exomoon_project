import matplotlib.pyplot as plt
#plt.rcParams['figure.figsize'] = [6, 4]
import numpy as np

import pandoramoon as pandora

import pickle
import math
import random
import os
import csv

from wotan import flatten, t14
from statsmodels.stats import stattools


def gen_separate_modelLCs(sim_params, data_time_interval):

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

    # w_bary = sim_params['Planet']['pomega']
    # params.w_bary = math.degrees(w_bary) if w_bary!=None else random.random()*2*np.pi  # [deg]
    params.w_bary = w_bary
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
    params.epochs = int(np.ceil(data_time_interval/(params.epoch_distance)))  # [int]
    params.supersampling_factor = 1  # [int]
    params.occult_small_threshold = 0.01  # [0..1]
    params.hill_sphere_threshold = 1.1

    num_moons = len(sim_params.keys()) - 1
    N_epoch = params.epoch_duration*params.cadences_per_day
    N = params.epochs* N_epoch
    TTs= np.ndarray((num_moons+1, params.epochs)) #array for Transit Times
    
    #Plug in moon(s) parameters and obtain their LCs
    flux_moons= np.ndarray((num_moons, N))
    # omega_moon = random.random() * 2*np.pi
    for i in range(num_moons):
        moon_key = list(sim_params.keys())[i+1]
        params.r_moon = sim_params[moon_key]['r']/params.R_star  # [R_star]
        params.per_moon = sim_params[moon_key]['P']/(24*3600) # [days]
        params.a_moon = sim_params[moon_key]['a']/params.R_star  # [R_star]
        params.tau_moon = (2*np.pi-sim_params[moon_key]['f'])/(2*np.pi)  # [0..1]
        # params.Omega_moon = math.degrees(omega_moon)  # [0..180]
        params.Omega_moon = omega_moon
        params.i_moon = (math.degrees(sim_params[moon_key]['inc']) + 90.0)%360  # [0..360]
        params.e_moon = sim_params[moon_key]['e']  # [0..1]
        params.w_moon = math.degrees(sim_params[moon_key]['pomega'])  # [deg]
        params.M_moon = sim_params[moon_key]['m']  #kg

        time = pandora.time(params).grid()
        model = pandora.moon_model(params)
        flux_total, flux_planet, flux_moons[i] = model.light_curve(time)

        for j in range(params.epochs):
            TTs[i+1][j] = time[j* N_epoch + np.nanargmin(flux_planet[j*N_epoch:(j+1)*N_epoch])]

    #calculate planet's LC without any effect of moons
    params.M_moon = 1e-8
    model = pandora.moon_model(params)
    flux_total, flux_planet_only, flux_moon = model.light_curve(time)

    for j in range(params.epochs):
        TTs[0][j] = time[j* N_epoch + np.nanargmin(flux_planet_only[j*N_epoch:(j+1)*N_epoch])]
    
    return time, flux_planet_only, flux_moons, params, num_moons, TTs, transit_duration


def shift(arr, delta=0, fill_value=1.0):
        result = np.empty_like(arr)
        if delta > 0:
            result[:delta] = fill_value
            result[delta:] = arr[:-delta]
        elif delta < 0:
            result[delta:] = fill_value
            result[:delta] = arr[-delta:]
        else:
            result[:] = arr
        
        return result

def calculateTTVs(params, num_moons, TTs, flux_planet_only, flux_moons):

    TTVs= np.ndarray((num_moons, params.epochs)) 

    x = [i for i in range(params.epochs)]
    for i in range(num_moons):
        line_fit = np.polyfit(x, y=TTs[i+1], deg=1)
        TTVs[i] = TTs[i+1] - np.polyval(line_fit, x)


    TTV_planet_allmoons = np.sum(TTVs, axis=0)
    TT_planet_allmoons = np.add(TTV_planet_allmoons, TTs[0])

    N_epoch = params.epoch_duration*params.cadences_per_day

    flux_planet_allmoons = np.empty_like(flux_planet_only)
    for i in range(params.epochs):
        transit = flux_planet_only[i*N_epoch:(i+1)*N_epoch]
        delta = int(np.round(TTV_planet_allmoons[i] * params.cadences_per_day))
        flux_planet_allmoons[i*N_epoch:(i+1)*N_epoch] = shift(transit, delta)

    
    flux_total = np.copy(flux_planet_allmoons) 
    for i in range(num_moons): flux_total+= flux_moons[i]
    flux_total-=(num_moons)

    return TT_planet_allmoons, flux_planet_allmoons, flux_total


sim_nos = np.ndarray((100))
dates = []
res_nores_list = []
wbary_arr = np.ndarray((100))
omegemoon_arr = np.ndarray((100))
f = open('/home/garvit/Downloads/Exomoon_project/most_detectable_systems_july13.txt', 'r')
lines = f.read().split('\n')
f.close()

for i in range(100): 
    sim_nos[i] = lines[i].split(' ')[0]
    dates.append(lines[i].split(' ')[1])
    res_nores_list.append(lines[i].split(' ')[2])
    wbary_arr[i] = lines[i].split(' ')[3]
    omegemoon_arr[i] = lines[i].split(' ')[4]


for n in range(100):

    sim_no = int(sim_nos[n])
    print(sim_no)
    res_nores = res_nores_list[n]
    w_bary = wbary_arr[n]
    omega_moon = omegemoon_arr[n]
    date = dates[n]

    f = open('/home/garvit/Downloads/Exomoon_project/sim_data/{}_{}/TTVsim{}_system_dictionary.pkl'.format(date, res_nores, sim_no), 'rb')
    sim_params= pickle.load(f)
    f.close()

    data_time_interval = 1471.0 

    time, flux_planet_only, flux_moons, params, num_moons, TTs, transit_duration = gen_separate_modelLCs(sim_params, data_time_interval)
    if(params.epochs<3): 
        print("less than 3 epochs! ",sim_no)

    TT_planet_allmoons, flux_planet_allmoons, flux_total = calculateTTVs(params, num_moons, TTs, flux_planet_only, flux_moons)

    f = open('/media/Datas/Exomoon_project/model_LCs_{}/{}_{}_{}.txt'.format(date, sim_no, res_nores, date), 'w')

    for j in range(params.epochs):  
        f.write('{}\t'.format(TT_planet_allmoons[j]))
    f.write('\n')

    for j in range(time.shape[0]):
        f.write('{}\t{}'.format(time[j], flux_planet_allmoons[j]))

        for k in range(num_moons):
            f.write('\t{}'.format(flux_moons[k][j]))

        f.write('\n')

    f.close()
