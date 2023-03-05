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

    w_bary = sim_params['Planet']['pomega']
    if(w_bary ==None): w_bary = random.random()*2*np.pi
    params.w_bary = math.degrees(w_bary)  # [deg]
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
    omega_moon = random.random() * 2*np.pi
    for i in range(num_moons):
        moon_key = list(sim_params.keys())[i+1]
        params.r_moon = sim_params[moon_key]['r']/params.R_star  # [R_star]
        params.per_moon = sim_params[moon_key]['P']/(24*3600) # [days]
        params.a_moon = sim_params[moon_key]['a']/params.R_star  # [R_star]
        params.tau_moon = (2*np.pi-sim_params[moon_key]['f'])/(2*np.pi)  # [0..1]
        params.Omega_moon = math.degrees(omega_moon)  # [0..180]
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


def moon_transits_details(params, flux_moons):

    N_transit = int(np.round(transit_duration * 48))

    transit_depths= np.ndarray(num_moons)
    moon_transits10 = ['' for i in range(num_moons)]
    N = int(flux_moons.shape[1]/params.epochs)
    for i in range(num_moons):
        transit_depths[i] = 1.0 - np.min(flux_moons[i])

        for j in range(params.epochs):
            if((flux_moons[i][j*N : (j+1)*N]==1.0).all()): moon_transits10[i]+='0'
            else: moon_transits10[i]+='1'

    return moon_transits10, transit_depths


for n in range(1000):

    res_nores = 'res'
    date= 'july13'
    
    f = open('/home/garvit/Downloads/Exomoon_project/sim_data/{}_{}/TTVsim{}_system_dictionary.pkl'.format(date, res_nores, n+1), 'rb')
    sim_params= pickle.load(f)
    f.close()

    
    data_time_interval = 1471.0

    time, flux_planet_only, flux_moons, params, num_moons, TTs, transit_duration = gen_separate_modelLCs(sim_params, data_time_interval)
    if(params.epochs<3): 
        continue
    TT_planet_allmoons, flux_planet_allmoons, flux_total = calculateTTVs(params, num_moons, TTs, flux_planet_only, flux_moons)
    
    moon_transits10, transit_depths = moon_transits_details(params, flux_moons)
    
    transit_fractions = []
    for j in range(num_moons):
        m = np.int32(list(moon_transits10[j])).sum()
        transit_fractions.append(m/params.epochs) 
        
    if( (np.array(transit_fractions) >= 0.5).all()):
    	print(n+1, date, params.w_bary, params.Omega_moon)

    #if( (np.array(transit_fractions) >= 0.5).all()):# and num_moons==5):
     #   good+=1 
      #  #print(n+1, good, num_moons, transit_fractions)
       # print(n+1, date)
