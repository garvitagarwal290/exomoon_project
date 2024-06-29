import matplotlib.pyplot as plt
# plt.rcParams['figure.figsize'] = [6, 4]
import numpy as np

import pandoramoon as pandora
from pandoramoon.helpers import ld_convert, ld_invert
import ultranest
import ultranest.stepsampler
from ultranest import ReactiveNestedSampler
from ultranest.plot import cornerplot

import pickle
import math
import random
import os
import csv
import sys
import json

from scipy.stats import norm,beta,truncnorm


def transform_uniform(x,lower,valrange):
	return lower + (valrange)*x

def transform_normal(x,mu,sigma):
	return norm.ppf(x,loc=mu,scale=sigma)

def transform_beta(x,lower,upper):
	return beta.ppf(x,lower,upper)

def transform_truncated_normal(x,mu,sigma,lower=0.,upper=1.):
	ar, br = (lower - mu) / sigma, (upper - mu) / sigma
	return truncnorm.ppf(x,ar,br,loc=mu,scale=sigma)


def prior_transform_planet_only(cube):
    p    = cube.copy()
    p[0] = per_bary* transform_normal(cube[0], 1.0 + perbary_shift, 0.01) # per_bary [days]

    p[1] = a_bary* transform_normal(cube[1], 1.0 + abary_shift, .1) # a_bary [R_star]

    p[2] = r_planet * transform_normal(cube[2], 1.0 + rp_shift, .1)  # r_planet [R_star]
    p[3] = cube[3] * 1             # b_bary [0..1]
    p[4] = cube[4] * 1			#ecc_bary
    p[5] = cube[5] * 360		#w_bary
    p[6] = transform_uniform(cube[6],-0.05, 0.1)     # t0_bary_offset [days]
    p[7] = cube[7]                  # LD q1 [0..1]
    p[8] = cube[8]                  # LD q2 [0..1]

    return p


def prior_transform_planet_moon(cube):
    p    = cube.copy()

    p[0] = per_bary* transform_normal(cube[0], 1.0 + perbary_shift, 0.01)  # per_bary [days]

    p[1] = a_bary* transform_normal(cube[1], 1.0 + abary_shift, .1) # a_bary [R_star]

    p[2] = r_planet * transform_normal(cube[2], 1.0 + rp_shift, .1) # r_planet [R_star]
    p[3] = cube[3] * 1             # b_bary [0..1]
    p[4] = cube[4] * 1                  #ecc_bary
    p[5] = cube[5] * 360                #w_bary
    p[6] = transform_uniform(cube[6],-0.05, 0.1)     # t0_bary_offset [days]
    
    p[7] = M_planet* transform_normal(cube[7], 1.0+ mplanet_shift, 0.05) #M_planet
    p[8] =  cube[8]*r_planet  #r_moon

    per_hill = (((r_hill)**3 * 4 * np.pi**2)/(6.674e-11 * M_planet))**.5
    per_hill /= (3600*24)
    per_lowerlim = ((( 2*systems_details['Planet']['Rp'])**3 * 4 * np.pi**2)/(6.674e-11 * M_planet))**.5
    per_lowerlim /= (3600*24)
    p[9] =  transform_uniform(cube[9], per_lowerlim, per_hill - per_lowerlim)   # per_moon

    p[10] =  cube[10] * 1 #tau
    p[11] =  cube[11] * 360 #Omega_moon
    p[12] =  cube[12] * 360 #i_moon
    p[13] =  M_planet * transform_uniform(cube[13], 1e-6, 1.0 - 1e-6)  #M_moon

    p[14] = cube[14]                  # LD q1 [0..1]
    p[15] = cube[15]                  # LD q2 [0..1]

    return p


def log_likelihood_planet_only(p):
    # Convert q priors to u LDs (Kipping 2013)
    q1 = p[7]
    q2 = p[8]
    u1, u2 = ld_convert(q1, q2)

    # Calculate pandora model with trial parameters
    _, _, flux_trial_total, _, _, _, _ = pandora.pandora(
        R_star = systems_details['Planet']['Rstar'],
        u1 = u1,
        u2 = u2,

        # Planet parameters
        per_bary = p[0],
        a_bary = p[1],
        r_planet = p[2],
        b_bary = p[3],
	ecc_bary = p[4],
        w_bary = p[5],
        t0_bary = float(system_line[18]),
        t0_bary_offset = p[6],   
        M_planet = systems_details['Planet']['m'],

        # Moon parameters
        r_moon = 1e-8,  # negligible moon size
        per_moon = 30,  # other moon params do not matter
        tau_moon = 0,
        Omega_moon = 0,
        i_moon = 0,
        ecc_moon = 0,
        w_moon = 0,
        M_moon = 1e-8,  # negligible moon mass

        # Other model parameters
        epoch_distance = systems_details['Planet']['Pp'],
        supersampling_factor = 1,
        occult_small_threshold = 0.01,
        hill_sphere_threshold=1.1,
        numerical_grid=25,
        time=timestamps,
        #cache=cache  # Can't use cache because free LDs
    )
    loglike = -0.5 * np.nansum(((flux_trial_total - flux) / noise)**2)

    return loglike


def log_likelihood_planet_moon(p):
    # Convert q priors to u LDs (Kipping 2013)
    q1 = p[14]
    q2 = p[15]
    u1, u2 = ld_convert(q1, q2)

    # Calculate pandora model with trial parameters
    _, _, flux_trial_total, _, _, _, _ = pandora.pandora(
	R_star =  R_star, #p[0],
        u1 = u1,
        u2 = u2,

        # Planet parameters
        per_bary = p[0],
        a_bary = p[1],
        r_planet = p[2],
        b_bary = p[3],
        ecc_bary = p[4],
	w_bary = p[5],
        t0_bary = float(system_line[18]),
        t0_bary_offset = p[6],   
        M_planet = p[7],

        # Moon parameters
        r_moon = p[8], 
        per_moon = p[9], 
        tau_moon = p[10],
        Omega_moon = p[11],
        i_moon = p[12],
        ecc_moon = 0,
        w_moon = 0,
        M_moon = p[13],  

        # Other model parameters
        epoch_distance = systems_details['Planet']['Pp'],
        supersampling_factor = 1,
        occult_small_threshold = 0.01,
        hill_sphere_threshold=1.1,
        numerical_grid=25,
        time=timestamps,
        #cache=cache  # Can't use cache because free LDs
    )
    loglike = -0.5 * np.nansum(((flux_trial_total - flux) / noise)**2)

    return loglike




index = int(sys.argv[1]) - 1
modeling_type = str(sys.argv[2])
if(modeling_type == ""): modeling_type = 'planet_only'

catalogue = open('/tiara/home/gagarwal/files_and_data/catalogue.txt', 'r').read().split('\n')
system_line = catalogue[index].split('\t')

system = system_line[0]
systems_details = json.load(open('/tiara/home/gagarwal/files_and_data/sim_data/final_unpacked_details/{}.json'.format(system), 'r'))

lline = open('/tiara/home/gagarwal/files_and_data/most_detectable_systems_final.txt', 'r').read().split('\n')[:-1][index-2]

filename = 'final_lightcurve_{}.csv'.format(system)
f = open("/tiara/home/gagarwal/files_and_data/final_lightcurves/{}".format(filename), 'r')
injected_lightcurve = f.readlines()
f.close()

kepler_file = injected_lightcurve[1].split(':')[1].strip()

N = len(injected_lightcurve) - 2
flux = np.ndarray(N)
timestamps = np.ndarray(N)

for j in range(N):
    line = injected_lightcurve[j+2].split(',')
    timestamps[j] = float(line[0])
    flux[j] = float(line[1])


f = open("/tiara/home/gagarwal/files_and_data/good_kepler_data/{}".format(kepler_file), 'r')
kepler_data = f.readlines()
f.close()

N = len(kepler_data) -1
flux_error = np.ndarray((N))
kepler_flux = np.ndarray((N))
kepler_timestamps = np.ndarray((N))
for j in range(N):
    line = kepler_data[j+1].split('\t')
    kepler_timestamps[j] = float(line[0])
    kepler_flux[j] = float(line[1])
    flux_error[j] = float(line[2])/float(line[1])
    
noise_arr = np.ndarray((timestamps.shape[0]))
for j in range(noise_arr.shape[0]):
    noise_arr[j] = flux_error[np.where(kepler_timestamps == timestamps[j])]

frac = float(sys.argv[3])
new_N = int(timestamps.shape[0] * frac)
timestamps = timestamps[:new_N]
flux = flux[:new_N]
noise_arr = noise_arr[:new_N]

noise = np.mean(noise_arr)

if(modeling_type == 'planet_moon'): 
	R_star = float(systems_details['Planet']['Rstar'])
	#print("R_star: {}\n".format(R_star))

per_bary = systems_details['Planet']['Pp']
a_bary = systems_details['Planet']['a']/systems_details['Planet']['Rstar']
r_planet = systems_details['Planet']['RpRstar']
w_bary = float(lline.split(' ')[3])
omega_moon = float(lline.split(' ')[4])

if(modeling_type == 'planet_only'):
	print('per_bary: {}\na_bary: {}\nr_planet: {}\nb_bary: {}\necc_bary: 0.0\nw_bary: {}\nq1: {}\nq2: {}'.format(per_bary, a_bary, r_planet, systems_details['Planet']['impact'],w_bary, systems_details['Planet']['q1'], systems_details['Planet']['q2']))

	parameters = ['per_bary','a_bary','r_planet','b_bary','ecc_bary', 'w_bary','t0_bary_offset','q1','q2']
	wrapped_params = [False,False,False,False,False, True,False,False,False]

elif(modeling_type == 'planet_moon'):
	print('per_bary: {}\na_bary: {}\nr_planet: {}\nb_bary: {}\necc_bary: 0.0\nw_bary: {}\n'.format(per_bary, a_bary, r_planet, systems_details['Planet']['impact'], w_bary))

	M_planet = systems_details['Planet']['m']
	r_moon = systems_details['I']['r']/R_star
	per_moon = systems_details['I']['P']/(24*3600)
	tau_moon = (2*np.pi-systems_details['I']['f'])/(2*np.pi)
	i_moon = (math.degrees(systems_details['I']['inc']) + 90.0)%360
	M_moon =  systems_details['I']['m']
	r_hill = systems_details['Planet']['RHill']
	a_moon = systems_details['I']['a']
	r_roche = systems_details['I']['r'] * (2*M_planet/M_moon)**(1/3)

	print('M_planet: {}\nr_moon: {}\nper_moon: {}\ntau_moon: {}\nomega_moon: {}\ni_moon: {}\nM_moon: {}\nq1: {}\nq2: {}\n'.format(M_planet, r_moon, per_moon, tau_moon,omega_moon, i_moon, M_moon, systems_details['Planet']['q1'], systems_details['Planet']['q2']))

	parameters = ['per_bary','a_bary','r_planet','b_bary','ecc_bary', 'w_bary','t0_bary_offset','M_planet','r_moon', 'per_moon','tau_moon', 'Omega_moon', 'i_moon', 'M_moon','q1','q2']
	wrapped_params = [False,False,False,False,False, True,False,False,False,False,True,True,True,False,False,False]


log_dir_planet_only="/tiara/home/gagarwal/modeling_results/{}_SR/planet_only/".format(system)
log_dir_planet_moon="/tiara/home/gagarwal/modeling_results/{}_SR/planet_moon/".format(system)

os.system('mkdir -p /tiara/home/gagarwal/files_and_data/prior_random_shifts/{}_SR/'.format(system))

if(len(sys.argv)==5 and sys.argv[4]=='resume'):
	resume = 'resume-similar'
	listt = open('/tiara/home/gagarwal/files_and_data/prior_random_shifts/{}_SR/prior_random_shifts_{}.txt'.format(system, modeling_type), 'r').readlines()
	perbary_shift= float(listt[0])
	abary_shift= float(listt[1])
	rp_shift = float(listt[2])
	mplanet_shift = float(listt[3])

else:
	resume = 'overwrite'
	perbary_shift = np.random.normal(scale=1e-5)
	abary_shift = np.random.normal(scale=1e-2)
	rp_shift = np.random.normal(scale=3* 1e-2)
	mplanet_shift = np.random.normal(scale=1e-2)

	f = open('/tiara/home/gagarwal/files_and_data/prior_random_shifts/{}_SR/prior_random_shifts_{}.txt'.format(system, modeling_type), 'w')
	f.write('{}\n{}\n{}\n{}'.format(perbary_shift, abary_shift,rp_shift, mplanet_shift))
	f.close()


if(modeling_type == 'planet_only'):
	sampler = ReactiveNestedSampler(
            parameters,
            log_likelihood_planet_only, 
            prior_transform_planet_only,
            wrapped_params=wrapped_params,
            resume= resume,
            log_dir = log_dir_planet_only
            )
elif(modeling_type == 'planet_moon'):
	sampler = ReactiveNestedSampler(
            parameters,
            log_likelihood_planet_moon, 
            prior_transform_planet_moon,
            wrapped_params=wrapped_params,
            resume = resume,
            log_dir=log_dir_planet_moon
            )


sampler.stepsampler = ultranest.stepsampler.RegionSliceSampler(
    nsteps=4000,
    adaptive_nsteps='move-distance',
    )  

#if(modeling_type == 'planet_only'):
#	result_planet_only = sampler.run(min_num_live_points=800)
#elif(modeling_type == 'planet_moon'):
#	result_planet_only = sampler.run(min_num_live_points=2000)

#result = sampler.run(min_num_live_points=1000)
result = sampler.run(region_class=ultranest.mlfriends.SimpleRegion, min_num_live_points=1000)


if(modeling_type == 'planet_only'):
        print('per_bary: {}\na_bary: {}\nr_planet: {}\nb_bary: {}\necc_bary: 0.0\nw_bary: {}\nq1: {}\nq2: {}'.format(per_bary, a_bary, r_planet, systems_details['Planet']['impact'],w_bary, systems_details['Planet']['q1'], systems_details['Planet']['q2']))


elif(modeling_type == 'planet_moon'):
        print('per_bary: {}\na_bary: {}\nr_planet: {}\nb_bary: {}\necc_bary: 0.0\nw_bary: {}\n'.format(per_bary, a_bary, r_planet, systems_details['Planet']['impact'], w_bary))
        print('M_planet: {}\nr_moon: {}\nper_moon: {}\ntau_moon: {}\nomega_moon: {}\ni_moon: {}\nM_moon: {}\nq1: {}\nq2: {}\n'.format(M_planet, r_moon, per_moon, tau_moon,omega_moon, i_moon, M_moon, systems_details['Planet']['q1'], systems_details['Planet']['q2']))

sampler.print_results()
