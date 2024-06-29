
import numpy as np
import pandoramoon as pandora
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [12, 8]
import json
import os
import csv


def system_index(system):
    for i in range(len(catalogue)):
        if(catalogue[i].find(system+'\t')!=-1): return i


def mask_indices(time_stamps):

    mask_size = 3.0*transit_duration
    if(mask_size/per_bary > 0.5): mask_size = 0.5 * per_bary

    masks = []
    for i in range(num_epochs):
        masks.append(np.where((time_stamps > TT_kepler[i] - mask_size/2) & (time_stamps < TT_kepler[i] + mask_size/2))[0])

    mask_idxs = np.concatenate(masks)

    return mask_idxs



# systems = open('/home/garvit/Downloads/Exomoon_project/massratio_weirdos.txt', 'r').read().split('\n')[:-1]
systems = open('/home/garvit/Downloads/Exomoon_project/completed_systems.txt', 'r').read().split('\n')[:-1]

catalogue = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'r').read().split('\n')

fractions = [float(line.split(',')[2]) for line in open('/home/garvit/Downloads/Exomoon_project/excl_regions.txt', 'r').read().split('\n')[1:101]]

bflist = []
biclist = []
knlist1 = []
knlist2 = []
anlist1 = []
anlist2 = []

for i in range(100):
    # index = int(systems[i])-1
    # index = system_index(systems[i])
    system_line = catalogue[i+2].split('\t')

    system = system_line[0]
    date, resnores, simno = system.split('_') 
    # print(system)

    systems_details = json.load(open('/home/garvit/Downloads/Exomoon_project/sim_data/final_unpacked_details/{}.json'.format(system), 'r'))

    filename = 'final_lightcurve_{}.csv'.format(system)
    # f = open("/home/garvit/Downloads/Exomoon_project/final_lightcurves_collect/{}".format(filename), 'r')
    f = open("/home/garvit/Downloads/Exomoon_project/final_lightcurves_SNR/{}".format(filename), 'r')
    injected_lightcurve = f.readlines()
    f.close()

    N = len(injected_lightcurve) - 2
    flux_og = np.ndarray(N)
    timestamps = np.ndarray(N)

    for j in range(N):
        line = injected_lightcurve[j+2].split(',')
        timestamps[j] = float(line[0])
        flux_og[j] = float(line[1])

        
    kepler_file = system_line[1]
    f = open("/home/garvit/Downloads/Exomoon_project/good_kepler_data/{}".format(kepler_file), 'r')
    kepler_data = f.readlines()
    f.close()

    N = len(kepler_data) -1
    flux_error = np.ndarray((N))
    kepler_timestamps = np.ndarray((N))
    for j in range(N):
        line = kepler_data[j+1].split('\t')
        kepler_timestamps[j] = float(line[0])
        flux_error[j] = float(line[2])/float(line[1])

    noise_arr = np.ndarray((timestamps.shape[0]))
    for j in range(noise_arr.shape[0]):
        noise_arr[j] = flux_error[np.where(kepler_timestamps == timestamps[j])]

    new_N = int(float(fractions[system_index(system)-2])* timestamps.shape[0])
    timestamps = timestamps[:new_N]
    flux_og = flux_og[:new_N]
    noise_arr = noise_arr[:new_N]

    kepler_noise = np.median(noise_arr)

#_------------------------------------
    transit_duration = float(system_line[17])
    per_bary = float(system_line[3])
    num_epochs = len([c for c in system_line[4]])

    num_moons = len(systems_details.keys()) - 1

    f = open("/media/Datas/Exomoon_project/final_model_LCs/{}_{}_{}.txt".format(simno, resnores, date), 'r')
    model_data = f.readlines()
    f.close()

    TTs_model = model_data[0].split('\t')[:num_epochs]
    TTs = np.ndarray((num_epochs))
    for j in range(num_epochs):
        TTs[j] = TTs_model[j]
    first_TT = float(system_line[18])
    TT_kepler = TTs + (first_TT - TTs[0])

    mask_idx = mask_indices(timestamps)
    actual_noise = np.std(np.delete(flux_og, mask_idx))
    
    print(system, actual_noise)

#------------------------------------------------------------

    # model_type = 'planet_only'

    # model_fit_params1 = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}/{}/results.json'.format(system, model_type), 'r'))

    # params = pandora.model_params()
    # # params.per_bary, params.a_bary, params.r_planet, params.b_bary, params.ecc_bary, params.w_bary, params.t0_bary_offset, q1, q2 = model_fit_params1['posterior']['median']
    # params.per_bary, params.a_bary, params.r_planet, params.b_bary, params.ecc_bary, params.w_bary, params.t0_bary_offset, q1, q2 = model_fit_params1['maximum_likelihood']['point']


    # params.u1 = 2*q2*np.sqrt(q1)
    # params.u2 = -2*(q2 - 0.5) * np.sqrt(q1)

    # params.R_star = systems_details['Planet']['Rstar']  # [m]

    # # Planet parameters
    # params.t0_bary = float(system_line[18])  # [days]

    # params.M_planet = systems_details['Planet']['m']  # [kg]
    # # Moon parameters
    # params.r_moon = 1e-8  # negligible moon size
    # params.per_moon = 30  # other moon params do not matter
    # params.tau_moon = 0
    # params.Omega_moon = 0
    # params.i_moon = 0
    # params.ecc_moon = 0
    # params.w_moon = 0
    # params.M_moon = 1e-8  # negligible moon mass
        

    # # Other model parameters
    # params.epochs = len([c for c in system_line[4]])  # [int]
    # params.epoch_duration = float(system_line[17])*10  # [days]
    # params.cadences_per_day = 48  # [int]
    # params.epoch_distance = params.per_bary   # [days]
    # params.supersampling_factor = 1  # [int]
    # params.occult_small_threshold = 0.01  # [0..1]
    # params.hill_sphere_threshold = 1.2

    # model_planetonly = pandora.moon_model(params)

    # flux_planet_only, flux_planet, flux_moon = model_planetonly.light_curve(timestamps)

    # # -----------------------------------------------------------------------------
    # model_type = 'planet_moon'

    # model_fit_params2 = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}/{}/results.json'.format(system, model_type), 'r'))

    # params = pandora.model_params()
    # # params.per_bary, params.a_bary, params.r_planet, params.b_bary, params.ecc_bary, params.w_bary, params.t0_bary_offset, params.M_planet, params.r_moon, params.per_moon, params.tau_moon, params.Omega_moon, params.i_moon, params.M_moon, q1, q2 = model_fit_params2['posterior']['median']
    # params.per_bary, params.a_bary, params.r_planet, params.b_bary, params.ecc_bary, params.w_bary, params.t0_bary_offset, params.M_planet, params.r_moon, params.per_moon, params.tau_moon, params.Omega_moon, params.i_moon, params.M_moon, q1, q2 = model_fit_params2['maximum_likelihood']['point']


    # params.u1 = 2*q2*np.sqrt(q1)
    # params.u2 = -2*(q2 - 0.5) * np.sqrt(q1)

    # params.R_star = systems_details['Planet']['Rstar']  # [m]

    # # Planet parameters
    # params.t0_bary = float(system_line[18])  # [days]

    # params.ecc_moon = 0

    # # Other model parameters
    # params.epochs = len([c for c in system_line[4]])  # [int]
    # params.epoch_duration = float(system_line[17])*10  # [days]
    # params.cadences_per_day = 48  # [int]
    # params.epoch_distance = params.per_bary   # [days]
    # params.supersampling_factor = 1  # [int]
    # params.occult_small_threshold = 0.01  # [0..1]
    # params.hill_sphere_threshold = 1.2

    # model_planetmoon = pandora.moon_model(params)

    # flux_planet_moon, flux_planet, flux_moon = model_planetmoon.light_curve(timestamps)



    # evidence1 = model_fit_params1['logz']
    # evidence2 = model_fit_params2['logz']

    # bayes_factor = evidence2 - evidence1

    # snr = sum([float(j) for j in catalogue[system_index(systems[i])].split('\t')[5:5+num_moons]])/num_moons


    # # if(bayes_factor>100 and snr < 12): 
    # if(actual_noise/kepler_noise> 2.0):# and system not in scuti_systems):
    #     # print(system_index(systems[i])+1,systems[i], bayes_factor, snr)

    #     # k1 = 9
    #     # k2 = 16

    #     # chisq1 = np.sum(np.power((flux_og - flux_planet_only)/np.full(flux_og.shape, actual_noise),2))
    #     # chisq2 = np.sum(np.power((flux_og - flux_planet_moon)/np.full(flux_og.shape, actual_noise),2))

    #     # bic1 = k1*np.log(new_N) + chisq1
    #     # bic2 = k2*np.log(new_N) + chisq2

    #     # bflist.append(bayes_factor)
    #     # biclist.append(bic1 - bic2)

    #     # knlist1.append(kepler_noise)
    #     # anlist1.append(actual_noise)

    #     # print(system, bayes_factor, int(bic1 - bic2), kepler_noise, actual_noise, actual_noise/kepler_noise)
    #     # print(system_index(system)+1, system, actual_noise/kepler_noise, kepler_noise, actual_noise)

    # else:
    #     knlist2.append(kepler_noise)
    #     anlist2.append(actual_noise)


# plt.scatter(np.arange(len(knlist1)), np.array(anlist1)/np.array(knlist1), s=50, color='red')
# plt.scatter(np.arange(len(knlist2)), np.array(anlist2)/np.array(knlist2), s=50, color='blue')
# plt.xticks([])
# plt.yscale('log')
# plt.show()



# plt.plot(bflist, 2*np.array(bflist), label='y = 2x line', color='black', zorder = 1)
# plt.scatter(bflist, biclist, s=50, zorder=2, color='orange', edgecolors='black')
# plt.xlabel('log(K)')
# plt.ylabel('$\Delta BIC$')
# plt.xscale('log')
# plt.yscale('log')
# plt.legend()
# plt.show()

    
