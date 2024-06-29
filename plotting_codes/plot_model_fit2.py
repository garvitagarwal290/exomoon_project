
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

# systems = open('/home/garvit/Downloads/Exomoon_project/massratio_weirdos.txt', 'r').read().split('\n')[:-1]
systems = open('/home/garvit/Downloads/Exomoon_project/completed_systems.txt', 'r').read().split('\n')[:-1]

catalogue = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'r').read().split('\n')

for i in range(19,len(systems)):
    # index = int(systems[i])-1
    index = system_index(systems[i])
    system_line = catalogue[index].split('\t')

    system = system_line[0]
    date, resnores, simno = system.split('_') 
    print(system)

    systems_details = json.load(open('/home/garvit/Downloads/Exomoon_project/sim_data/final_unpacked_details/{}.json'.format(system), 'r'))

    filename = 'final_lightcurve_{}.csv'.format(system)
    # os.system('cp /media/Datas/Exomoon_project/final_model_LCs/{}_{}_{}.txt /media/Datas/model_lightcurves/{}.txt'.format(simno, resnores, date, system))
    os.system('cp /home/garvit/Downloads/Exomoon_project/final_lightcurves_collect/{} /home/garvit/Downloads/injected_lightcurves/{}.txt'.format(filename, system))   
    f = open("/home/garvit/Downloads/Exomoon_project/final_lightcurves_collect/{}".format(filename), 'r')
    injected_lightcurve = f.readlines()
    f.close()

    N = len(injected_lightcurve) - 2
    flux_og = np.ndarray(N)
    timestamps = np.ndarray(N)

    for j in range(N):
        line = injected_lightcurve[j+2].split(',')
        timestamps[j] = float(line[0])
        flux_og[j] = float(line[1])

    full_time = np.arange(timestamps[0], timestamps[-1], .0208)


    model_type = 'planet_only'

    model_fit_params = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}/{}/results.json'.format(system, model_type), 'r'))

    params = pandora.model_params()
    # params.per_bary, params.a_bary, params.r_planet, params.b_bary, params.ecc_bary, params.w_bary, params.t0_bary_offset, q1, q2 = model_fit_params['posterior']['median']
    params.per_bary, params.a_bary, params.r_planet, params.b_bary, params.ecc_bary, params.w_bary, params.t0_bary_offset, q1, q2 = model_fit_params['maximum_likelihood']['point']


    params.u1 = 2*q2*np.sqrt(q1)
    params.u2 = -2*(q2 - 0.5) * np.sqrt(q1)

    params.R_star = systems_details['Planet']['Rstar']  # [m]

    # Planet parameters
    params.t0_bary = float(system_line[18])  # [days]

    params.M_planet = systems_details['Planet']['m']  # [kg]
    # Moon parameters
    params.r_moon = 1e-8  # negligible moon size
    params.per_moon = 30  # other moon params do not matter
    params.tau_moon = 0
    params.Omega_moon = 0
    params.i_moon = 0
    params.ecc_moon = 0
    params.w_moon = 0
    params.M_moon = 1e-8  # negligible moon mass
        

    # Other model parameters
    params.epochs = len([c for c in system_line[4]])  # [int]
    params.epoch_duration = float(system_line[17])*10  # [days]
    params.cadences_per_day = 48  # [int]
    params.epoch_distance = params.per_bary   # [days]
    params.supersampling_factor = 1  # [int]
    params.occult_small_threshold = 0.01  # [0..1]
    params.hill_sphere_threshold = 1.2

    model_planetonly = pandora.moon_model(params)

    flux_planet_only, flux_planet, flux_moon = model_planetonly.light_curve(timestamps)
    flux_planet_only_full, flux_planet, flux_moon = model_planetonly.light_curve(full_time)

    # -----------------------------------------------------------------------------
    model_type = 'planet_moon'

    model_fit_params = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}/{}/results.json'.format(system, model_type), 'r'))

    params = pandora.model_params()
    # params.per_bary, params.a_bary, params.r_planet, params.b_bary, params.ecc_bary, params.w_bary, params.t0_bary_offset, params.M_planet, params.r_moon, params.per_moon, params.tau_moon, params.Omega_moon, params.i_moon, params.M_moon, q1, q2 = model_fit_params['posterior']['median']
    params.per_bary, params.a_bary, params.r_planet, params.b_bary, params.ecc_bary, params.w_bary, params.t0_bary_offset, params.M_planet, params.r_moon, params.per_moon, params.tau_moon, params.Omega_moon, params.i_moon, params.M_moon, q1, q2 = model_fit_params['maximum_likelihood']['point']


    params.u1 = 2*q2*np.sqrt(q1)
    params.u2 = -2*(q2 - 0.5) * np.sqrt(q1)

    params.R_star = systems_details['Planet']['Rstar']  # [m]

    # Planet parameters
    params.t0_bary = float(system_line[18])  # [days]

    params.ecc_moon = 0

    # Other model parameters
    params.epochs = len([c for c in system_line[4]])  # [int]
    params.epoch_duration = float(system_line[17])*10  # [days]
    params.cadences_per_day = 48  # [int]
    params.epoch_distance = params.per_bary   # [days]
    params.supersampling_factor = 1  # [int]
    params.occult_small_threshold = 0.01  # [0..1]
    params.hill_sphere_threshold = 1.2

    model_planetmoon = pandora.moon_model(params)

    flux_planet_moon, flux_planet, flux_moon = model_planetmoon.light_curve(timestamps)
    flux_planet_moon_full, flux_planet, flux_moon = model_planetmoon.light_curve(full_time)

    # plt.scatter(timestamps, flux_og, s=2, color='black', label='simulated lightcurve')
    # # plt.plot(timestamps, flux_planet_only, color='red', linewidth=1, label='planet-only fit')
    # # plt.plot(timestamps, flux_planet_moon, color='blue',linewidth=1, label='planet + one moon fit')
    # plt.plot(full_time, flux_planet_only_full, color='red', linewidth=1, label='planet-only fit')
    # plt.plot(full_time, flux_planet_moon_full, color='blue',linewidth=1, label='planet + one moon fit')
    # plt.xlabel('Time (days)')
    # plt.ylabel('Flux')
    # plt.legend()
    # plt.show()

    f = open('/home/garvit/Downloads/model_fit_lightcurves/{}_planet_only.csv'.format(system), 'w')
    f.write('time,flux\n')
    writer1 = csv.writer(f, delimiter=',')
    writer1.writerows(zip(timestamps, flux_planet_only))
    f.close()

    f = open('/home/garvit/Downloads/model_fit_lightcurves/{}_planet_moon.csv'.format(system), 'w')
    f.write('time,flux\n')
    writer1 = csv.writer(f, delimiter=',')
    writer1.writerows(zip(timestamps, flux_planet_moon))
    f.close()

    os.system('cp /media/Datas/Exomoon_project/modeling_results/{}/planet_only/results.json /home/garvit/Downloads/modeling_results/{}_planet_only.json'.format(system, system))
    os.system('cp /media/Datas/Exomoon_project/modeling_results/{}/planet_moon/results.json /home/garvit/Downloads/modeling_results/{}_planet_moon.json'.format(system, system))

    px, py, mx, my = model_planetonly.coordinates(timestamps)

    extent = 2.5

    good_px_idxs = np.where((px >= -extent) & (px <= extent))[0]
    good_py_idxs = np.where((py >= -extent) & (py <= extent))[0]
    #### good planet indices are ones where both x and y coordinate requirements are satisfied
    good_planet_idxs = np.intersect1d(good_px_idxs, good_py_idxs)

    good_mx_idxs = np.where((mx >= -extent) & (mx <= extent))[0]
    good_my_idxs = np.where((my >= -extent) & (my <= extent))[0]
    #### good moon indices are ones where both x and y coordinate requirements are satisfied
    good_moon_idxs = np.intersect1d(good_mx_idxs, good_my_idxs)

    #### want to animate when either the planet or the moon is within the bounds
    animate_idxs = np.union1d(good_planet_idxs, good_moon_idxs)

    #### now update px, py, mx, my
    model_planetonly.px = px[animate_idxs]
    model_planetonly.py = py[animate_idxs]
    model_planetonly.mx = mx[animate_idxs]
    model_planetonly.my = my[animate_idxs]
    timestamps_anim = timestamps[animate_idxs]

    # Create video
    video = model_planetonly.video(
    time=timestamps_anim,
    limb_darkening=True,
    teff=3000,
    planet_color="blue",
    moon_color="white",
    ld_circles=100
    )
    # Save video to disk
    video.save(filename="/media/Datas/Exomoon_project/modeling_results/animations/{}_planet_only.mp4".format(system), fps=5, dpi=200)





    px, py, mx, my = model_planetmoon.coordinates(timestamps)

    extent = 2.5

    good_px_idxs = np.where((px >= -extent) & (px <= extent))[0]
    good_py_idxs = np.where((py >= -extent) & (py <= extent))[0]
    #### good planet indices are ones where both x and y coordinate requirements are satisfied
    good_planet_idxs = np.intersect1d(good_px_idxs, good_py_idxs)

    good_mx_idxs = np.where((mx >= -extent) & (mx <= extent))[0]
    good_my_idxs = np.where((my >= -extent) & (my <= extent))[0]
    #### good moon indices are ones where both x and y coordinate requirements are satisfied
    good_moon_idxs = np.intersect1d(good_mx_idxs, good_my_idxs)

    #### want to animate when either the planet or the moon is within the bounds
    animate_idxs = np.union1d(good_planet_idxs, good_moon_idxs)

    #### now update px, py, mx, my
    model_planetmoon.px = px[animate_idxs]
    model_planetmoon.py = py[animate_idxs]
    model_planetmoon.mx = mx[animate_idxs]
    model_planetmoon.my = my[animate_idxs]
    timestamps_anim = timestamps[animate_idxs]

    if(model_planetmoon.mx.shape[0]==0):continue

    video = model_planetmoon.video(timestamps_anim,
    limb_darkening=True,
    teff=3000,
    planet_color="blue",
    moon_color="white",
    ld_circles=100
    )
    video.save(filename="/media/Datas/Exomoon_project/modeling_results/animations/{}_planet_moon.mp4".format(system), fps=5, dpi=200)

