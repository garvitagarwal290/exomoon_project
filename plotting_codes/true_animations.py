
import numpy as np
import pandoramoon as pandora
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [12, 8]
import json
import random
import math

sys_no = 44
index = sys_no - 1

catalogue = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'r').read().split('\n')
system_line = catalogue[index].split('\t')

system = system_line[0]
print(system)

sim_params = json.load(open('/home/garvit/Downloads/Exomoon_project/sim_data/final_unpacked_details/{}.json'.format(system), 'r'))

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

params.cadences_per_day = 500 #48 # [int]         #Using short cadence to generate model light curves and caluclating TTVs
params.epoch_distance = params.per_bary   # [days]
params.epochs = int(np.ceil(1400.0/(params.epoch_distance)))  # [int]
params.supersampling_factor = 1  # [int]
params.occult_small_threshold = 0.01  # [0..1]
params.hill_sphere_threshold = 1.1

N_epoch = params.epoch_duration*params.cadences_per_day
N = params.epochs* N_epoch

#Plug in moon parameters and obtain their LCs
omega_moon = random.random() * 2*np.pi
moon_key = 'I'
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


video = model.video(time,
    limb_darkening=True,
    teff=3000,
    planet_color="black",
    moon_color="black",
    ld_circles=200
)
video.save(filename="/home/garvit/Downloads/video3.mp4", fps=10, dpi=200)

