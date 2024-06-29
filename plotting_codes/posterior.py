
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [12, 8]
import json
from astropy.constants import M_earth, R_earth, au

planetormoon = 'moon'

params = ['per_bary', 'a_bary', 'r_planet', 'b_bary', 'ecc_bary', 'w_bary', 't0_bary_offset', 'M_planet','r_moon', 'per_moon', 'tau_moon', 'Omega_moon', 'i_moon', 'w_moon', 'M_moon', 'q1', 'q2']
param_toplot = 'per_moon'

params2 = ['Pp', 'a', 'Rp', 'impact', 'e', 'pomega', 'm', 'r', 'P', 'f', 'inc', 'q1', 'q2']
param_toplot2 = 'P'
pindex= params.index(param_toplot)

sys_no = 80
index = sys_no - 1

catalogue = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'r').read().split('\n')
system_line = catalogue[index].split('\t')

system = system_line[0] 
print(system)

systems_details = json.load(open('/home/garvit/Downloads/Exomoon_project/sim_data/final_unpacked_details/{}.json'.format(system), 'r'))

const = 1.0
if(param_toplot=='M_planet' or param_toplot == 'M_moon'): const = 1/M_earth.value
elif(param_toplot=='r_planet' or param_toplot == 'r_moon'): const = systems_details['Planet']['Rstar']/R_earth.value
elif(param_toplot=='a_bary'): const = systems_details['Planet']['Rstar']/au.value

const2 = 1.0
if(param_toplot2=='a'): const2 = 1/au.value
elif(param_toplot2=='Rp' or param_toplot2 == 'r'): const2 = 1/R_earth.value
elif(param_toplot2=='m'): const2 = 1/M_earth.value
elif(param_toplot2=='P'): const2 = 1/(3600*24)


f = open('/home/garvit/Downloads/Exomoon_project/modeling_results/{}/planet_moon/equal_weighted_post.txt'.format(system), 'r')
samples = np.loadtxt(f, skiprows=1)

posterior = samples[:,pindex]*const

weights = np.ones_like(posterior) / len(posterior)
hist = plt.hist(posterior, bins=20, weights=weights, color='red')[0]#, density=True)[0]
# print(hist)
plt.xlabel(param_toplot)
plt.title('posterior distribution of {}'.format(param_toplot))

if(planetormoon=='planet'):
    actual_param = systems_details['Planet'][param_toplot2]*const2
    plt.vlines(actual_param, ymin=0.0, ymax=np.max(hist), colors='blue', label='True Planet value', linewidth=3)

elif(planetormoon=='moon'):
    num_moons = len(systems_details.keys()) - 1
    keys = ['I', 'II', 'III', 'IV', 'V']
    colors = ['magenta', 'green', 'orange', 'pink', 'cyan']
    actual_params = np.ndarray((num_moons))
    for i in range(num_moons):
        actual_params[i] = systems_details[keys[i]][param_toplot2]*const2

        plt.vlines(actual_params[i], ymin=0.0, ymax=np.max(hist), colors=colors[i], label='Moon {}'.format(keys[i]), linewidth=3)

# plt.vlines(2.645, ymin=0.0, ymax=np.max(hist), colors='red', linewidth=3)

plt.legend()
plt.show()

