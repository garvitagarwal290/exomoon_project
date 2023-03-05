
import numpy as np
import pickle
import json


sim_nos = np.ndarray((100))
dates = []
res_nores_list = []
f = open('/home/garvit/Downloads/Exomoon_project/most_detectable_systems_july13.txt', 'r')

lines = f.read().split('\n')
f.close()

for i in range(100): 
    sim_nos[i] = lines[i].split(' ')[0]
    dates.append(lines[i].split(' ')[1])
    res_nores_list.append(lines[i].split(' ')[2])


for n in range(100):
    sim_no = int(sim_nos[n])
    res_nores = res_nores_list[n]
    date = dates[n]

    f = open('/home/garvit/Downloads/Exomoon_project/sim_data/{}_{}/TTVsim{}_system_dictionary.pkl'.format(dates, res_nores, sim_no), 'rb')
    sim_params= pickle.load(f)
    f.close()

    f = open('/home/garvit/Downloads/Exomoon_project/sim_data/unpacked_details_{}/{}_{}_{}.json'.format(dates, dates, res_nores, sim_no), 'w')
    json.dump(sim_params, f)
    # for key, value in sim_params.items():
    #     f.write(str(key)+': '+str(value)+'\n')
    # f.close()