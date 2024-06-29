import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [12, 8]
font = {'size'   : 15}

matplotlib.rc('font', **font)
import json


systems = open('/home/garvit/Downloads/Exomoon_project/completed_systems.txt', 'r').read().split('\n')[:-1]

derived_ratio = [[] for i in range(5)]
actual_ratio = [[] for i in range(5)]

for i in range(len(systems)):

    system_details = json.load(open('/home/garvit/Downloads/Exomoon_project/sim_data/final_unpacked_details/{}.json'.format(systems[i]), 'r'))
    num_moons = len(system_details.keys()) - 1

    results = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}/planet_moon/results.json'.format(systems[i]), 'r'))
    derived_ratio[num_moons-1].append(results['maximum_likelihood']['point'][13]/results['maximum_likelihood']['point'][7])
    # derived_ratio[num_moons-1].append(results['posterior']['median'][13]/results['posterior']['median'][7])

    total_moonmass = 0
    keys = ['I', 'II', 'III', 'IV', 'V']
    for j in range(int(num_moons)):
        total_moonmass += system_details[keys[j]]['m']

    actual_ratio[num_moons-1].append(total_moonmass/system_details['Planet']['m'])


colors = ['red', 'green', 'blue', 'magenta', 'Gold']

plt.plot([0,0.015], [0,0.015], color='black', label='x = y line')
for i in range(5):
    plt.scatter(actual_ratio[i], derived_ratio[i], s=50, c=colors[i], label='No of Moons: {}'.format(i+1))

plt.xlabel('actual $(total\,\,M_M)/M_P$ ratio')
plt.ylabel('derived $M_M/M_P$ ratio')
plt.legend()
plt.show()