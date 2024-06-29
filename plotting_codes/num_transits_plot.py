import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [8,12]
font = {'size'   : 15, 'family':'serif'}
matplotlib.rc('font', **font)

def system_index(system):
    for i in range(len(catalogue)):
        if(catalogue[i].find(system+'\t')!=-1): return i

systems = open('/home/garvit/Downloads/Exomoon_project/completed_systems.txt', 'r').read().split('\n')[:-1]
catalogue = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'r').read().split('\n')[:-1]
SNRs = [[] for i in range(5)]
bayes_factors = [[] for i in range(5)]
transits_arr = [[] for i in range(5)]
derived_ratios = [[] for i in range(5)]

fracs = open('/home/garvit/Downloads/Exomoon_project/fraction&bad_detrending.txt','r').readlines()
exceptions = {'july13_nores_421':16, 'july13_nores_88':9}

for i in range(len(systems)):

    sim_data = json.load(open('/home/garvit/Downloads/Exomoon_project/sim_data/final_unpacked_details/{}.json'.format(systems[i]), 'r'))
    num_moons = len(sim_data.keys()) - 1

    results = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}/planet_only/results.json'.format(systems[i]), 'r'))
    evidence1 = results['logz']

    results = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}/planet_moon/results.json'.format(systems[i]), 'r'))
    evidence2 = results['logz']

    bayes_factor = evidence2 - evidence1
    bayes_factors[num_moons-1].append(bayes_factor)

    snr = sum([float(j) for j in catalogue[system_index(systems[i])].split('\t')[5:5+num_moons]])/num_moons
    SNRs[num_moons-1].append(snr)

    derived_ratio = results['maximum_likelihood']['point'][13]/results['maximum_likelihood']['point'][7]
    derived_ratios[num_moons-1].append(derived_ratio)

    epochs = len([j for j in catalogue[system_index(systems[i])].split('\t')[4]])

    if(systems[i] not in exceptions.keys()): num_transits = float(fracs[system_index(systems[i])].split(' ')[4]) * epochs
    else: num_transits = exceptions[systems[i]]
    transits_arr[num_moons-1].append(num_transits)

    # if(bayes_factor<0.0): 
    #     print(systems[i], snr, bayes_factor, num_transits)



colors = ['red', 'green', 'blue', 'magenta', 'Gold']
markers = ['P', 'o', 'v', 'X', 's']
for i in range(5):
    plt.scatter(transits_arr[i], bayes_factors[i], s=50,marker=markers[i], color=colors[i], label='No of Moons: {}'.format(i+1))

plt.xlabel('Number of transit epochs modeled')
plt.ylabel('log(Bayes factor)')
plt.legend()
plt.yscale('symlog')
plt.show()






