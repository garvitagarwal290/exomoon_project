import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# plt.rcParams['figure.figsize'] = [8,12]
font = {'size'   : 15}
matplotlib.rc('font', **font)

def system_index(system):
    for i in range(len(catalogue)):
        if(catalogue[i].find(system+'\t')!=-1): return i

systems = open('/home/garvit/Downloads/Exomoon_project/completed_systems.txt', 'r').read().split('\n')[:-1]
catalogue = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'r').read().split('\n')[:-1]
SNRs = [[] for i in range(5)]
bayes_factors = [[] for i in range(5)]

n=0
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

    # total_moonmass = 0
    # keys = ['I', 'II', 'III', 'IV', 'V']
    # for j in range(int(num_moons)):
    #     total_moonmass += sim_data[keys[j]]['m']

    # if(bayes_factor<0): 
    #     print(system_index(systems[i])+1,systems[i], bayes_factor, snr)
        # print(snr, bayes_factor)
        # print(results['posterior']['median'][7]/results['posterior']['median'][13])

    # derived_ratio = results['maximum_likelihood']['point'][7]/results['maximum_likelihood']['point'][13]   
    # if(derived_ratio>10 and bayes_factor<0):# and bayes_factor>0):
    #     n+=1
        # print(system_index(systems[i])+1, systems[i], snr, bayes_factor, derived_ratio)
        # print(system_index(systems[i])+1, systems[i], bayes_factor)
    print(systems[i], bayes_factor)


print(n)
colors = ['red', 'green', 'blue', 'magenta', 'Gold']
for i in range(5):
    plt.scatter(SNRs[i], bayes_factors[i], s=50,marker='+', color=colors[i], label='No of Moons: {}'.format(i+1))

plt.xlabel('Average SNR of Moons in the system')
plt.ylabel('log(Bayes factor)')
# plt.ylim(-100, threshold1)
# plt.xlim(0,threshold2)
plt.legend()
plt.show()