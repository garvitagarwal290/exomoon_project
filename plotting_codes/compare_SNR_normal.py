import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [8,12]
font = {'size'   : 15, 'family':'serif'}
matplotlib.rc('font', **font)

def system_index(system):
    for i in range(len(catalogue)):
        if(catalogue[i].find(system)!=-1): return i


catalogue = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'r').read().split('\n')[:-1]
systems = open('/home/garvit/Downloads/Exomoon_project/completed_systems2.txt', 'r').read().split('\n')[:-1]
bayes_factors = [[] for i in range(5)]
bayes_factors_SNR = [[] for i in range(5)]
SNRs = [[] for i in range(5)]


n=0
for i in range(len(systems)):

    sim_data = json.load(open('/home/garvit/Downloads/Exomoon_project/sim_data/final_unpacked_details/{}.json'.format(systems[i]), 'r'))
    num_moons = len(sim_data.keys()) - 1

    results = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}/planet_only/results.json'.format(systems[i]), 'r'))
    evidence1 = results['logz']

    results = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}/planet_moon/results.json'.format(systems[i]), 'r'))
    evidence2 = results['logz']

    bayes_factor1 = evidence2 - evidence1
    bayes_factors[num_moons-1].append(bayes_factor1)


    results = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}_SNR/planet_only/results.json'.format(systems[i]), 'r'))
    evidence1 = results['logz']

    results = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}_SNR/planet_moon/results.json'.format(systems[i]), 'r'))
    evidence2 = results['logz']

    bayes_factor2 = evidence2 - evidence1
    bayes_factors_SNR[num_moons-1].append(bayes_factor2)

    snr = sum([float(j) for j in catalogue[system_index(systems[i])].split('\t')[5:5+num_moons]])/num_moons
    SNRs[num_moons-1].append(snr)

    print(system_index(systems[i])+1,systems[i], bayes_factor1, bayes_factor2, snr )

    # print(systems[i], bayes_factor)

print(n)
# plt.hlines(0.0, 0, 6, color='black')

# for i in range(5):
    # plt.scatter([i+1 for j in range(len(bayes_factors[i]))], bayes_factors[i], s=50,marker='X', color='blue', label='normal')
    # plt.scatter([i+1 for j in range(len(bayes_factors_SNR[i]))], bayes_factors_SNR[i], s=50,marker='o', color='red', label='SNR')

bayes_factors_flat = [i for sublist in bayes_factors for i in sublist]
bayes_factors_SNR_flat = [i for sublist in bayes_factors_SNR for i in sublist]
snr_flat = [i for sublist in SNRs for i in sublist]

plt.plot(snr_flat[0], bayes_factors_flat[0],markersize=10,marker='X', color='blue', label='original runs')
plt.plot(snr_flat[0], bayes_factors_SNR_flat[0],markersize=10,marker='o', color='red', label='common SNR runs')
for i in range(1,len(systems)):
    plt.plot(snr_flat[i], bayes_factors_flat[i],markersize=10,marker='X', color='blue')
    plt.plot(snr_flat[i], bayes_factors_SNR_flat[i],markersize=10,marker='o', color='red')

# plt.xlabel('Number of Moons')
plt.xlabel('Avg moon SNR in original run')
plt.ylabel('log(Bayes factor)')
# plt.xticks([1,2,3, 4, 5])
# plt.xticks([])
# plt.ylim(top=threshold, bottom=-100)
plt.yscale('symlog')
plt.legend()
plt.show()