
import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [12, 8]
font = {'size'   : 15}
matplotlib.rc('font', **font)

def system_index(system):
    for i in range(len(catalogue)):
        if(catalogue[i].find(system+'\t')!=-1): return i

systems = open('/home/garvit/Downloads/Exomoon_project/completed_systems.txt', 'r').read().split('\n')[:-1]
bayes_factors = [[] for i in range(5)]

catalogue = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'r').read().split('\n')
new = open('/home/garvit/Downloads/catalogue.txt', 'w')

n=0
for i in range(len(systems)):

    index = system_index(systems[i])

    sim_data = json.load(open('/home/garvit/Downloads/Exomoon_project/sim_data/final_unpacked_details/{}.json'.format(systems[i]), 'r'))
    num_moons = len(sim_data.keys()) - 1

    results = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}/planet_only/results.json'.format(systems[i]), 'r'))
    evidence1 = results['logz']

    results = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}/planet_moon/results.json'.format(systems[i]), 'r'))
    evidence2 = results['logz']

    bayes_factor = evidence2 - evidence1
    bayes_factors[num_moons-1].append(bayes_factor)

    new.write(catalogue[index]+'\t'+str(bayes_factor)+'\n')


    # print(systems[i], bayes_factor)

print(n)
plt.hlines(0.0, 0, 6, color='black')

for i in range(5):
    plt.scatter([i+1 for j in range(len(bayes_factors[i]))], bayes_factors[i], s=50,marker='+', color='blue')

plt.xlabel('Number of Moons')
plt.ylabel('log(Bayes factor)')
plt.xticks([1,2,3, 4, 5])
# plt.ylim(top=threshold, bottom=-100)
plt.yscale('symlog')
plt.show()