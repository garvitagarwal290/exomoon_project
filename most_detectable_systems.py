import pickle
import numpy as np

top_n_moons = 10

f = open('/home/garvit/Downloads/Exomoon_project/most_detectable_systems_july13.txt', 'w')

greatest_rmrs_ratio_nores = [[] for i in range(5)]
greatest_rmrs_ratio_res = [[] for i in range(5)]
system_indices_res=[[0 for i in range(top_n_moons)] for j in range(5)]
system_indices_nores=[[0 for i in range(top_n_moons)] for j in range(5)]


def cal_ratios(a):
    rmrs_median_ratios= [[] for i in range(5)]
    all_median_ratios = []
    
    lines = open('/home/garvit/Downloads/Exomoon_project/one_time_use/best_transit_fraction_systems_{}.txt'.format(a), 'r').read().split('\n')[:-1]

    for i in range(len(lines)):
        date = lines[i].split(' ')[1]
        sim_no = lines[i].split(' ')[0]
        f = open('/home/garvit/Downloads/Exomoon_project/sim_data/{}_{}/TTVsim{}_system_dictionary.pkl'.format(date, a, sim_no), 'rb')

        sim_params= pickle.load(f)

        
        num_moons= len(sim_params.keys()) - 1

        rmrs_ratios=np.ndarray(num_moons)
        for j in range(num_moons):
            moon_key = list(sim_params.keys())[j+1]

            rmrs_ratios[j] = sim_params[moon_key]['r']/sim_params['Planet']['Rstar']

        median_ratio = np.median(rmrs_ratios)
        rmrs_median_ratios[num_moons-1].append(median_ratio)
        all_median_ratios.append(median_ratio)

        '''
        biggest_rmrs_ratio=0.0
        for j in range(num_moons):
            moon_key = list(sim_params.keys())[j+1]

            rmrs_ratios = sim_params[moon_key]['r']/sim_params['Planet']['Rstar']
            if(rmrs_ratio> biggest_rmrs_ratio): biggest_rmrs_ratio = rmrs_ratio
        '''
            
        #rmrs_ratios[num_moons-1].append(biggest_rmrs_ratio)
        #all_ratios[i] = biggest_rmrs_ratio

    return rmrs_median_ratios, all_median_ratios, lines



# f.write('no resonance systems\n')
# rmrs_ratios_nores, ratios, lines_nores = cal_ratios('nores')

# for i in range(5): 
#     rmrs_ratios_nores[4-i].sort(reverse=True)
#     greatest_rmrs_ratio_nores[4-i] = rmrs_ratios_nores[4-i][:top_n_moons]
#     greatest_rmrs_ratio_nores[4-i].sort()

#     f.write('\nnum moons = {}\n'.format(5-i))
#     for j in range(top_n_moons):
#         system_indices_nores[4-i][j] = ratios.index(greatest_rmrs_ratio_nores[4-i][j])
#         f.write(lines_nores[system_indices_nores[4-i][j]].split(' ')[0]+' '+lines_nores[system_indices_nores[4-i][j]].split(' ')[1]+" nores\n")


# f.write('\n\nresonance systems\n')
# rmrs_ratios_res, ratios, lines_res = cal_ratios('res')

# for i in range(5): 
#     rmrs_ratios_res[4-i].sort(reverse=True)
#     greatest_rmrs_ratio_res[4-i] = rmrs_ratios_res[4-i][:top_n_moons]
#     greatest_rmrs_ratio_res[4-i].sort()

#     f.write('\nnum moons = {}\n'.format(5-i))
#     for j in range(top_n_moons):
#         system_indices_res[4-i][j] = ratios.index(greatest_rmrs_ratio_res[4-i][j])
#         f.write(lines_res[system_indices_res[4-i][j]].split(' ')[0]+' '+lines_res[system_indices_res[4-i][j]].split(' ')[1]+" res\n")


rmrs_ratios_res, ratios, lines_res = cal_ratios('res')

for i in range(5): 
    rmrs_ratios_res[4-i].sort(reverse=True)
    greatest_rmrs_ratio_res[4-i] = rmrs_ratios_res[4-i][:top_n_moons]
    greatest_rmrs_ratio_res[4-i].sort()

    for j in range(top_n_moons):
        system_indices_res[4-i][j] = ratios.index(greatest_rmrs_ratio_res[4-i][j])


rmrs_ratios_nores, ratios, lines_nores = cal_ratios('nores')

for i in range(5): 
    rmrs_ratios_nores[4-i].sort(reverse=True)
    greatest_rmrs_ratio_nores[4-i] = rmrs_ratios_nores[4-i][:top_n_moons]
    greatest_rmrs_ratio_nores[4-i].sort()

    for j in range(top_n_moons):
        system_indices_nores[4-i][j] = ratios.index(greatest_rmrs_ratio_nores[4-i][j])


for i in range(5):
    for j in range(top_n_moons):
        f.write(lines_res[system_indices_res[4-i][j]].split(' ')[0]+' '+lines_res[system_indices_res[4-i][j]].split(' ')[1]+" res "+lines_res[system_indices_res[4-i][j]].split(' ')[2]+' '+ lines_res[system_indices_res[4-i][j]].split(' ')[3]+"\n")
    for j in range(top_n_moons):
        f.write(lines_nores[system_indices_nores[4-i][j]].split(' ')[0]+' '+lines_nores[system_indices_nores[4-i][j]].split(' ')[1]+" nores "+lines_nores[system_indices_nores[4-i][j]].split(' ')[2]+' '+ lines_nores[system_indices_nores[4-i][j]].split(' ')[3]+"\n")