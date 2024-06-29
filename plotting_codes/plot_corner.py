import numpy as np
import corner
import json
import math
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from astropy.constants import M_earth, R_earth, au
from PIL import Image
import os

def system_index(system):
    for i in range(len(catalogue)):
        if(catalogue[i].find(system+'\t')!=-1): return i

systems = open('/home/garvit/Downloads/Exomoon_project/completed_systems.txt', 'r').read().split('\n')[:-1]
catalogue = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'r').read().split('\n')

for k in range(95,len(systems)):
    model_types = ['planet_only', 'planet_moon']
    for model_type in model_types:
        index = int(system_index(systems[k]))

        system_line = catalogue[index].split('\t')

        system = system_line[0]
        print(system)

        systems_details = json.load(open('/home/garvit/Downloads/Exomoon_project/sim_data/final_unpacked_details/{}.json'.format(system), 'r'))
        w_bary = float(open('/home/garvit/Downloads/Exomoon_project/most_detectable_systems_final.txt', 'r').readlines()[index-2].split(' ')[3])
        omega_moon = float(open('/home/garvit/Downloads/Exomoon_project/most_detectable_systems_final.txt', 'r').readlines()[index-2].split(' ')[4])

        if(model_type=='planet_only'): 
            num_moons = 0
            true_params_planet = np.array([systems_details['Planet']['Pp'], systems_details['Planet']['a']/systems_details['Planet']['Rstar'], systems_details['Planet']['RpRstar'], systems_details['Planet']['impact'], 0.0, w_bary, systems_details['Planet']['q1'], systems_details['Planet']['q2']])

            paramnames = ['$P_P$', '$a_{bary}$', '$R_P$', '$b_{bary}$', '$e_P$', '$w_{bary}$','q1', 'q2']

        else: 
            num_moons = len(systems_details.keys()) - 1
            true_params_planet = np.array([systems_details['Planet']['Pp'], systems_details['Planet']['a']/systems_details['Planet']['Rstar'], systems_details['Planet']['RpRstar'], systems_details['Planet']['impact'], 0.0,  w_bary, systems_details['Planet']['m'], systems_details['Planet']['q1'], systems_details['Planet']['q2']])

            true_params_moons = np.ndarray((num_moons, 15))
            keys = ['I', 'II', 'III', 'IV', 'V']
            for i in range(num_moons):
                #true_params_moons[i] = systems_details[keys[i]]['r']/systems_details['Planet']['Rstar'],systems_details[keys[i]]['P']/(24*3600), (2*np.pi-systems_details[keys[i]]['f'])/(2*np.pi),omega_moon, (math.degrees(systems_details[keys[i]]['inc']) + 90.0)%360 , math.degrees(systems_details[keys[i]]['pomega']), systems_details[keys[i]]['m']
                
                true_params_moons[i] = np.array([systems_details['Planet']['Pp'], systems_details['Planet']['a']/systems_details['Planet']['Rstar'], systems_details['Planet']['RpRstar'], systems_details['Planet']['impact'],0.0, w_bary, systems_details['Planet']['m'],systems_details[keys[i]]['r']/systems_details['Planet']['Rstar'],systems_details[keys[i]]['P']/(24*3600), (2*np.pi-systems_details[keys[i]]['f'])/(2*np.pi),omega_moon, (math.degrees(systems_details[keys[i]]['inc']) + 90.0)%360, systems_details[keys[i]]['m'], systems_details['Planet']['q1'], systems_details['Planet']['q2']])

            paramnames = ['$P_P$', '$a_{bary}$', '$R_P$', '$b_{bary}$', '$e_P$', '$w_{bary}$', '$M_P$', '$R_M$', '$P_M$', '$tau_M$', '$Omega_M$', '$inc_M$', '$M_M$', 'q1', 'q2']



        model_fit_params = json.load(open('/media/Datas/Exomoon_project/modeling_results/{}/{}/results.json'.format(system, model_type), 'r'))
        # paramnames = model_fit_params['paramnames'].copy()
        # paramnames.remove('t0_bary_offset')

        maxl_params = np.delete(np.array(model_fit_params["maximum_likelihood"]["point"]), 6)


        f = open('/media/Datas/Exomoon_project/modeling_results/{}/{}/equal_weighted_post.txt'.format(system, model_type), 'r')
        samples = np.loadtxt(f, skiprows=1)
        cols = [i for i in range(len(model_fit_params['paramnames']))]
        cols.remove(6)
        samples = samples[:, cols]
        samples[:,1]*= (systems_details['Planet']['Rstar']/au.value)
        samples[:,2]*= (systems_details['Planet']['Rstar']/R_earth.value)
        if(model_type=='planet_moon'):
            samples[:,6]/= M_earth.value
            samples[:,12]/= M_earth.value
            samples[:,7]*= (systems_details['Planet']['Rstar']/R_earth.value)

        figure = corner.corner(samples, labels=paramnames, show_titles=True)

        # corner.overplot_lines(figure, maxl_params, color="red")
        # corner.overplot_points(figure, maxl_params[None], marker="s", color="red")

        patches = []
        # patches.append(mpatches.Patch(color='red', label='maximum likelihood parameters'))

        if(model_type=='planet_only'):
            # corner.overplot_lines(figure, true_params_planet, color="blue")
            corner.overplot_points(figure, true_params_planet[None], marker="s", color="blue")

            # # Loop over the diagonal
            ndim = 8
            axes = np.array(figure.axes).reshape((ndim, ndim))
            for i in range(ndim):
                ax = axes[i, i]
                
                # ax.axvline(maxl_params[i], color='red')
                ax.axvline(true_params_planet[i], color="blue")
                
            patches.append(mpatches.Patch(color='blue', label='true parameters: planet'))

            # corner.overplot_lines(figure, maxl_params, color="red")


        else:
            ndim= 15
            colors = ['magenta', 'green', 'orange', 'pink', 'cyan']

            for i in range(num_moons):
                corner.overplot_points(figure, true_params_moons[i][None], marker="s", color=colors[i])
                patches.append(mpatches.Patch(color=colors[i], label='true parameters: moon {}'.format(i+1)))

            axes = np.array(figure.axes).reshape((ndim, ndim))

            planet_indices = [0,1,2,3,4,5,6,13,14]
            # # Loop over the diagonal
            for i in range(ndim):
                ax = axes[i, i]
                
                # ax.axvline(maxl_params[i], color='red')
                for j in range(num_moons):
                    ax.axvline(true_params_moons[j][i], color=colors[j])

                if(i in planet_indices):
                    if(i ==13 or i==14):ax.axvline(true_params_planet[i-6], color="blue")
                    else: ax.axvline(true_params_planet[i], color="blue")

            # # Loop over the histograms
            for y in range(ndim):
                for x in range(y):
                    ax = axes[y, x]

                    xindex = x
                    yindex =y
                    
                    if(xindex ==13 or xindex==14): xindex -=6
                    if(yindex ==13 or yindex==14): yindex -=6

                    if(y in planet_indices and x in planet_indices): ax.plot(true_params_planet[xindex], true_params_planet[yindex], 'sb')

                    # ax.tick_params(which=u'both', color=[1,1,1,0.5])

            patches.append(mpatches.Patch(color='blue', label='true parameters: planet'))

        figure.legend(handles=patches, loc='upper center')
        plt.savefig('/media/Datas/Exomoon_project/modeling_results/corner_plots/{}_{}.png'.format(system, model_type))
        plt.close()
        # im = Image.open('/media/Datas/Exomoon_project/modeling_results/corner_plots/{}_{}.png'.format(system, model_type))
        # im.show()

