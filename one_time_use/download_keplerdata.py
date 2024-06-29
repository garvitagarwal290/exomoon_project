from moonpy import * 
import numpy as np
import os

# swp = np.genfromtxt('/home/garvit/Downloads/Exomoon_project/stars_without_planets_jun2022.txt')
# kics, mags = swp.T[0], swp.T[2]
# pct1 = np.nanpercentile(mags, 1) #### returns 12.381
# best_mag_idxs = np.where(mags <= pct1)[0]  #### returns 19190 star indices
# best_stars = kics[best_mag_idxs].astype(int)
# random_kics = np.random.choice(best_stars, size=700) 

# for i in range(random_kics.shape[0]):
#     print(random_kics[i])






# kics = open('/home/garvit/Downloads/Exomoon_project/one_time_use/new_kepler_ids3.txt', 'r').read().split('\n')

# for i in range(451,len(kics)):
#     try:
#         kic_lc = MoonpyLC(targetID='KIC'+str(kics[i]))
#     except (IndexError, RecursionError, EOFError):
#         continue


# os.system('ls /home/garvit/Downloads/MoonPy-master/Central_Data/Kepler_lightcurves/ > temp')

# folders = open('temp', 'r').read().split('\n')

# for i in range(len(folders)):
#     os.system('cp /home/garvit/Downloads/MoonPy-master/Central_Data/Kepler_lightcurves/{}/*.tsv /home/garvit/Downloads/Exomoon_project/no_planet_lightcurves3/'.format(folders[i]))