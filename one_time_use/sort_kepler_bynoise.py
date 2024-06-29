import numpy as np
import os

f = open('/home/garvit/Downloads/Exomoon_project/good_kepler_data.txt', 'r')
lightcurves_files = f.read().split('\n')[:-1]
f.close()

noise_dict={}

for i in range(len(lightcurves_files)):
    #print(i+1)
    f = open("/home/garvit/Downloads/Exomoon_project/good_kepler_data/{}".format(lightcurves_files[i]), 'r')
    kepler_data = f.readlines()
    f.close()

    N = len(kepler_data) -1
    kepler_flux = np.ndarray((N))
    timestamps = np.ndarray((N))
    flux_error = np.ndarray((N))

    for j in range(N):
        timestamps[j] = float(kepler_data[j+1].split('\t')[0])
        kepler_flux[j] = float(kepler_data[j+1].split('\t')[1])
        flux_error[j] = float(kepler_data[j+1].split('\t')[2])

        # if(j!=0): print(timestamps[j] - timestamps[j-1], 1/48)

    noise_arr = flux_error/kepler_flux

    noise_dict[lightcurves_files[i]] = np.mean(noise_arr[np.logical_not(np.isnan(noise_arr))])


new_list = list(dict(sorted(noise_dict.items(), key=lambda item: item[1])).keys())

# print(new_list)

#f = open('/home/garvit/Downloads/Exomoon_project/good_kepler_data_best.txt', 'w')
f = open('/home/garvit/Downloads/Exomoon_project/temp.txt', 'w')
for i in range(len(new_list)):
    f.write(new_list[i]+'\n')
    print(noise_dict[new_list[i]])
f.close()




##### BELOW CODE WAS USED TO FIND THE BEST KEPLER DATA FROM THE 3 FOLDERS THAT I HAD

# f = open('/home/garvit/Downloads/Exomoon_project/good_kepler_data3.txt', 'r')
# lightcurves_files = f.read().split('\n')[:-1]
# f.close()

# lc_list=[]

# for i in range(len(lightcurves_files)):
#     #print(i+1)
#     f = open("/home/garvit/Downloads/Exomoon_project/no_planet_lightcurves3/{}".format(lightcurves_files[i]), 'r')
#     kepler_data = f.readlines()
#     f.close()

#     N = len(kepler_data) -1
#     kepler_flux = np.ndarray((N))
#     timestamps = np.ndarray((N))
#     flux_error = np.ndarray((N))

#     for j in range(N):
#         timestamps[j] = float(kepler_data[j+1].split('\t')[0])
#         kepler_flux[j] = float(kepler_data[j+1].split('\t')[1])
#         flux_error[j] = float(kepler_data[j+1].split('\t')[2])

#     noise_arr = flux_error/kepler_flux
#     noise = np.mean(noise_arr[np.logical_not(np.isnan(noise_arr))])

#     if(noise < 1e-4):
#         digits = [c for c in lightcurves_files[i] if c.isdigit()]
#         kid = ''.join(digits)
#         os.system('cp /home/garvit/Downloads/Exomoon_project/no_planet_lightcurves3/{} /home/garvit/Downloads/Exomoon_project/no_planet_lightcurves_best/{}.txt'.format(lightcurves_files[i], kid)) 

#         print(noise)
        #lc_list.append(kid)



# f = open('/home/garvit/Downloads/Exomoon_project/good_kepler_data_best.txt', 'a')
# for i in range(len(lc_list)):
#     f.write(lc_list[i]+'.txt\n')
# f.close()


