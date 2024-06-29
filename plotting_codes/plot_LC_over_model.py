
import numpy as np
import matplotlib.pyplot as plt


catalogue = open('/home/garvit/Downloads/Exomoon_project/catalogue.txt', 'r').read().split('\n')

sysno = 89
catalogue_line = catalogue[sysno -1].split('\t')
system = catalogue_line[0]
num_moons = int(catalogue_line[2])

filename = 'final_lightcurve_{}.csv'.format(system)
f = open("/home/garvit/Downloads/Exomoon_project/final_lightcurves_SNR/{}".format(filename), 'r')
# f = open("/home/garvit/Downloads/Exomoon_project/final_lightcurves_collect/{}".format(filename), 'r')
injected_lightcurve = f.readlines()
f.close()

N = len(injected_lightcurve) - 2
flux = np.ndarray(N)
timestamps = np.ndarray(N)

for j in range(N):
    line = injected_lightcurve[j+2].split(',')
    timestamps[j] = float(line[0])
    flux[j] = float(line[1])


f = open("/media/Datas/Exomoon_project/shifted_finalmodelLCs/{}.txt".format(system), 'r')
model_data = f.readlines()
f.close()

flux_planet_allmoons = np.ndarray((len(model_data) - 1))
time = np.ndarray((len(model_data)-1))
flux_moons = np.ndarray((num_moons, len(model_data)-1))

for j in range(len(model_data)-1):
    time[j] = float(model_data[j+1].split(',')[0])
    flux_planet_allmoons[j] = float(model_data[j+1].split(',')[1])
    for k in range(num_moons):
        flux_moons[k][j] = float(model_data[j+1].split(',')[2+k])


labels = ['I', 'II', 'III', 'IV', 'V']
# plt.plot(time, flux_planet_allmoons, label='planet')
# for j in range(num_moons):
#     plt.plot(time, flux_moons[j], label=labels[j])
plt.scatter(timestamps, flux, s=5, color='black')
# plt.legend()
plt.show()