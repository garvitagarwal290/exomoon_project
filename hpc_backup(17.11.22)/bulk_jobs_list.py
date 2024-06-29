import os
import sys
import time

system_nos = [3, 4, 5, 8, 11, 13, 14, 16, 17, 18, 21, 22, 25, 26, 30, 36, 46, 57, 6, 7, 9, 10, 12, 15, 20, 44, 49]

fracs = []
for line in open('/tiara/home/gagarwal/files_and_data/excl_regions.txt','r').readlines()[1:101]: fracs.append(float(line.split(',')[2]))


for i in range(len(system_nos)):
	os.system('python create_job_script.py {} planet_only {}'.format(system_nos[i], fracs[system_nos[i]-3]))
	time.sleep(2)

	os.system('qsub /tiara/home/gagarwal/job_script.sh')
	time.sleep(5)

	os.system('python create_job_script.py {} planet_moon {}'.format(system_nos[i], fracs[system_nos[i]-3]))
	time.sleep(2)

	os.system('qsub /tiara/home/gagarwal/job_script.sh')
	time.sleep(5)



