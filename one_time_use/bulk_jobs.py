import os
import sys
import time

start = int(sys.argv[1])
end= int(sys.argv[2])

fracs = []
for line in open('/tiara/home/gagarwal/files_and_data/fraction&bad_detrending.txt','r').readlines()[start-1:end]: fracs.append(float(line.split(' ')[4]))

j=0
for i in range(start, end+1):
    os.system('python create_job_script_SNR.py {} planet_only {}'.format(i, fracs[j]))
    time.sleep(2)

    os.system('qsub /tiara/home/gagarwal/job_script.sh')
    time.sleep(5)

    os.system('python create_job_script_SNR.py {} planet_moon {}'.format(i, fracs[j]))
    time.sleep(2)

    os.system('qsub /tiara/home/gagarwal/job_script.sh')
    time.sleep(5)

    j+=1