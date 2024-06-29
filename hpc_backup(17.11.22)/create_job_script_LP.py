import sys
import os

systemno = int(sys.argv[1])
modeling_type = str(sys.argv[2])
frac = float(sys.argv[3])
if(len(sys.argv) == 5 and str(sys.argv[4]) == 'resume'): resume = True
else: resume = False

index = systemno -1
catalogue = open('/tiara/home/gagarwal/files_and_data/catalogue.txt', 'r').read().split('\n')
system_line = catalogue[index].split('\t')
system = system_line[0]
system_split = system.split('_')
simno = system_split[-1]

lines = open("/tiara/home/gagarwal/template_job_script.sh", 'r').readlines()

job_line = lines[2].split(' ')
job_name = simno +'_LP_'+ modeling_type  + '_' + job_line[-1]
job_line[-1] = job_name
lines[2] = ' '.join(job_line)

outline = lines[4].split('/')
outline[-1] = "output_"+ simno + '_LP_' + system_split[0]+'_'+system_split[1] +"_"+modeling_type
lines[4] = '/'.join(outline)+'\n'

errline = lines[5].split('/')
errline[-1] = "error_"+ simno + '_LP_' + system_split[0]+'_'+system_split[1] +"_"+modeling_type
lines[5] = '/'.join(errline)+'\n'


if(resume):
	lines[21] = 'python /tiara/home/gagarwal/code/modeling_LP.py {} {} {} resume\n'.format(systemno, modeling_type, frac)
else:
	lines[21] = 'python /tiara/home/gagarwal/code/modeling_LP.py {} {} {}\n'.format(systemno, modeling_type, frac)

f = open('/tiara/home/gagarwal/job_script.sh', 'w')
f.write(''.join(lines))


