# -*- coding: utf-8 -*-
import os 
import shutil 
import numpy

for virus_birth in numpy.arange(0.1, 0.5, 0.1): # from 1 to 5, in steps of 1, please adjust to your liking
		for virus_death in numpy.arange(0.1, 0.5, 0.1): # from 0.1 to 2 in steps of 0.1
			DIRNAME = "data_b_" + str(virus_birth) + "_d_" + str(virus_death)
			os.mkdir(DIRNAME)

			source = 'testmodel.out'
			destination = './' + DIRNAME + '/testmodel.out'
			shutil.copyfile(source, destination)
			os.chdir(DIRNAME)
			os.popen('chmod 777 testmodel.out')    										# this gives administrator rights

			c=open("config.ini",'w')													# the config file!!
			c.write('maximum_time = 1000\n')
			c.write('time_adding_virus = 168\n')										# 7 days of 24 hours = 168 hours, following Berg et al.
			c.write('initial_number_cancer_cells = 500\n')								# following Berg et al.
			c.write('initial_number_normal_cells = 500\n')
			
			c.write('birth_normal = 0.5\n')
			c.write('death_normal = 0.2\n')

			c.write('birth_cancer = 1.0\n')
			c.write('death_cancer = 0.1\n')

			c.write('birth_infected = ' + str(virus_birth) + '\n')
			c.write('death_infected = ' + str(virus_death) + '\n')

			c.write('percent_infected = 0.1\n')											# please verify if Berg et al. also used this value!
			c.write('sq_num_cells = 1000\n')												# for testing purposes, 100 should be better
			c.write('infection_type = center_infection\n')								# please verify with Berg et al that they used this (I think so)
			c.write('start_setup = converge\n')											# this is a new start setup, that follows Berg et al's approach more closely


			c.write('birth_resistant = 0.9\n')											# not used in this version
			c.write('death_resistant = 0.1\n')											# not used in this version
			c.write('prob_normal_infection = 0.0\n')									# This is our own addition, not used by Berg et al.
			c.write('freq_resistant = 0.0\n')											# This is our own addition, not used by Berg et al.
			c.write('distance_infection_upon_death = 0\n')								# This is our own addition, not used by Berg et al.
			c.write('prob_infection_upon_death = 0\n')									# This is our own addition, not used by Berg et al.
			c.write('time_adding_cancer = 1000\n')										# this parameter is not used
			
			
			os.chdir('../')


			# please find below al job script as I would use on our cluster here:

			job_name= 'jobfile_' + str(virus_birth) + '_' + str(virus_death) + '.sh'
			f=open(job_name,'w') 
			f.write('#!/bin/bash\n')
			f.write('#SBATCH -J Range_' + str(virus_birth) + str(virus_death) + '\n')
			f.write('#SBATCH -N 1\n')
			f.write('#SBATCH -n 1\n')
			f.write('#SBATCH --mem 2048\n')
			f.write('#SBATCH -t 3-00:00  # time (D-HH:MM)\n')
			f.write('#SBATCH -o slurm.%j.out # STDOUT\n')
			f.write('#SBATCH -e slurm.%j.err # STDERR\n')
			f.write('module load GCC/8.3.0\n')
			f.write('cd ' + DIRNAME + '\n')
			f.write('./testmodel.out')
			f.close();
