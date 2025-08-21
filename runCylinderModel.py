#!/usr/bin/env python

#Andrew Burt - a.burt.12@ucl.ac.uk

import argparse
import os

def runMe(cloud,dmin_range,dmin_step,rcov_diff,nmin,lcyl,runs,workers,matlab_path):

	###
	dmin0 = 0.2
	rcov0 = dmin0+dmin0*rcov_diff
	nmin0 = 4
	noground = 1
	###
	model_names = []
	cloud_name = cloud.split('/')[len(cloud.split('/'))-1].split('.')[0]
	dmin = dmin_range[0]
	while(dmin < dmin_range[1] + (dmin_step/2)):
		rcov = dmin + dmin * rcov_diff
		for i in range(runs):
			m_name = cloud_name+'-'+str(dmin0)+'-'+str(rcov0)+'-'+str(nmin0)+'-'+str(dmin)+'-'+str(rcov)+'-'+str(nmin)+'-'+str(lcyl)+'-'+str(noground)+'-'+str(i)+'.mat'
			if(os.path.isfile(m_name) == False):
				model_names.append(m_name)
		dmin += dmin_step
	###
	if(len(model_names) > 0):
		models_string = '{'
		for i in range(len(model_names)):
			if(i != len(model_names)-1):
				models_string += "'" + model_names[i] + "'" + ','
			else:
				models_string += "'" + model_names[i] + "'" + '}'
		matlab_run_string = matlab_path + ' -nodisplay -r "runCylinderModel('+"'"+cloud+"',"+models_string+","+str(workers)+')"'
		print(matlab_run_string)
		os.system(matlab_run_string)
	###

if __name__ == "__main__":

        parser = argparse.ArgumentParser()
        parser.add_argument('-i','--clouds',nargs='*',default=False,help='ASCII xyx cloud')
        parser.add_argument('-p','--multiprocess',default=1,help='thread count')
        parser.add_argument('-m','--matlab_path',default='/usr/local/MATLAB/R2021b/bin/matlab',help='path to MATLAB binary')
        args = parser.parse_args()
        for i in range(len(args.clouds)):
             runMe(args.clouds[i],[0.035,0.095],0.005,0.1,2,5,10,args.multiprocess,args.matlab_path)
#		runMe(args.clouds[i],[0.045,0.055],0.005,0.1,2,5,1,args.multiprocess,args.matlab_path)
