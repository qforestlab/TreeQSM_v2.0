#!/usr/bin/env python

#Andrew Burt - a.burt.12@ucl.ac.uk

import argparse
import os
import glob
import numpy
import math
import scipy.io
import scipy.stats
import shutil
import subprocess


#import plotTrees

def getResults(model_names,minradius):

	lid = model_names[0].split('/')[len(model_names[0].split('/'))-1].split('-')[0]
	volume = []
	for i in range(len(model_names)):
		try:
			model = scipy.io.loadmat(model_names[i])
			model_volume = 0
			for j in range(len(model['Rad'])):
				if(model['Rad'][j][0] >= minradius):
					model_volume += math.pi * model['Rad'][j][0] * model['Rad'][j][0] * model['Len'][j][0]
			volume.append(model_volume)
		except:
			continue
	volume = numpy.array(volume)
	mean = numpy.mean(volume)
	stddev = numpy.std(volume)
	idx = numpy.argmin(numpy.abs(volume[:]-mean))
	best_model = model_names[idx]
	results = numpy.zeros(1,dtype=[('lid','S256'),('lv',float),('lv_ru',float)])
	results['lid'][0] = lid
	results['lv'][0] = mean
	relative_uncertainty = 2 * stddev / mean
	results['lv_ru'][0] = relative_uncertainty
	return results,best_model

def getVolume(model_names,minradius):

	volume = []
	for i in range(len(model_names)):
		try:
			model = scipy.io.loadmat(model_names[i])
			model_volume = 0
			for j in range(len(model['Rad'])):
				if(model['Rad'][j][0] >= minradius):
					model_volume += math.pi * model['Rad'][j][0] * model['Rad'][j][0] * model['Len'][j][0]
			volume.append(model_volume)
		except:
			continue
	volume = numpy.array(volume)
	mean = numpy.mean(volume)
	stddev = numpy.std(volume)
#	print stddev
#	print mean
	coeff_variation = stddev / mean
#	print coeff_variation
	return mean,coeff_variation

def getTrunkCylinderFromModel(model,z_coordinate):

	z = 0
	for i in range(len(model['CiB'][0][0])):
		z += model['Len'][(model['CiB'][0][0][i][0]-1)]
		if(z >= z_coordinate):
			cylinder = model['CiB'][0][0][i][0]
			break
	return cylinder

def getCloudFromModel(cloud,model,cylinder):

	cylinder = cylinder-1 # -1 PYTHON CONVENTION
	x_t = model['Sta'][cylinder][0]+(model['Len'][cylinder][0]*model['Axe'][cylinder][0])
	y_t = model['Sta'][cylinder][1]+(model['Len'][cylinder][0]*model['Axe'][cylinder][1])
	z_t = model['Sta'][cylinder][2]+(model['Len'][cylinder][0]*model['Axe'][cylinder][2])
	x_b = model['Sta'][cylinder][0]
	y_b = model['Sta'][cylinder][1]
	z_b = model['Sta'][cylinder][2]
	x_min = min(x_t-model['Rad'][cylinder],x_b-model['Rad'][cylinder])
	x_max = max(x_t+model['Rad'][cylinder],x_b+model['Rad'][cylinder])
	y_min = min(y_t-model['Rad'][cylinder],y_b-model['Rad'][cylinder])
	y_max = max(y_t+model['Rad'][cylinder],y_b+model['Rad'][cylinder])
	z_min = min(z_t,z_b)
	z_max = max(z_t,z_b)
	out = []
	#f1 = open('test.xyz','w')
	for i in range(len(cloud)):
		if(cloud[i][0] >= x_min and cloud[i][0] <= x_max and cloud[i][1] >= y_min and cloud[i][1] <= y_max and cloud[i][2] >= z_min and cloud[i][2] <= z_max):
			out.append([cloud[i][0],cloud[i][1],cloud[i][2]])
			#tmp_out = str(cloud[i][0])+' '+str(cloud[i][1])+' '+str(cloud[i][2])+'\n'
			#f1.write(tmp_out)
	#f1.close()
	cloud = numpy.array(out)
	return cloud

def getDiameterFromCloud(cloud):

	try:
		dbh_pts = cloud
		x=dbh_pts[:,0]
		y=dbh_pts[:,1]
		z=dbh_pts[:,2]
		x_m = x.mean()
		y_m = y.mean()
		u = x - x_m
		v = y - y_m
		Suv = sum(u*v)
		Suu = sum(u**2)
		Svv = sum(v**2)
		Suuv = sum(u**2*v)
		Suvv = sum(u*v**2)
		Suuu = sum(u**3)
		Svvv = sum(v**3)
		A = numpy.array([[ Suu, Suv ],[Suv, Svv]])
		B = numpy.array([Suuu+Suvv,Svvv+Suuv])/2.0
		uc, vc = numpy.linalg.solve(A, B)
		xc_1 = x_m+uc
		yc_1 = y_m+vc
		Ri_1 = numpy.sqrt((x-xc_1)**2 +(y-yc_1)**2)
		r_final = numpy.mean(Ri_1)
		diameter = r_final * 2
		#residu_1 = sum((Ri_1-r_final)**2)
		return diameter,x_m,y_m
	except:
		return numpy.nan,numpy.nan,numpy.nan

def getTrunkCloudModelComparison(cloud_name,model_names):

	TRUNK_POSITIONS = [0.075,0.1,0.125,0.15]
	cloud = numpy.loadtxt(cloud_name)
	tree_diff = []
	for i in range(len(model_names)):
		#this try/catch is extremely stupid - hack for bug in scipy.io.loadmat giving occasional error
		try:
			model = scipy.io.loadmat(model_names[i])
			trunk_length = model['BLen'][0][0]
			position_diff = []
			for j in range(len(TRUNK_POSITIONS)):
				cylinder = getTrunkCylinderFromModel(model,trunk_length*TRUNK_POSITIONS[j]) ### MATLAB CONVENTION
				model_diameter = (model['Rad'][(cylinder-1)][0])*2 ### -1 FOR PYTHON CONVENTION 
				trunk_cloud = getCloudFromModel(cloud,model,cylinder)
				cloud_diameter,cloud_x,cloud_y = getDiameterFromCloud(trunk_cloud)
				diff = min(cloud_diameter,model_diameter)/max(cloud_diameter,model_diameter)
				position_diff.append(diff)
			tree_diff.append(position_diff)
		except:
			continue
	tree_diff = numpy.array(tree_diff)
	result = tree_diff[~numpy.isnan(tree_diff)].mean()
	return result

def optimise(model_id,cloud_dir,model_dir,radiusmin):

	cloud_search = cloud_dir+model_id+'.txt'
	cloud = glob.glob(cloud_search)[0]
	dmin_search = model_dir+model_id+"-*.mat"
	fnames = glob.glob(dmin_search)
	dmin_range = []
	for i in range(len(fnames)):
		dmink = float(fnames[i].split("-")[4])
		dmin_range.append(dmink)
	dmin_range = numpy.array(dmin_range)
	dmin_range = numpy.unique(dmin_range)
	dmin_range = numpy.sort(dmin_range)
	metadata = []
	for j in range(len(dmin_range)):
		model_search = model_dir+model_id+'-*-*-*-'+str(dmin_range[j])+'-*-*-*-*-*.mat'
		models = glob.glob(model_search)
		if(len(models) >= 3):
			mean,coeff_variation = getVolume(models,radiusmin)
			trunk_difference = getTrunkCloudModelComparison(cloud,models)
			metadata.append([dmin_range[j],mean,coeff_variation,trunk_difference])
	metadata = numpy.array(metadata)
	metadata_name = model_id + '.old.opt'
#	numpy.savetxt('opt_param/'+metadata_name,metadata,fmt='%.4f')
	##		
	dmin_opt = numpy.nan
	##
	indices = [i for (i,v) in enumerate(metadata[:,1]) if v==0] #remove values with na - not enough samples to calculate mean?
	metadata=numpy.delete(metadata, indices, 0)
	metadata_name = model_id + '.opt'
	numpy.savetxt('opt_param/'+metadata_name,metadata,fmt='%.4f')
#	for m in range(len(metadata)):
#		if(metadata[m][2] < min_cov * 5 and metadata[m][3] > max_conf * 0.875):
#			dmin_opt = metadata[m][0]
#			break
	min_cov = numpy.min(metadata[:,2])
	max_conf = numpy.max(metadata[:,3])
	for m in range(len(metadata)):
		if(metadata[m][2] < min_cov * 5 and metadata[m][3] > max_conf * 0.95):
			dmin_opt = metadata[m][0]
			break
	if(numpy.isnan(dmin_opt) == True):
		idx = numpy.argmin(metadata[:,2])
		dmin_opt = metadata[idx][0]
	##
	model_search = model_dir+model_id+'-*-*-*-'+str(dmin_opt)+'-*-*-*-*-*.mat'
	models = glob.glob(model_search)
	results,best_model = getResults(models,radiusmin)
	#numpy.savetxt(model_id+'.dat',results,fmt='%s %.3f %.3f')
	#shutil.copy(best_model,model_id+'.mat')
	image_name = model_id+'.pdf'
	#plotTrees.plotCloudsModels([cloud],[model_id+'.mat'],0,0,radiusmin,image_name)
	for i in range(len(models)):
		string1 = 'cp '+models[i]+' opt/'
		string2 = 'cp results/ModelData_'+models[i].split("/")[1]+' opt/'
		string3 = 'cp results/cyl_data_'+models[i].split("/")[1]+'.txt opt/'
		os.system(string1)
		os.system(string2)
		#os.system(string3)
	string1 = 'cp `ls opt/'+model_id+'*mat | head -n 1` opt_1mod'
	os.system(string1)
	string2 = 'cp `ls opt/ModelData_'+model_id+'*mat | head -n 1` opt_1mod'
	os.system(string2)
	string3 = 'cp `ls opt/cyl_data_'+model_id+'*txt | head -n 1` opt_1mod'
	#os.system(string3)

if __name__ == "__main__":

	def func_star(a_b):
	    return optimise(*a_b)

	parser = argparse.ArgumentParser()
	parser.add_argument('-c','--cloud_dir',help='path to cloud directory')
	parser.add_argument('-m','--model_dir',help='path to model directory')
	parser.add_argument('-r','--radius_min',type=float,default=0.025,help='min model rad')
	args = parser.parse_args()
	model_fnames = glob.glob(str(args.model_dir)+"*-*.mat")
	model_id = []
	for i in range(len(model_fnames)):
		mid = model_fnames[i].split('/')[len(model_fnames[i].split('/'))-1].split('-')[0]
		model_id.append(mid)
	model_id = numpy.array(model_id)
	model_id = numpy.unique(model_id)
	cloud_fnames = glob.glob(str(args.cloud_dir)+"*.txt")
	cloud_id = []
	for j in range(len(cloud_fnames)):
		cid = cloud_fnames[j].split('/')[len(cloud_fnames[j].split('/'))-1].split('.')[0]
		print(cid)
		cloud_id.append(cid)
	tid = list(set(cloud_id).intersection(model_id))
	print(tid)
	string = 'mkdir opt'
	os.system(string)
	string = 'mkdir opt_1mod'
	os.system(string)
	string = 'mkdir opt_param'
	os.system(string)
	for k in range(len(tid)):
		print(tid[k])
		tracker = "opt/"+tid[k]+".tracker"
		if not os.path.isfile(tracker):
			subprocess.call(['touch', tracker])
			optimise(tid[k],args.cloud_dir,args.model_dir,args.radius_min)
