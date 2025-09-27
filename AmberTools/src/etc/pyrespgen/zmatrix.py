#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
from pyrespgen.constant import *
from collections import Counter
from pyrespgen.constant_radii import tbl
#import sys
#sys.setrecursionlimit(99999)

def flatten_list(_2d_list):
	flat_list = []
	# Iterate through the outer list
	for element in _2d_list:
		# if the element is of type list, iterate through the sublist
		if type(element) is list:
			for item in element:
				flat_list.append(item)
		else:
			flat_list.append(element)
	return flat_list

def readconnectTBL(filename='CONNECT.TPL'):
	TBL = {'elem_syl':[],'elem_num':[],'elem_radii':[]}
	with open(filename,'r') as f:
		while True:
			line = f.readline()
			if not line:
				break
			else:
				if '//' not in line[:5]:
					TBL['elem_syl'].append(line.split()[0])
					TBL['elem_num'].append(int(line.split()[1]))
					TBL['elem_radii'].append(float(line.split()[2]))
		else:
			print('ERROR: could not find **CONNECT.TBL**',end='\n')
	
	TBL = pd.DataFrame(data=TBL)
	return TBL


def getINFO(espdat):
	info = {'natoms':[],'ngrids':[]}
	with open(espdat,'r') as f:
		while True:
			line = f.readline()
			if not line:
				break
			else:
				info['natoms'] = int(line[0:5])
				info['ngrids'] = int(line[5:10])
				break
	return info


def getXYZ(espdat):
	xyz = {'ida':[],'syl':[],'x':[],'y':[],'z':[]}
	with open(espdat,'r') as f:
		while True:
			line=f.readline()
			if not line:
				break
			else:
				if len(line[:16].strip()) ==0:
					xyz['ida'].append(int(line.split()[3]))
					xyz['syl'].append(line.split()[4])
					xyz['x'].append(float(line.split()[0])/ang2bohr)
					xyz['y'].append(float(line.split()[1])/ang2bohr)
					xyz['z'].append(float(line.split()[2])/ang2bohr)
	xyz = pd.DataFrame(data=xyz)
	return xyz

def distMatrix(xyz):
	dmat = distance_matrix(xyz[['x','y','z']].values, xyz[['x','y','z']].values)
	#print('dmat type',type(dmat))
	dmat = np.round(dmat,decimals=3)
	dmat = pd.DataFrame(data=dmat, columns=xyz['syl'].values, index=xyz['ida'].values)
	return dmat

#def getBOND(xyz,bcut=1.55,fconnecttbl='./CONNECT.TPL'):
#def getBOND(xyz,fconnecttbl='CONNECT.TPL'): # No need for cutoff again
def getBOND(xyz,dtbl=tbl): # No need for cutoff again
	"""
		get bond number of a certain molecule
		bcut: cutoff in the unit of angstrom
				default set to be 1.94 angstrom
				which resembles the C-Br single bond
		bond length is borrowed from the following website:
		http://www.csb.yale.edu/userguides/datamanip/uppsala/manuals.20020128/typical_bonds.html
	"""
	Nbond = 0
	#Nbond_new = 0      # for testing the new function
	dmat = distMatrix(xyz)
	# read the atom radii from CONNECT.TPL
	#tpl = readconnectTBL(fconnecttbl)
	tpl = pd.DataFrame(data=dtbl,columns=['elem_syl','elem_num','elem_radii'])
	ref_range = [0.0, 1.5, 1.9, 2.05]
	bnl_offset_fact = [0.15, 0.11, 0.09, 0.08]

	a,b = dmat.shape[0], dmat.shape[1]
	#print(a,b)
	for i in range(a):
		for j in range(b):
			if i > j:
				## the following logic was copied directly from Yong Duan's pol_resp code
				ref = tpl.loc[tpl['elem_num']==dmat.index[i]]['elem_radii'].values[0] + tpl.loc[tpl['elem_num']==dmat.index[j]]['elem_radii'].values[0]
				for num_i in range(len(ref_range)):
					if num_i < len(ref_range) -1:
						if ref > ref_range[num_i] and ref <= ref_range[num_i+1]:
							offsetfact = bnl_offset_fact[num_i]
							offset = ref*offsetfact
							if dmat.iloc[i][j] < ref + offset and dmat.iloc[i][j] > ref * 0.5:
								#print('Atom {} - Atom {}'.format(i,j))
								Nbond+=1
					else:
						if ref > ref_range[num_i]:
							offsetfact = bnl_offset_fact[num_i]
							offset = ref*offsetfact
							if dmat.iloc[i][j] < ref + offset and dmat.iloc[i][j] > ref * 0.5:
								#print('Atom {} - Atom {}'.format(i,j))
								Nbond+=1
				
			#	# criterion for C-C single bond
			#	if dmat.index[i] ==6 and dmat.index[j] == 6:
			#			if dmat.iloc[i][j] < 1.55:
			#				Nbond +=1
			#	# criterion for N-N single bond
			#	elif dmat.index[i] ==7 and dmat.index[j] == 7:
			#			if dmat.iloc[i][j] < 1.36:
			#				Nbond +=1
			#	# criterion for S-S single bond
			#	elif dmat.index[i] == 16 and dmat.index[j] == 16:
			#			if dmat.iloc[i][j] < 2.08:
			#				Nbond +=1
			#	# criterion for C-O or O-C single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 8) or (dmat.index[i] ==8 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.44:
			#				Nbond +=1
			#	# criterion for C-N or N-C single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 7) or (dmat.index[i] ==7 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.48:
			#				Nbond +=1
			#	# criterion for X-H, where X is not H
			#	elif (dmat.index[i] !=1 and dmat.index[j] == 1) or (dmat.index[i] ==1 and dmat.index[j] != 1) :
			#			if dmat.iloc[i][j] < 1.36:
			#				Nbond +=1
			#	# criterion for O-N or N-O single bond
			#	elif (dmat.index[i] ==8 and dmat.index[j] == 7) or (dmat.index[i] ==7 and dmat.index[j] == 8) :
			#			if dmat.iloc[i][j] < 1.25:
			#				Nbond +=1
			#	# criterion for O-P or P-O single bond
			#	elif (dmat.index[i] ==8 and dmat.index[j] == 15) or (dmat.index[i] ==15 and dmat.index[j] == 8) :
			#			if dmat.iloc[i][j] < 1.64:
			#				Nbond +=1
			#	# criterion for C-Br single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 35) or (dmat.index[i] ==35 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.95:
			#				Nbond +=1
			#	# criterion for C-F single bond 
			#	#           the length threshold was increased from 1.38 to 1.39
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 9) or (dmat.index[i] ==9 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.39:
			#				Nbond +=1
			#	# criterion for C-Cl single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 17) or (dmat.index[i] ==17 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.77:
			#				Nbond +=1
			#	# criterion for C-I single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 53) or (dmat.index[i] ==53 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 2.15:
			#				Nbond +=1
			#	# criterion for C-P single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 15) or (dmat.index[i] ==15 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.85:
			#				Nbond +=1
			#	# criterion for C-S single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 16) or (dmat.index[i] ==16 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.83:
			#				Nbond +=1
			#	# for missing cases!
			#	else:
			#		if dmat.iloc[i][j] < bcut:
			#			Nbond+=1
	#print('Old Bond Num: {}, New: {}'.format(Nbond, Nbond_new))
	return Nbond


#def getCONN(xyz,bcut=1.55,fconnecttbl='./CONNECT.TPL'):
def getCONN(xyz,dtbl=tbl): # No need for cutoff again
	"""
		get bond number of a certain molecule
		bcut: cutoff in the unit of angstrom
				default set to be 1.94 angstrom
				which resembles the C-Br single bond
		bond length is borrowed from the following website:
		http://www.csb.yale.edu/userguides/datamanip/uppsala/manuals.20020128/typical_bonds.html
	"""
	# initializing connect map
	conIdx = {}
	connect = {}	# storing connected atom index in the esp data file, such as 0,1,2,3,4, etc.
	mapAtmIdx = {}
	for i,idx in enumerate(xyz['ida'].values):
		connect[i] = []
		conIdx[i] = []
	# mapping atomtype to index
	conAtm = {} 	# storing connected atom type such as nh, hn, c3 etc.
	mapAtmType = {}	# mapping index to the atom symbol or the atom type utilized in the amber force field
	for i,sym in enumerate(xyz['syl'].values):
		mapAtmType[i] = sym
		conAtm[i] = []

	for syl,idx in zip(xyz['syl'].values,xyz['ida'].values):
		mapAtmIdx[syl] = idx	# mapping symbol to the atom number, such as hn -> 1 | c3 -> 6
	#print(connect)
	# constructing the distance matrix
	dmat = distMatrix(xyz)
	a,b = dmat.shape[0], dmat.shape[1]

	# read the atom radii from CONNECT.TPL
	#tpl = readconnectTBL(fconnecttbl)
	tpl = pd.DataFrame(data=dtbl,columns=['elem_syl','elem_num','elem_radii'])
	ref_range = [0.0, 1.5, 1.9, 2.05]
	bnl_offset_fact = [0.15, 0.11, 0.09, 0.08]
	
	#print(a,b)
	for i in range(a):
		for j in range(b):
			if i != j:
				## the following logic was copied directly from Yong Duan's pol_resp code
				ref = tpl.loc[tpl['elem_num']==dmat.index[i]]['elem_radii'].values[0] + tpl.loc[tpl['elem_num']==dmat.index[j]]['elem_radii'].values[0]
				for num_i in range(len(ref_range)):
					if num_i < len(ref_range)-1:
						if ref > ref_range[num_i] and ref <= ref_range[num_i+1]:
							offsetfact = bnl_offset_fact[num_i]
							offset = ref*offsetfact
							if dmat.iloc[i][j] < ref + offset and dmat.iloc[i][j] > ref * 0.5:
								connect[i].append(j)
								conAtm[i].append(mapAtmType[j])
								conIdx[i].append(mapAtmIdx[mapAtmType[j]])
					else:
						if ref > ref_range[num_i]:
							offsetfact = bnl_offset_fact[num_i]
							offset = ref*offsetfact
							if dmat.iloc[i][j] < ref + offset and dmat.iloc[i][j] > ref * 0.5:
								connect[i].append(j)
								conAtm[i].append(mapAtmType[j])
								conIdx[i].append(mapAtmIdx[mapAtmType[j]])
							
							#Nbond+=1
			#	# criterion for C-C single bond
			#	if dmat.index[i] ==6 and dmat.index[j] == 6:
			#			if dmat.iloc[i][j] < 1.55:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#	# criterion for N-N single bond
			#	elif dmat.index[i] ==7 and dmat.index[j] == 7:
			#			if dmat.iloc[i][j] < 1.36:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# criterion for S-S single bond
			#	elif dmat.index[i] == 16 and dmat.index[j] == 16:
			#			if dmat.iloc[i][j] < 2.08:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# criterion for C-O or O-C single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 8) or (dmat.index[i] ==8 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.44:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# criterion for C-N or N-C single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 7) or (dmat.index[i] ==7 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.48:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# criterion for X-H, where X is not H
			#	elif (dmat.index[i] !=1 and dmat.index[j] == 1) or (dmat.index[i] ==1 and dmat.index[j] != 1) :
			#			if dmat.iloc[i][j] < 1.36:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# criterion for O-N or N-O single bond
			#	elif (dmat.index[i] ==8 and dmat.index[j] == 7) or (dmat.index[i] ==7 and dmat.index[j] == 8) :
			#			if dmat.iloc[i][j] < 1.25:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# criterion for O-P or P-O single bond
			#	elif (dmat.index[i] ==8 and dmat.index[j] == 15) or (dmat.index[i] ==15 and dmat.index[j] == 8) :
			#			if dmat.iloc[i][j] < 1.64:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# criterion for C-Br single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 35) or (dmat.index[i] ==35 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.95:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# criterion for C-F single bond
			#	#              the bond length increased from 1.38 to 1.39
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 9) or (dmat.index[i] ==9 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.39:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# criterion for C-Cl single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 17) or (dmat.index[i] ==17 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.77:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# criterion for C-I single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 53) or (dmat.index[i] ==53 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 2.15:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# criterion for C-P single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 15) or (dmat.index[i] ==15 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.85:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# criterion for C-S single bond
			#	elif (dmat.index[i] ==6 and dmat.index[j] == 16) or (dmat.index[i] ==16 and dmat.index[j] == 6) :
			#			if dmat.iloc[i][j] < 1.83:
			#				connect[i].append(j)
			#				conAtm[i].append(mapAtmType[j])
			#				conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#				#conIdx[i].append(mapAtmIdx[j])
			#	# for missing cases!
			#	else:
			#		if dmat.iloc[i][j] < bcut:
			#			connect[i].append(j)
			#			conAtm[i].append(mapAtmType[j])
			#			conIdx[i].append(mapAtmIdx[mapAtmType[j]])
			#			#conIdx[i].append(mapAtmIdx[j])
	return connect,conAtm,conIdx

def getEQU(connect,xyz,depth=2):
	"""
		input results extracted from func getCONN
			- connect;
			- conAtm;
			- xyz
		depth: specifying how deep the search should be
	"""
	connectAtm = {}
	# constructing tree depth
	for i in range(1,depth):
		connect_tmp = {}
		for i in connect.keys():
			tmp = [connect[x] for x in connect[i]]
			tmp_flatten = flatten_list(tmp)
			connect_tmp[i] = [item for item in tmp_flatten if item !=i] # trim the center atom from list
		connect = connect_tmp
	
	# count the keys bonded to the center atom with specified depth
	for x in connect.keys():
		connectAtm[x] = Counter([xyz['syl'][idx] for idx in connect[x]])
	
	# build map to find equilvalent atom
	equAtm = np.zeros((len(xyz['syl'].values),len(xyz['syl'].values))) 

	for i in range(len(xyz['syl'])):
		for j in range(len(xyz['syl'])):
			# if two atoms are equivalent, two criteria should be satisfied:
			#           - these two atom type defined in the amber force field should be the same
			#           - the connected atom in the defined depth (dafault set to be 3) should be the same
			if connectAtm[i] == connectAtm[j] and xyz['syl'][i] == xyz['syl'][j]: 
				equAtm[i][j] = 1.0
			else:
				equAtm[i][j] = 0.0
	
	equAtm = pd.DataFrame(data=equAtm, index= xyz['ida'].values, columns=xyz['ida'].values)

	return equAtm

def getCH2CH3(conIdx,xyz):
	'''
		find index of -CH2- (methylene)  and -CH3 (methyl)
	'''
	idx_ch3 = []
	idx_ch2 = []
	tpl_ch3 = [6,1,1,1]
	tpl_ch2 = [6,6,1,1]
	for key,values in conIdx.items():
		if len(values) ==4 and Counter(values)[1]==3 and xyz['ida'][key]==6:
			idx_ch3.append(key)
		elif len(values) ==4 and Counter(values)[1]==2 and xyz['ida'][key]==6:
			idx_ch2.append(key)
		elif len(values) ==4 and Counter(values)[1]==4 and xyz['ida'][key]==6:
			# condition for CH4
			idx_ch3.append(key)

	return idx_ch2,idx_ch3

def getXHX(conIdx, xyz):
	"""
		find index of EQUILVALENT HYDROGEN function groups except -CH2-(methylene) and -CH3(methyl)
	"""
	pass




def getSymIdx(equAtm,xyz,idx_ch2,idx_ch3,conn):
	"""
		get the equilvalent atom index:
			- symidx_hvy: heavy atom index which are equilvalent
			- symidx_ch2: hydrogen equilvalent index for fragment -CH2-
			- symidx_ch3: hydrogen equilvalent index for fragment -CH3
			- symidx_hxx: hydrogen equilvalent index for function groups except -CH2- and -CH3
	"""
	symidx_hvy = []
	symidx_ch2 = []
	symidx_ch3 = []
	symidx_hxx = []
	for i in range(equAtm.shape[0]):
		for j in range(equAtm.shape[1]):
			if i < j and equAtm.values[i][j] == 1.0:
				#print('{} | {} | {}'.format(i,j,equAtm.values[i][j]))
			#if equAtm.values[i][j] == 1.0: # enumerate all
				# 1st stage: find equivalent heavy atoms including the carbon in -CH3 and -CH2-
				if xyz['ida'][i] !=1:
					if j not in flatten_list(symidx_hvy):
						symidx_hvy.append([i,j])
				# find equivalent hydrogen atoms for -CH2-, -CH3 and others
				else:
					if len(idx_ch2) !=0 and len(idx_ch3) !=0:
						for ich2 in idx_ch2:
							if i in conn[ich2]:
								if j not in flatten_list(symidx_ch2):
									symidx_ch2.append([i,j])
							elif i not in flatten_list([conn[ich2_tmp] for ich2_tmp in idx_ch2+idx_ch3]):
								#print(i,flatten_list([conn[ich2_tmp] for ich2_tmp in idx_ch2+idx_ch3]))
								# make sure i not in -CH2- connect list
								if j not in flatten_list(symidx_hxx):
									symidx_hxx.append([i,j])
							#else:
							#	if j not in flatten_list(symidx_hxx):
							#		symidx_hxx.append([i,j])
						for ich3 in idx_ch3:
							if i in conn[ich3]:
								if j not in flatten_list(symidx_ch3):
									symidx_ch3.append([i,j])
							elif i not in flatten_list([conn[ich3_tmp] for ich3_tmp in idx_ch3+idx_ch2]):
								if j not in flatten_list(symidx_hxx):
									symidx_hxx.append([i,j])
							#else:
							#	if j not in flatten_list(symidx_hxx):
							#		symidx_hxx.append([i,j])
					elif len(idx_ch2) == 0 and len(idx_ch3) !=0:
						for ich3 in idx_ch3:
							if i in conn[ich3]:
								if j not in flatten_list(symidx_ch3):
									symidx_ch3.append([i,j])
							elif i not in flatten_list([conn[ich3_tmp] for ich3_tmp in idx_ch3]):
								#print(i,conn[ich3])
								if j not in flatten_list(symidx_hxx):
									symidx_hxx.append([i,j])
							#else:
							#	if j not in flatten_list(symidx_hxx):
							#		symidx_hxx.append([i,j])
					elif len(idx_ch2) !=0 and len(idx_ch3) ==0:
						for ich2 in idx_ch2:
							if i in conn[ich2]:
								if j not in flatten_list(symidx_ch2):
									symidx_ch2.append([i,j])
							elif i not in flatten_list([conn[ich2_tmp] for ich2_tmp in idx_ch2]):
								if j not in flatten_list(symidx_hxx):
									symidx_hxx.append([i,j])
							#else:
							#	if j not in flatten_list(symidx_hxx):
							#		symidx_hxx.append([i,j])
					else:
						if j not in flatten_list(symidx_hxx):
							symidx_hxx.append([i,j])
	return symidx_hvy,symidx_ch2,symidx_ch3,symidx_hxx

def dumpSymIdx(symidx_hvy, symidx_ch2, symidx_ch3, symidx_hxx, idx_ch2, idx_ch3, xyz, conn,ptype='perm',strategy=2):
	in1_idx = np.zeros(len(xyz['ida']))
	in2_idx = [-1 for x in range(len(xyz['ida']))]
	in1_idx_dipole = [0 for k,v in conn.items() for x in v]
	in2_idx_dipole = [-1 for k,v in conn.items() for x in v]
	if ptype=='chg' or ptype=='ind':
		# define the equivalent index in the 1st file 
		#		for hydrogen atoms not involved in -CH2- or -CH3
		for h in symidx_hxx:
			in1_idx[h[0]] = 0
			in1_idx[h[1]] = h[0]+1
			# treating the whole function group to change together
			#	such as the NH2, which not only put the equivalent label to H
			#	but also let charge of atom N to chage freely.
			for k,v in conn.items():
				if h[0] in v and h[1] in v: # if two hydrogen are in the connected, we guess the center atom is also in the function group
					in1_idx[k] =0
		## define the equivalent index in the 2nd file ; remain the label to be -1 to be fixed
		##		for only hydrogen atoms
		#for h in symidx_hxx:
		#	in2_idx[h[0]] = 0
		#	in2_idx[h[1]] = h[0]+1
		#	# treating the whole function group to change together
		#	#	such as the NH2, which not only put the equivalent label to H
		#	#	but also let charge of atom N to chage freely.
		#	for k,v in conn.items():
		#		if h[0] in v and h[1] in v: # if two hydrogen are in the connected, we guess the center atom is also in the function group
		#			in2_idx[k] =0
		if len(idx_ch2) !=0:
			# let charge of atom C to change freely
			for c_ch2 in idx_ch2:
				in2_idx[c_ch2] = 0
			for h_ch2 in symidx_ch2:
			# let hydrogen of -CH2- to maintance the equivalent
				in2_idx[h_ch2[0]] = 0
				in2_idx[h_ch2[1]] = h_ch2[0]+1
		if len(idx_ch3) !=0:
			# let charge of atom C to change freely
			for c_ch3 in idx_ch3:
				in2_idx[c_ch3] = 0
			# let hydrogen of -CH3 to maintance the equivalent
			for h_ch3 in symidx_ch3:
				in2_idx[h_ch3[0]] = 0
				in2_idx[h_ch3[1]] = h_ch3[0]+1

		# define the equivalent index in the 1st file
		#		for heavy atoms
		for hvy in symidx_hvy:
			in1_idx[hvy[1]] = hvy[0]+1
	
		# set the equivalent index in the 2nd file,  remain the label to be -1 to be fixed
		#			for heavy atoms
		for hvy in symidx_hvy:
			#print('TEST',hvy, symidx_hvy)
			if hvy[0] in idx_ch2 or hvy[0] in idx_ch3:
				in2_idx[hvy[1]] = hvy[0]+1
		return in1_idx,in2_idx, in1_idx_dipole, in2_idx_dipole

	elif ptype=='perm':
		# define the equivalent index in the 1st file 
		#		for hydrogen atoms not involved in -CH2- -CH3
		for h in symidx_hxx:
			in1_idx[h[0]] = 0
			in1_idx[h[1]] = h[0]+1
			# treating the whole function group to change together
			#	such as the NH2, which not only put the equivalent label to H
			#	but also let charge of atom N to chage freely.
			for k,v in conn.items():
				if h[0] in v and h[1] in v:
					in1_idx[k] =0
		## define the equivalent index in the 2nd file ; remain the label to be -1 to be fixed
		##		for only hydrogen atoms
		#for h in symidx_hxx:
		#	in2_idx[h[0]] = 0
		#	in2_idx[h[1]] = h[0]+1
		#	# treating the whole function group to change together
		#	#	such as the NH2, which not only put the equivalent label to H
		#	#	but also let charge of atom N to chage freely.
		#	for k,v in conn.items():
		#		if h[0] in v and h[1] in v:
		#			in2_idx[k] =0
		if len(idx_ch2) !=0:
			# let charge of atom C to change freely
			for c_ch2 in idx_ch2:
				in2_idx[c_ch2] = 0
			for h_ch2 in symidx_ch2:
			# let hydrogen of -CH2- to maintance the equivalent
				in2_idx[h_ch2[0]] = 0
				in2_idx[h_ch2[1]] = h_ch2[0]+1
		if len(idx_ch3) !=0:
			# let charge of atom C to change freely
			for c_ch3 in idx_ch3:
				in2_idx[c_ch3] = 0
			# let hydrogen of -CH3 to maintance the equivalent
			for h_ch3 in symidx_ch3:
				in2_idx[h_ch3[0]] = 0
				in2_idx[h_ch3[1]] = h_ch3[0]+1
	
		# define the equivalent index in the 1st file
		#		for heavy atoms
		for hvy in symidx_hvy:
			in1_idx[hvy[1]] = hvy[0]+1
		# set the equivalent index in the 2nd file
		#			for heavy atoms
		for hvy in symidx_hvy:
			if hvy[0] in idx_ch2 or hvy[0] in idx_ch3:
				in2_idx[hvy[1]] = hvy[0]+1
		
		# construct bond info, the info. is arranged in the following order:
		#     bond index,    center atom idx.,   bonded atom idx
		bond_info = {}
		bond_idx = 0
		for k,v in conn.items():
			for x in v:
				bond_info[bond_idx] = [k,x]
				bond_idx +=1
		#print('## bond_info')
		#print(bond_info)

		# find bond equivalent info. with the help of symidx_hvy, symidx_ch2, symidx_ch3, and symidx_hxx
		#   the assumption is that, if any two is equivalent, then bond formed between them should be equivalent
		# split bond_equl to the following categories:
		#  -CH2- : 	x C->H  -CH3: o C->H  oth: x HVY-HVY
		#			X H->C		  o H->C       x HVY-HXX
		#						               x HXX-HVY
		bond_equl_ch2 = []
		bond_equl_h2c = []
		bond_equl_ch3 = []
		bond_equl_h3c = []
		bond_equl_oth = []
		for k1,v1 in bond_info.items():
			for k2,v2 in bond_info.items():
				if k1 <k2:
					# find bond equivalent for fragment such as -ch2-, -ch3, -nh2, and -nh3+
					if pair_cmp_list(v1,v2,symidx_hvy,symidx_hvy) or \
						pair_cmp_list(v1,v2,symidx_hxx,symidx_hxx) or \
						pair_cmp_list(v1,v2,symidx_hvy,symidx_hxx) or \
						pair_cmp_list(v1,v2,symidx_hxx,symidx_hvy):
##					if (v1[0] in flatten_list(symidx_hvy) and v2[0] in flatten_list(symidx_hvy) and v1[1] == v2[1]) or \
##						(v1[1] in flatten_list(symidx_hvy) and v2[1] in flatten_list(symidx_hvy) and v1[0] == v2[0]) or\
##						(v1[0] in flatten_list(symidx_ch2) and v2[0] in flatten_list(symidx_ch2) and v1[1] == v2[1]) or\
##						(v1[1] in flatten_list(symidx_ch2) and v2[1] in flatten_list(symidx_ch2) and v1[0] == v2[0]) or\
##						(v1[0] in flatten_list(symidx_ch3) and v2[0] in flatten_list(symidx_ch3) and v1[1] == v2[1]) or\
##						(v1[1] in flatten_list(symidx_ch3) and v2[1] in flatten_list(symidx_ch3) and v1[0] == v2[0]) or\
##						(v1[0] in flatten_list(symidx_hxx) and v2[0] in flatten_list(symidx_hxx) and v1[1] == v2[1]) or\
##						(v1[1] in flatten_list(symidx_hxx) and v2[1] in flatten_list(symidx_hxx) and v1[0] == v2[0]):
##						bond_equl.append([k1,k2])


					# find bond equivalent for -c_ah2,-c_bh2,  where c_a and c_b and h are equivalent
##					if ((v1[0] in flatten_list(symidx_hvy) and v2[0] in flatten_list(symidx_hvy)) \
##						or (v1[0] in flatten_list(symidx_ch2) and v2[0] in flatten_list(symidx_ch2)) \
##						or (v1[0] in flatten_list(symidx_ch3) and v2[0] in flatten_list(symidx_ch3)) \
##						or (v1[0] in flatten_list(symidx_hxx) and v2[0] in flatten_list(symidx_hxx))) \
##						and ((v1[1] in flatten_list(symidx_hvy) and v2[1] in flatten_list(symidx_hvy)) \
##						or (v1[1] in flatten_list(symidx_ch2) and v2[1] in flatten_list(symidx_ch2)) \
##						or (v1[1] in flatten_list(symidx_ch3) and v2[1] in flatten_list(symidx_ch3)) \
##						or (v1[1] in flatten_list(symidx_hxx) and v2[1] in flatten_list(symidx_hxx))):

				#	if ((v1[0] in flatten_list(symidx_hvy) or v2[0] in flatten_list(symidx_hvy)) \
				#		and (v1[0] in flatten_list(symidx_hvy) or v2[0] in flatten_list(symidx_hvy))) \
				#		or ((v1[0] in flatten_list(symidx_ch2) or v2[0] in flatten_list(symidx_ch2)) \
				#		and (v1[0] in flatten_list(symidx_ch2) or v2[0] in flatten_list(symidx_ch2))) \
				#		or ((v1[0] in flatten_list(symidx_ch3) or v2[0] in flatten_list(symidx_ch3)) \
				#		and (v1[0] in flatten_list(symidx_ch3) or v2[0] in flatten_list(symidx_ch3))) \
				#		or ((v1[0] in flatten_list(symidx_hxx) or v2[0] in flatten_list(symidx_hxx)) \
				#		and (v1[0] in flatten_list(symidx_hxx) or v2[0] in flatten_list(symidx_hxx))):
						
						#if k1 not in flatten_list(bond_equl) and k2 not in flatten_list(bond_equl):
						bond_equl_oth.append([k1,k2])
					elif pair_cmp_list(v1,v2,symidx_hvy,symidx_ch2):
						bond_equl_ch2.append([k1,k2])
					elif pair_cmp_list(v1,v2,symidx_ch2,symidx_hvy):
						bond_equl_h2c.append([k1,k2])
					elif pair_cmp_list(v1,v2,symidx_hvy,symidx_ch3):
						bond_equl_ch3.append([k1,k2])
					elif pair_cmp_list(v1,v2,symidx_ch3,symidx_hvy):
						bond_equl_h3c.append([k1,k2])
					#elif 
		#print('## bond_equl')
		#print(bond_equl)
		b_sym_ch2 = {}
		b_sym_h2c = {}
		b_sym_ch3 = {}
		b_sym_h3c = {}
		b_sym_oth = {}
		for b_eq in bond_equl_ch2:
			b_sym_ch2[b_eq[0]]=[]
		for b_eq in bond_equl_ch2:
			b_sym_ch2[b_eq[0]].append(b_eq[1])

		for b_eq in bond_equl_h2c:
			b_sym_h2c[b_eq[0]]=[]
		for b_eq in bond_equl_h2c:
			b_sym_h2c[b_eq[0]].append(b_eq[1])

		for b_eq in bond_equl_ch3:
			b_sym_ch3[b_eq[0]]=[]
		for b_eq in bond_equl_ch3:
			b_sym_ch3[b_eq[0]].append(b_eq[1])

		for b_eq in bond_equl_h3c:
			b_sym_h3c[b_eq[0]]=[]
		for b_eq in bond_equl_h3c:
			b_sym_h3c[b_eq[0]].append(b_eq[1])

		for b_eq in bond_equl_oth:
			b_sym_oth[b_eq[0]]=[]
		for b_eq in bond_equl_oth:
			b_sym_oth[b_eq[0]].append(b_eq[1])

		#print('BOND SYM',b_sym)
		# some equivalent info is redundant, so trim the info. 
		b_sym_trim_ch2 = {}
		for k,v in b_sym_ch2.items():
			if k not in flatten_list(b_sym_ch2.values()):
				b_sym_trim_ch2[k]=v

		b_sym_trim_h2c = {}
		for k,v in b_sym_h2c.items():
			if k not in flatten_list(b_sym_h2c.values()):
				b_sym_trim_h2c[k]=v

		b_sym_trim_ch3 = {}
		for k,v in b_sym_ch3.items():
			if k not in flatten_list(b_sym_ch3.values()):
				b_sym_trim_ch3[k]=v

		b_sym_trim_h3c = {}
		for k,v in b_sym_h3c.items():
			if k not in flatten_list(b_sym_h3c.values()):
				b_sym_trim_h3c[k]=v

		b_sym_trim_oth = {}
		for k,v in b_sym_oth.items():
			if k not in flatten_list(b_sym_oth.values()):
				b_sym_trim_oth[k]=v
		#print('BOND SYM TRIM',b_sym_trim_ch2, b_sym_trim_ch3, b_sym_trim_oth)
		#idx_dip = 0
		for i,in1_d in enumerate(in1_idx_dipole):
			for k,v in b_sym_trim_ch2.items():
				if i==k:
					for x in v:
						in1_idx_dipole[x] = k+1
			for k,v in b_sym_trim_ch3.items():
				if i==k:
					for x in v:
						in1_idx_dipole[x] = k+1
			for k,v in b_sym_trim_oth.items():
				if i==k:
					for x in v:
						in1_idx_dipole[x] = k+1
			#****************************Cautions******************************#
			# following two loops just for the test the accuracy of strategy II
			# where in 1st and 2nd strategy, equivalence of C-H and H-C of -CH3
			# and -CH2- are both maintained
			# when test is done, this command should be commented out!!!!
			#-Test show that Accuracy of Strategy II is better than Strategy I.
			#-In following version, the default strategy is II.
			#******************************************************************#
			if strategy==2:
				#print('*****THIS VERSION IS A TEST VERSION*****')
				for k,v in b_sym_trim_h2c.items():
					if i==k:
						for x in v:
							in1_idx_dipole[x] = k+1
				for k,v in b_sym_trim_h3c.items():
					if i==k:
						for x in v:
							in1_idx_dipole[x] = k+1

		for i,in2_d in enumerate(in2_idx_dipole):
			for k,v in b_sym_trim_ch2.items():
				if i==k:
					in2_idx_dipole[i] = 0
					for x in v:
						in2_idx_dipole[x] = k+1
			for k,v in b_sym_trim_h2c.items():
				if i==k:
					in2_idx_dipole[i] = 0
					for x in v:
						in2_idx_dipole[x] = k+1
			for k,v in b_sym_trim_ch3.items():
				if i==k:
					in2_idx_dipole[i] = 0
					for x in v:
						in2_idx_dipole[x] = k+1
			for k,v in b_sym_trim_h3c.items():
				if i==k:
					in2_idx_dipole[i] = 0
					for x in v:
						in2_idx_dipole[x] = k+1
				
		in1_idx_dipole_reshape = dipole_idx_reshape(in1_idx_dipole,conn)		
		in2_idx_dipole_reshape = dipole_idx_reshape(in2_idx_dipole,conn)		

		return in1_idx, in2_idx, in1_idx_dipole_reshape, in2_idx_dipole_reshape
	else:
		print('WRONG P-TYPE WAS SPECIFIED')
		
#def pair_cmp_list(v1,v2,symidx_a, symidx_b):
#	if v1[0] != v2[0]:
#		for a_list in symidx_a:
#			if (v1[0] in a_list and v2[0] in a_list):
#				for b_list in symidx_b:
#					if (v1[1] in b_list and v2[1] in b_list) or (v1[1] == v2[1]):
#						return True
#	else:
#		for b_list in symidx_b:
#			if (v1[1] in b_list and v2[1] in b_list) or (v1[1] == v2[1]):
#				return True
#	return False

def pair_cmp_list(v1,v2,symidx_a, symidx_b):
	if v1[0] != v2[0]:
		if len(symidx_a) !=0:
			#for a_list in merge_lst(symidx_a):
			for a_list in merge_list(symidx_a):
				if (v1[0] in a_list and v2[0] in a_list):
					if v1[1] == v2[1]:
						return True
					else:
						if len(symidx_b) !=0:
							#for b_list in merge_lst(symidx_b):
							for b_list in merge_list(symidx_b):
								if (v1[1] in b_list and v2[1] in b_list):
									return True
	else:
		if len(symidx_b) !=0:
			#for b_list in merge_lst(symidx_b):
			for b_list in merge_list(symidx_b):
				if (v1[1] in b_list and v2[1] in b_list):
					return True
		
	return False
			
def found_idx(lst,k):
	'''
	lst: [[1,2,3],[2,1],[3,4],[5],[6,7]]
	k = 3
	return: [0,2], [1,2,3,4]
	'''
	idx = []
	for i, js in enumerate(lst):
		for j in js:
			if j == k:
				idx.append(i)
				continue
	mg_lst = []
	for x in idx:
		mg_lst += lst[x]
	return idx, list(set(mg_lst))

def merge_lst(lst):
	'''
	lst: [[1,2,3],[2,1],[3,4],[5],[6,7]]
	return: [[1,2,3,4],[5],[6,7]]
	'''
	outlst = []
	for i,js in enumerate(lst):
		for j in js:
			idx, mg_lst = found_idx(lst,j)
			if mg_lst not in outlst:
				outlst.append(list(mg_lst))
	if outlst == lst:
		return outlst
	else:
		return merge_lst(outlst)


def merge_list(L):
	# convert list to set for utilizing API union
	L = [set(i) for i in L]
	lenth = len(L)
	for i in range(1,lenth):
		for j in range(i):
			if L[i] == {0} or L[j] == {0}:
				continue
			x = L[i].union(L[j])
			y = len(L[i]) + len(L[j])
			if len(x) < y:
				L[i] = x
				L[j] = {0}

	return [list(i) for i in L if i != {0}]



def dipole_idx_reshape(in1_idx_dipole,conn):
	output = {}
	idx = 0
	for k,vs in conn.items():
		output[k] = []
	for k,vs in conn.items():
		for value in vs:
			output[k].append(in1_idx_dipole[idx])
			idx +=1

	return output

def dumpXYZ(filename,data,unit='angstrom'):
	'''
		dump the xyz coordinate of ESP data to xyz file format for visualization.
		- input: 
			+ filename of xyz
			+ data return by getXYZ
				
	'''
	with open(filename,'w') as f:
		f.write("{:10d}\n".format(len(data['ida'].values)))
		f.write('\n')
		for i in range(len(data['ida'].values)):
			if unit == 'au':
				f.write("{:4d} {:10.4f} {:10.4f} {:10.4f}\n".format(data['ida'][i],data['x'][i]/ang2bohr,data['y'][i]/ang2bohr,data['z'][i]/ang2bohr))
			elif unit == 'angstrom':
				f.write("{:4d} {:10.4f} {:10.4f} {:10.4f}\n".format(data['ida'][i],data['x'][i],data['y'][i],data['z'][i]))




if __name__ == '__main__':
	#filename = 'example/C2H6_b3lyp_321g_esp.dat'
	#filename = '../C5NH5_b3lyp_321g_esp.dat'
	#filename = '../PO3CH5_b3lyp_321g_esp.dat'
	#filename = '../C4H10_b3lyp_321g_esp.dat'
	#filename = 'example/CH3S2CH3_b3lyp_321g_esp.dat'
	#filename = 'example/CH3SO2CH3_mp2_a4z_esp.dat'
	#filename = 'example/C6H6_b3lyp_321g_esp.dat'
	#filename = 'example/CH3NH2_mp2_a4z_esp.dat'
	#filename = 'example/CHNHOH_mp2_a4z_esp.dat'
	#filename = 'example/qiang_test_case_esp_b321.dat'
	#filename = 'example/esp.dat'
	#filename = 'example/CH3F_ccsd_a4z_esp.dat'
	#filename = 'example/bis_esp.dat'
	#filename = '../../test/CH2O_ccsd_a4z_esp.dat'
	filename = 'example/C3H8_mp2_a4z_esp.dat'
	xyz = getXYZ(filename)
	print(xyz[['x','y','z']])
	dMat = distMatrix(xyz)
	print('# DMatrix')
	print(dMat)
	print(dMat.columns)
	print(dMat.index)
	Nbond = getBOND(xyz)
	print('Num. of bonds is: {}'.format(Nbond))
	conn,conAtm,conIdx = getCONN(xyz)
	#conn,conAtm = getCONN(xyz)
	print('#Variable $conn$ content:')
	print(conn)
	print('#Variable $conAtm$ content:')
	print(conAtm)
	print('#Variable $conIdx$ content:')
	print(conIdx)
	equAtmMap = getEQU(conn,xyz,depth=3)
	idx_ch2,idx_ch3 = getCH2CH3(conIdx,xyz)
	print(equAtmMap)
	print(idx_ch2, idx_ch3)
	print(type(equAtmMap.values))
	sym_hvy, sym_ch2, sym_ch3, sym_hxx = getSymIdx(equAtmMap, xyz, idx_ch2, idx_ch3,conn)
	print('HVY atom \t H(-CH2-) \t H(-CH3) \t Other H')
	print(sym_hvy,sym_ch2,sym_ch3, sym_hxx )
	in1_idx,in2_idx,in1_idx_dipole, in2_idx_dipole = dumpSymIdx(sym_hvy, sym_ch2, sym_ch3, sym_hxx,idx_ch2, idx_ch3,xyz,conn,ptype='perm')
	print(in1_idx)
	print(in2_idx)
	print(in1_idx_dipole)
	print(in2_idx_dipole)
	#print(equAtmMap.shape[0],equAtmMap.shape[1])
	dumpXYZ('test.xyz',xyz,unit='angstrom')
