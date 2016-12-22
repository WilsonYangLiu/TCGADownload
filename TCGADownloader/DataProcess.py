from __future__ import print_function

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import os, gzip
import numpy as np
import pandas as pd
from pandas import Series, DataFrame

def filter(obj, filt=None):
	'''
	Parameters:
		obj: a DataFrame object. row are observation, columns are variables.
		filt: a self-definded operation. format: {variable: (op, value)}
			variable is one of the variables in obj. value is the posible value of variable
			op is the operation between variable and value, currently support '==', '!=', 'in', 'not in'
			
	Returns:
		a nameList that contains the list of filenames. each file correspond to one sample
	'''
	for key in filt.keys():
		if filt[key][0] == '==':
			obj = obj[obj[key] == filt[key][1] ]
			
		elif filt[key][0] == '!=':
			obj = obj[obj[key] != filt[key][1] ]
			
		elif filt[key][0] == 'in':
			col = obj[key].values
			row = []
			for item in col:
				if item in filt[key][1]:
					row.append(True)
				else:
					row.append(False)
					
			obj = obj[row]
			
		elif filt[key][0] == 'not in':
			col = obj[key].values
			row = []
			for item in col:
				if item in filt[key][1]:
					row.append(False)
				else:
					row.append(True)
					
			obj = obj[row]
	
	return obj

def integration(nameList, DIR='./'):
	'''
	intergrate selected samples together.
	
	Parameters:
		nameList: the list of filenames. each file correspond to one sample
		
	Returns:
		allSample: a DataFrame object. the intergrated data
	'''
	allSample = DataFrame()
	for name in nameList:
		aDict = {}
		with gzip.open(''.join([DIR, name] ), 'rb') as handle:
			for line in handle:
				line = line.split('\t')
				aDict[line[0] ] = line[1]
			
		oneSample = Series(aDict, dtype=np.float)
		if allSample is None:
			oneSample.name = name.split('.')[0]
			oneSample.index.name = 'Gene'
			allSample = DataFrame(oneSample)
		else:
			allSample[name.split('.')[0] ] = oneSample
		
	return allSample	
	
if __name__ == '__main__':
	pass
