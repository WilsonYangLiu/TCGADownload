from __future__ import print_function

#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Copyright (c) 2016 Wei-Xin Liu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

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
