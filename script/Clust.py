from __future__ import print_function

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

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

import os
import numpy as np
import pandas as pd
from pandas import DataFrame, Series

def ssnDistMat(fileList):
	'''
	Calculate the dissimilarity matrix (n * n) base on single-sample network
	
	Parameters:
		fileList: List. a list contained the files of each filtered SSN.
		
	Returns:
		Matrix: ndarray. the dissimilarity matrix
	'''
	num = len(fileList)
	Matrix = np.zeros(shape=(num, num), dtype=np.int)
	for i in range(num):
		print(r'[WORKING] {} ...'.format(fileList[i] ) )
		for j in range(i+1, num):
			s1 = set(pd.read_csv(fileList[i], index_col=0, sep='\t').index)
			s2 = set(pd.read_csv(fileList[j], index_col=0, sep='\t').index)
			Matrix[i, j] = Matrix[j, i] = len(s1 & s2)

		print(r'[DONE]')
	
	Matrix = np.max(Matrix) - Matrix
	for i in range(num):
		Matrix[i, i] = 0
		
	return Matrix

if __name__ == '__main__':
	#os.chdir(r'E:/Project_G/db.TCGA/SSN/Clust')
	os.chdir(r'/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_G/db.TCGA/SSN/Clust')
	
	MSI_L = [''.join(['../filteredSSN_MSI-L/', f] ) for f in os.listdir('../filteredSSN_MSI-L') if f.endswith('filtered.tsv') ]
	MSI_H = [''.join(['../filteredSSN_MSI-H/', f] ) for f in os.listdir('../filteredSSN_MSI-H') if f.endswith('filtered.tsv') ]
	allFile = MSI_L + MSI_H
	Matrix = ssnDistMat(allFile)
