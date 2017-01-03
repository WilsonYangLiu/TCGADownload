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
import pandas as pd
import numpy as np
from pandas import DataFrame, Series

def ssnFilter(ssn, threshold=0.05):
	'''
	Parameters: 
		ssn: a DataFrame containing all the edges extracted from single-sample network meathod.
			each row contains [node1, node2, PCC_n, PCC_n+1, p_val]
		threshold: num. the cutoff of the significant edges
	
	Returns:
		a DataFrame, its index was the edge, for example, 'ENSG00000262814-ENSG00000265681'. 
			each contains [PCC_n, PCC_n+1, p_val]
	'''
	filtered = ssn[ssn['p_val'] < threshold]
	
	ssnDict = {}
	for row in range(filtered.shape[0] ):
		oneRow = filtered.iloc[row, :]
		rowKey = '-'.join(sorted([oneRow[0], oneRow[1] ] ) )
		if ssnDict.has_key(rowKey):
			assert ssnDict[rowKey][2] == oneRow[4]
		else:
			ssnDict[rowKey] = oneRow[2:]
			
	return DataFrame(ssnDict).T

if __name__ == '__main__':
	#os.chdir(r'E:/Project_G/db.TCGA/SSN/rawSSN')
	os.chdir(r'/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_G/db.TCGA/SSN/rawSSN_MSI-L')
	
	ssnFile = [f for f in os.listdir('./') if f.endswith('Tue_Jan_.tsv') ]

	for f in ssnFile:
		print(r'[WORKING] {} ...'.format(f) )
		ssn = pd.read_csv(f, sep='\t')
		filtered = ssnFilter(ssn)
		filtered.to_csv('.'.join(f.split('.')[:-1]+['filtered.tsv'] ), sep='\t')
		print(r'[DONE] save to ', '.'.join(f.split('.')[:-1]+['filtered.tsv'] ) )
	
	