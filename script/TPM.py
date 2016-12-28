from __future__ import print_function, division

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

import os, csv
import pandas as pd
import numpy as np
from pandas import DataFrame, Series
from itertools import islice
from FPKM import calcFPKM

def calcTPM(Table, GeneLenDict, from_RAW=True):
	'''
	Parameters:
		Table: DataFrame. It stores the read count table for each gene in each sample
		GeneLenDict: a Dict stores the length of each gene
		
	Returns:
		tmpTable: DataFrame. It stores the TPM table for each gene in each sample
	'''
	if from_RAW:
		Table = calcFPKM(Table, GeneLenDict)
		if isinstance(Table, DataFrame):
			Table = Table.iloc[:-5,:]
		elif isinstance(Table, Series):
			Table = Table[:-5]
		else:
			print(r'the type of data should be DataFrame or Series')
			raise Exception

	if isinstance(Table, DataFrame):
		tpmTable = DataFrame(Table)
		Lib = np.sum(Table, axis=0)
		for col in Table.columns:
			for gene in GeneLenDict.keys():
				tpmTable.ix[gene, col] = 10e6 * Table.ix[gene, col] / Lib[col]
				
	elif isinstance(Table, Series):
		tpmTable = Series(Table)
		Lib = np.sum(Table, axis=0)
		for gene in GeneLenDict.keys():
			tpmTable[gene] = 10e6 * Table[gene] / Lib
	
	else:
		print(r'the type of data should be DataFrame or Series')
		raise Exception
			
	return tpmTable
	
if __name__ == '__main__':
	os.chdir(r'E:/Project_G/db.TCGA/TCGADownloader/script')
	
	Table = pd.read_csv(r'../data/COADREAD_trans_GeneExp_Counts.csv', index_col=0)
	
	with open(r'../data/gencode.v22.annotation.used4FPKM.csv', 'rb') as csvfile:
		spamreader = csv.reader(csvfile)
		GeneLenDict = {line[0]:int(line[2]) for line in islice(spamreader, 1, None) }
		
	Table = calcTPM(Table.iloc[:,0], GeneLenDict)
	#Table.to_csv(r'../data/COADREAD_trans_GeneExp_Counts2TPM.csv')
	
	