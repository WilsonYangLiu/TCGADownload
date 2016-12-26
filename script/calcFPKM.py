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

def FPKM(Table, GeneLenDict):
	'''
	Parameters:
		Table: DataFrame. It stores the read count table for each gene in each sample
		GeneLenDict: a Dict stores the length of each gene
		
	Returns:
		Table: DataFrame. It stores the FPKM table for each gene in each sample
	'''
	Lib = np.sum(Table, axis=0)
	for col in Table.columns:
		for gene in GeneLenDict.keys():
			Table.ix[gene, col] = 10e9 * Table.ix[gene, col] / (GeneLenDict[gene] * Lib[col] )
			
	return Table
	
if __name__ == '__main__':
	os.chdir(r'E:/Project_G/db.TCGA/TCGADownloader/script')
	
	Table = pd.read_csv(r'../data/COADREAD_trans_GeneExp_Counts.csv', index_col=0)
	
	with open(r'../data/gencode.v22.annotation.used4FPKM.exon.csv', 'rb') as csvfile:
		spamreader = csv.reader(csvfile)
		GeneLenDict = {line[0]:int(line[2]) for line in islice(spamreader, 1, None) }
		
	Table = FPKM(Table, GeneLenDict)
	Table = Table.iloc[:-5,:]
	Table.to_csv(r'../data/COADREAD_trans_GeneExp_Counts2FPKM.exon.csv')
	
	
	
	
	
	
	