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
from scipy.stats import pearsonr
from scipy.stats import norm
from time import ctime

class SSNClient(object):
	def __init__(self, reference, experiment, background, silent=True):
		'''
		Parameters:
			reference: DataFrame (N * m). the expression dataframe of reference network
			experiment: Series (N * 1). the expression of one single sample. 
				The index was the same as reference.
			background: List (N * 2), each item (row) is a paired links (such as two gene). the background network. 
				Usually recheived from prebuilt interactive network, such as KEGG, STRING and so on.
			silent: Boolean. default to True, means the result will write to file directly. 
				Otherwise, the calcutale result will save to memory, which was not recommended. 
				However, when True, you can still use the meathod 'onePair' to calcutale single-sample network of a pair one by one
		'''
		if not np.all(reference.index == experiment.index):
			print(r'[ERROR] Index between reference and experiment must be same')
			raise IndexError
			
		self.indexes = reference.index	# genes that contains expression data
			
		if not experiment.name:
			experiment.name = 'expt'
		
		self.reference = reference
		self.experiment = experiment
		self.background = background
		self.silent = silent
		
		self.num = self.reference.shape[1]	# the number of reference samples
		if self.num < 3:
			print(r'[ERROR] reference must contains more than or equel 3 samples')
			raise Exception
			
		if not silent:
			self.Result = {}
			
		self.nodeNotInIndex = set()		# the nodes dose not presented in self.indexes
	
	def _calc_pval(self, num, delta, pcc):
		'''
		Parameters:
			num: the number of sample in reference
			delta: PCC_n+1 - PCC_n
			pcc: PCC_n
			
		Returns:
			p_value: the p value of the selected link
				If pcc is 1, which means the two nodes are exactly the same, the p value would be set 1.0
		'''
		if pcc ** 2 != 1:
			z_scores = (num - 1) * delta / (1 - pcc ** 2 )
			return norm.sf(np.abs(z_scores) ) * 2
		else:
			return 1.0
		
	def onePair(self, pr):
		'''
		Parameters: 
			pr: a 2-element list or tuple. the selected gene pair
			
		Returns:
			p_value or -1: -1 means cannot get the p_value due to some Error exist. See the printed messages
		'''
		if not isinstance(pr, list):
			pr = list(pr)
		
		if (pr[0] in self.indexes) and (pr[1] in self.indexes):
			# calcutale PCC_n
			gene1 = self.reference.ix[pr[0] ]; gene2 = self.reference.ix[pr[1] ]
			pcc1 = pearsonr(gene1, gene2)[0]
			# calcutale PCC_n+1
			gene1['expt'] = self.experiment[pr[0] ]; gene2['expt'] = self.experiment[pr[1] ]
			pcc2 = pearsonr(gene1, gene2)[0]
			# Get the p value of the selected link
			return pcc1, pcc2, self._calc_pval(self.num, pcc2-pcc1, pcc1)

		else:
			with open(r'{}.log'.format(self.experiment.name), 'a+') as handle:
				if pr[0] not in self.indexes:
					self.nodeNotInIndex.add(pr[0] )
					handle.write('[WARNING] node [{}] in indexes\n'.format(pr[0] ) )
				if pr[1] not in self.indexes:
					self.nodeNotInIndex.add(pr[1] )
					handle.write('[WARNING] node [{}] in indexes\n'.format(pr[1] ) )
			
			return -1,-1,-1
			
	def calcAllPair(self):
		'''
		
		'''
		if self.silent:
			with open(r'{}.{}.tsv'.format(self.experiment.name, '_'.join(ctime().replace(':', '_').split(' ')[:3] ) ), 'wb') as handle:
				handle.write('#node1\tnode2\tPCC_n\tPCC_n+1\tp_val\n')
				for pr in self.background:
					pcc1, pcc2, p_value = self.onePair(pr)
					if p_value != -1:
						handle.write('{}\t{}\t{}\t{}\t{}\n'.format(pr[0], pr[1], pcc1, pcc2, p_value) )
		
		else:
			for pr in self.background:
				if not self.Result.has_key(tuple(sorted(pr) ) ):
					_, _, p_value = self.onePair(pr)
					if p_value != -1:
						self.Result[tuple(sorted(pr) ) ] = p_value
				
	
def remove_vesion(idxes, sep='.'):
	'''
	Parameters:
		idxes: Index object. the index in idxes contains the suffix as version number, sepetated by 'sep', default to '.'. 
		sep: the delimeter between id and version number. Default to '.'
		
	Returns:
		newIdx: Index object with the order of idxes. the item in newIdx does not contain version suffix
	'''
	newIdx = pd.Index([item.split(sep)[0] for item in idxes], dtype=object)
	return newIdx

def work(filename):
	#os.chdir(r'E:/Project_G/db.TCGA/SSN')
	os.chdir(r'/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_G/db.TCGA/SSN')
	
	# Get the background network
	with open(r'link.higher9.ENSG.tsv', 'rb') as handle:
		links = [line.strip().split('\t') for line in handle]
	
	# Get the expression dataframe of reference network. This dataframe was obtained by performing the normalization presedure
	refExprDF = pd.read_csv(r'COADREAD_trans_GeneExp_TPM_Normal_selected.csv', index_col=0)
	refExprDF.index = remove_vesion(refExprDF.index)
	
	# Get the expression dataframe of all testable experiment sample.
	exptExprDF = pd.read_csv(filename, index_col=0)
	exptExprDF.index = remove_vesion(exptExprDF.index)
	exptExprDF = exptExprDF.reindex(refExprDF.index)
	
	for i in range(exptExprDF.shape[1] ):
		ssn = SSNClient(refExprDF, exptExprDF.iloc[:,i], links, silent=True)
		print('[WORKING] {}...'.format(exptExprDF.iloc[:,i].name) )
		ssn.calcAllPair()
		print('[DONE] {}.\n   {} of nodes does not exist in index'.format(exptExprDF.iloc[:,i].name, len(ssn.nodeNotInIndex) ) )
		
if __name__ == '__main__':
	work(r'COADREAD_trans_GeneExp_TPM_MSI-H.csv')
