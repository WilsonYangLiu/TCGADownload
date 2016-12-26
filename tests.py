from __future__ import print_function

"""
    Unit tests for TCGADownloader.
"""

import os
import pandas as pd
from pandas import DataFrame, Series
from TCGADownloader.gdcFileClient import *
from TCGADownloader import DataProcess

os.chdir(r'E:/Project_G/db.TCGA/TCGADownloader')
#os.chdir(r'/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_G/db.TCGA/TCGADownloader')

def expression():
	#print('Example 1:')
	Downloader = gdcFileClient(r'./data/tmp_manifest.trans_GeneExp.tsv', DIR=r'./data',\
								uuid_info_file=r'uuid_info_file.trans_GeneExp.tsv')
	#Downloader.fileInfo()
	Downloader.fileDownload()
	
def CNV():
	Downloader = gdcFileClient(r'./data/tmp_manifest.CNV_masked.tsv', DIR=r'./data',\
								uuid_info_file=r'uuid_info_file.CNV_masked.tsv')
	#Downloader.fileInfo()
	Downloader.fileDownload()
	
def SNV():
	#print('Example 2:')
	Downloader = gdcSNVClient(r'./data/tmp_manifest.snv_mSomatic.tsv', DIR=r'./data',\
								uuid_info_file=r'uuid_info_file.snv_mSomatic.json')
	Downloader.fileInfo()
	Downloader.fileDownload()

def clinical():
	Downloader = gdcFileClient('./data/tmp_manifest.clinical.tsv', DIR=r'./data',\
								uuid_info_file=r'uuid_info_file.clinical.tsv')
	Downloader._fields_retieved_str = ['file_id', 'cases_0_submitter_id', 'file_name']
	Downloader.set_fields_retieved()
	Downloader.fileDownload()
	
def dataProcess():
	# Load the info file first
	uuid_info_DF = pd.read_table(r'./data/uuid_info_file.trans_GeneExp.tsv', sep='\t')
	# Filter the info by using the `DataProcess.filter` meathod
	obj = DataProcess.filter(uuid_info_DF, filt={'cases_0_samples_0_sample_type': ('==', 'Primary Tumor'),\
										'analysis_workflow_type': ('==', 'HTSeq - Counts') } )
	obj = obj.ix[:, ['file_name', 'cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id'] ]
	# Combine the info of obj to construct the file name defined in `gdcFileClient`: [barcode].[origin filename]
	nameList = []
	for idx in range(obj.shape[0]):
		nameList.append('.'.join([obj.iloc[idx][1], obj.iloc[idx][0] ] ) )
	# Integrate the selected file
	itData = DataProcess.integration(nameList, DIR='./data/')
	itData.to_csv(r'./data/COADREAD_trans_GeneExp_Counts.csv')

if __name__ == '__main__':
	#expression()
	#SNV()
	#clinical()
	#dataProcess()
	pass
