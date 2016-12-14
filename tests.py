from __future__ import print_function

"""
    Unit tests for TCGADownloader.
"""

import os
from TCGADownloader.gdcFileClient import *

os.chdir(r'E:/Project_G/db.TCGA/TCGADownloader')
#os.chdir(r'/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_G/db.TCGA/TCGADownloader')

def expression():
	#print('Example 1:')
	Downloader = gdcFileClient(r'./data/tmp_manifest.trans_GeneExp.tsv', DIR=r'./data',\
								uuid_info_file=r'uuid_info_file.trans_GeneExp.tsv')
	#Downloader.fileInfo()
	Downloader.fileDownload()
	
def SNV():
	#print('Example 2:')
	Downloader = gdcSNVClient(r'./data/tmp_manifest.snv_mSomatic.tsv', DIR=r'./data',\
								uuid_info_file=r'uuid_info_file.snv_mSomatic.json')
	#Downloader.fileInfo()
	Downloader.fileDownload()

def clinical():
	Downloader = gdcFileClient('./data/tmp_manifest.clinical.tsv', DIR=r'./data',\
								uuid_info_file=r'uuid_info_file.clinical.tsv')
	Downloader._fields_retieved_str = ['file_id', 'cases_0_submitter_id', 'file_name']
	Downloader.set_fields_retieved()
	Downloader.fileDownload()

if __name__ == '__main__':
	#expression()
	#SNV()
	#clinical()
	pass
