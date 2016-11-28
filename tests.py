from __future__ import print_function

"""
    Unit tests for TCGADownloader.
"""

import os
from TCGADownloader.gdcFileClient import *

os.chdir(r'E:/Project_G/db.TCGA/TCGADownloader')
#os.chdir(r'/media/wilson/b776f228-366c-4e52-acd6-65df5b458e8c/Project_G/db.TCGA/TCGADownloader')

def example_1():
	print('Example 1:')
	Downloader = gdcFileClient(r'./data/tmp_manifest.trans_GeneExp.tsv', DIR=r'./data',\
								uuid_info_file=r'uuid_info_file.trans_GeneExp.tsv')
	#Downloader.fileInfo()
	Downloader.fileDownload()
	
def example_2():
	print('Example 2:')
	Downloader = gdcSNVClient(r'./data/tmp_manifest.snv_mSomatic.tsv', DIR=r'./data',\
								uuid_info_file=r'uuid_info_file.snv_mSomatic.json')
	#Downloader.fileInfo()
	Downloader.fileDownload()

if __name__ == '__main__':
	example_1()
	example_2()