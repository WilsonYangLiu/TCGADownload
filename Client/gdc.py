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

import os
import requests
import json
from itertools import islice
from time import sleep

__all__ = ["BaseClient", "gdcSNVClient", "gdcFileClient"]

class DownloadError(KeyError):
	pass
	
class BaseClient(object):
	def __init__(self, manifest_file, DIR='.'):
		'''
		Parameters:
			manifest_file: file containing the file id (uuid) from the TCGA portal is the one needed in my class
			DIR: the directory used for save the data. Default to currently directory -- '.'
		'''
		self.manifest_file = manifest_file
		self.DIR = DIR
		
		self.uuid = []
		self.uuid2filename = {}
		
	def set_manifest(self, manifest_file):
		self.manifest_file = manifest_file
				
	def clean_uuid(self):
		#self.clean_uuid2filename()
		self.uuid = []

	def clean_uuid2filename(self):
		#self.clean_uuid()
		self.uuid2filename = {}
	
	def get_uuid(self):
		'''	
		Get file id (uuid) from manifest_file
		'''
		self.clean_uuid()
		with open(self.manifest_file, 'rb') as f:
			for line in islice(f, 1, None):
				self.uuid.append(line.strip().split('\t')[0] )	
				
	def _SplitList(self, uuid, interval):
		'''
		Split the whole set of file id (uuid) into small subsets with the numbers of interval
		
		Parameters:
			uuid: the file id (uuid) set
			interval: the number of file id in the subset
			
		Returns:
			newID: a list containing all subsets
		'''
		LEN = len(uuid)
		BIN = range(0, LEN, interval); BIN.append(LEN)
		LEFT = BIN[0]; RIGHT = BIN[1]
		newID = [uuid[LEFT:RIGHT] ]
		for idx in islice(BIN, 2, None):
			LEFT = RIGHT
			RIGHT = idx
			newID.append(uuid[LEFT:RIGHT] )
			
		return newID
		
	def responseInfo(self, idList=None, format="TSV", query_fields=None):
		'''
		get the response from TCGA server
		
		Parameters:
			idList: file id (uuid) to be sended to TCGA server
			format: the return format from the TCGA server. Default to "TSV", posssible value: "TSV", "json"
			query_fields: the searched fields corresponding to file id (uuid) on the TCGA server
				More fields can be checked on the server: https://gdc-docs.nci.nih.gov/API/Users_Guide/Appendix_A_Available_Fields/
		
		Returns:
			response: the result from the server
		'''
		if not idList:
			self.get_uuid()
			idList = self.uuid
			
		if not query_fields:
			query_fields = "file_id,file_name,cases.submitter_id,cases.case_id,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id,analysis.workflow_type"
		
		file_endpt = 'https://gdc-api.nci.nih.gov/files'
		
		filt = {"op":"in",
				"content":{
					"field": "files.file_id",
					"value": idList
				}
		}
		params = {'filters':json.dumps(filt), 
				  "format":format,
				  "fields":query_fields,
				  "size":"100"
		}

		responseStatus = True
		response = None
		while responseStatus:
			connectionStatus = True
			while connectionStatus:
				try:
					response = requests.post(file_endpt, params=params, headers={'Content-Type': 'application/json'} )
					connectionStatus = False
				except Exception, e:
					print('[fileInfo Conncetion] Could not establish connection: {}'.format(e) )
					print('   Wait 5 seconds, reconnecting...')
					sleep(5)
				
			if response.status_code == requests.codes.ok:
				responseStatus = False
				print('[fileInfo Response, # of uuids downloading: {}] OK: code[{}], continue.'.format(len(idList), response.status_code) )
			else:
				print('[fileInfo Response] Response Error: code[{}], try again... '.format(response.status_code) )
				sleep(5)
				
		return response
	
	def fileInfo(self):
		'''
		Get the infomation corresponding to file id (uuid) from the TCGA server
		'''
		pass
		
	def gdcDownload(self, ID, name):
		'''
		Download the file based on the file id (uuid) from the server
		
		Parameters:
			ID: file id (uuid)
			name: save file name
		'''
		data_endpt = 'https://gdc-api.nci.nih.gov/data/'
		
		responseStatus = True
		while responseStatus:
			connectionStatus = True
			response = None
			while connectionStatus:
				try:
					response = requests.get(data_endpt + ID)
					connectionStatus = False
				except Exception, e:
					print('[gdcDownload Conncetion] Could not establish connection: {}'.format(e) )
					print('   [ID: {}]'.format(ID) )
					print('   Wait 5 seconds, reconnecting...')
					sleep(5)
				
			if response.status_code == requests.codes.ok:
				responseStatus = False
				print('[gdcDownload Response] OK: code[{}]'.format(response.status_code) )
			else:
				print('[gdcDownload Response] Response Error: {}, try again... '.format(response.status_code) )
				print('   [{}]'.format(ID) )
				sleep(5)
				continue
				
			with open('/'.join([self.DIR, name] ), 'wb') as handle:
				for block in response.iter_content(1024):
					handle.write(block)
					
	def fileDownload(self):
		''' 
		Download a set of files base on the uuids
		'''
		pass
	
	def _set_uuid2filename(self):
		'''
		Construct the display file name in disk from the '*_file' for the corresponding file id (uuid)
		'''
		pass
					
class gdcSNVClient(BaseClient):
	def __init__(self, manifest_file, DIR = '.', uuid_info_file=r'uuid_info_file.json', format="json", query_fields=None):
		'''
		Parameters:
			manifest_file: file containing the file id (uuid) from the TCGA portal is the one needed in my class
			DIR: the directory used for save the data. Default to currently directory -- '.'
			uuid_info_file: the file stores the infomation corresponding to file id (uuid). Default to 'uuid_info_file.tsv'
			query_fields: the searched fields corresponding to file id (uuid) on the TCGA server
				More fields can be checked on the server: https://gdc-docs.nci.nih.gov/API/Users_Guide/Appendix_A_Available_Fields/
		'''
		BaseClient.__init__(self, manifest_file, DIR)
		
		self.uuid_info_file = uuid_info_file
		self.format = format
		if query_fields:
			self.query_fields = query_fields
		else:
			self.query_fields = "file_id,file_name,cases.submitter_id,cases.case_id,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id,analysis.workflow_type"
		
	def fileInfo(self):
		'''
		Get the infomation corresponding to file id (uuid) from the TCGA server
		'''
		self.get_uuid()
		
		response = self.responseInfo(idList=self.uuid, format=self.format, query_fields=self.query_fields)
		parsed = json.loads(response.text)
		print('[fileInfo Saving] Save to "{}"'.format('/'.join([self.DIR, self.uuid_info_file] )) )
		with open('/'.join([self.DIR, self.uuid_info_file] ), 'wb') as handle:
			handle.write(json.dumps(parsed, indent=4, sort_keys=False) )
			
	def gdcDownload(self, ID, name):
		'''
		Download the file based on the file id (uuid) from the server
		
		Parameters:
			ID: file id (uuid)
			name: save file name
		'''
		BaseClient.gdcDownload(self, ID, name)
		
	def fileDownload(self, END=None):
		''' 
		Download a set of files base on the uuids
		'''
		self._set_uuid2filename()
		for ID in self.uuid:
			self.gdcDownload(ID, self.uuid2filename[ID] )
			sleep(1)
			
	def _set_uuid2filename(self):
		'''
		Construct the display file name in disk from the 'manifest_file' for the corresponding file id (uuid)
		'''
		self.get_uuid()
		with open(self.manifest_file, 'rb') as f:
			for line in islice(f, 1, None):
				line = line.strip().split('\t')
				line[1] = '.'.join([item for i, item in enumerate(line[1].split('.') ) if i != 3 ] )
				self.uuid2filename[line[0] ] = line[1]
	
class gdcFileClient(BaseClient):
	def __init__(self, manifest_file, DIR='.', uuid_info_file=r'uuid_info_file.tsv', format="TSV", query_fields=None,\
				interval=100, fields_retieved=[7, 8, 3], fields_retieved_str=None):
		'''
		Parameters:
			manifest_file: file containing the file id (uuid) from the TCGA portal is the one needed in my class
			DIR: the directory used for save the data. Default to currently directory -- '.'
			uuid_info_file: the file stores the infomation corresponding to file id (uuid). Default to 'uuid_info_file.tsv'
			query_fields: the searched fields corresponding to file id (uuid) on the TCGA server
				More fields can be checked on the server: https://gdc-docs.nci.nih.gov/API/Users_Guide/Appendix_A_Available_Fields/
			interval: the number of uuids searched at a time. Default to 100
			fields_retieved: the fields that retieved from the uuid_info_file. Number means column number (0-based)
				Note that first element must be the file id (uuid). Default to [7, 3, 8], meaning that the 7th file id (uuid),
				3rd file name on the TCGA server and 8th the aliquots submitter id (cases_samples_portions_analytes_aliquots_submitter_id)
			fields_retieved_str: the STRING verfion of the fields that retieved from the uuid_info_file.
				For clinic data: ['file_id', 'cases_0_submitter_id', 'file_name']
				for trancripome data and CNV data: ['file_id', 'cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id', 'file_name']
				Default to trancripome data and CNV data
		'''
		BaseClient.__init__(self, manifest_file, DIR)
		
		self.uuid_info_file = uuid_info_file
		self.format = format
		if query_fields:
			self.query_fields = query_fields
		else:
			self.query_fields = "file_id,file_name,cases.submitter_id,cases.case_id,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id,analysis.workflow_type"
		self.interval = interval
		self.fields_retieved = fields_retieved
		if fields_retieved_str:
			self._fields_retieved_str = fields_retieved_str
		else:
			self._fields_retieved_str = ['file_id', 'cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id', 'file_name']
		
		self._uuid_info_file_header = None

	def fileInfo(self, END=None):
		'''
		Get the infomation corresponding to file id (uuid) from the TCGA server
		'''
		self.get_uuid()
		
		split_uuid = self._SplitList(self.uuid, self.interval)
		Data = ''
		for ids in islice(split_uuid, 0, END):
			response = self.responseInfo(idList=ids, format=self.format, query_fields=self.query_fields)
					
			# 将 colname 与 subject 分开
			# Be careful for '\r', Cause may be other condition besed on the system
			colname_endidx = response.text.index('\r')
			subject = response.text[colname_endidx:-2]
			if not Data:
				self._uuid_info_file_header = response.text[:colname_endidx]
				Data = response.text[:colname_endidx]
			
			Data = ''.join([Data, subject] )
			
		print('[fileInfo Saving] Save to "{}"'.format('/'.join([self.DIR, self.uuid_info_file] ) ) )
		with open('/'.join([self.DIR, self.uuid_info_file] ), 'wb') as handle:
			handle.write(Data)
	
	def gdcDownload(self, ID, name):
		'''
		Download the file based on the file id (uuid) from the server
		
		Parameters:
			ID: file id (uuid)
			name: save file name
		'''
		BaseClient.gdcDownload(self, ID, name)
	
	def fileDownload(self, END=None, uuidFromFile=True):
		''' 
		Download a set of files base on the uuids
		
		Parameters:
			uuidFromFile: Does the uuid from the uuid_info_file? Default to `True`
				If `False`, it means you can set your own uuid instead from uuid_info_file
		'''
		if not self.uuid:
			print('[Warning] this instance does not contain any file id (uuid)\n\tFile id would be extracted from [{}]'.format(self.manifest_file))
			uuidFromFile = True
		
		self._set_uuid2filename(uuidFromFile=uuidFromFile)
			
		for ID in self.uuid:
			self.gdcDownload(ID, self.uuid2filename[ID] )
			sleep(1)

	def set_fields_retieved(self, fields_retieved_str=None):
		'''
		Determine the indexes of fields_retieved_str in the 'uuid_info_file' header (self._uuid_info_file_header)
		
		Parameters:
			fields_retieved_str: the STRING verfion of the fields that retieved from the uuid_info_file.
		'''
		if not fields_retieved_str:
			fields_retieved_str = self._fields_retieved_str
		
		if not self._uuid_info_file_header:
			if not os.path.isfile('/'.join([self.DIR, self.uuid_info_file] ) ):
				self.fileInfo()
				
			with open('/'.join([self.DIR, self.uuid_info_file] ), 'rb') as f:
				for line in islice(f, 0, 1):
					self._uuid_info_file_header = line.strip()
					
		header = self._uuid_info_file_header.strip().split('\t')
		try:
			self.fields_retieved = [header.index(str) for str in fields_retieved_str]
		except Exception, e:
			print(e)
			print('''
[Suggest] If you try to download "clinic" data, set: 
`your instance name`._fields_retieved_str = ['file_id', 'cases_0_submitter_id', 'file_name'].
Ortherwise, check your uuid_info_file and set your own `your instance name`._fields_retieved_str
''')

	def _set_uuid2filename(self, END=None, uuidFromFile=True):
		'''
		Construct the display file name in disk from the 'uuid_info_file' for the corresponding file id (uuid)
		
		Parameters:
			uuidFromFile: Does the uuid from the uuid_info_file? Default to `True`
				If `False`, it means you must set your own uuid instead from uuid_info_file
		'''
		if uuidFromFile:
			self.get_uuid()
			
		if not os.path.isfile('/'.join([self.DIR, self.uuid_info_file] ) ):
			self.fileInfo()
			
		with open('/'.join([self.DIR, self.uuid_info_file] ), 'rb') as f:
			try:
				for line in islice(f, 1, END):
					line = line.strip().split('\t')
					self.uuid2filename[line[self.fields_retieved[0] ] ] = '.'.join([line[self.fields_retieved[1] ], line[self.fields_retieved[2] ] ] )
			except Exception, e:
				print(e)
				print('''
[Suggest] If you try to download "clinic" data, set: 
`your instance name`._fields_retieved_str = ['file_id', 'cases_0_submitter_id', 'file_name'].
Then run `your instance name`.set_fields_retieved().
Ortherwise, check your uuid_info_file and set your own `your instance name`._fields_retieved_str 
and run `your instance name`.set_fields_retieved().
''')

