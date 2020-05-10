#This program is used for calculating the mean number of amino-acid-changing mutations (except the five robustness-improving mutations) in YFP molecules which carried at least one of robustness-improving mutations (F47L, F65L, K102E, V164A or I172V) and in those YFP molecules which didn't carry any of robustness-improving mutations in each replicate population during phase I evolution
import os
import csv
import sys
import re

#the following codes are used for grouping sequences based on generation 
GList=['I-1','I-2','I-3','I-4','Anc1','Anc2','Con1','Con2']
GSea=['phase-I_1st','phase-I_2nd','phase-I_3rd','phase-I_4th',
'YFP\(ancestor\)_replicate1','YFP\(ancestor\)_replicate2','>Control_replicate1','>Control_replicate2']

def repSeq_File(filepath):
	f1 = open(filepath, "r")
	SNPlines = f1.readlines()

	Numseq=[]
	Seqlist=[]
	for i in range (len(GList)):
		Numseq.append(0) 
		Seqlist.append([])  
	for line in SNPlines:
		if not line.strip(): continue
		if re.search(">",line):
			seqName=line
		else:
			for gs in range (len(GSea)):
				if re.search(GSea[gs],seqName):
					Seqlist[gs].append(line.strip())
					Numseq[gs]=Numseq[gs]+1
	return Seqlist

#ref indicates YFP (ancestor) protein sequence
ref="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKLHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK*"

#the following codes are used for calculating the mean number of amino-acid-changing mutations (excluding the five robustness-improving mutations) in YFP molecules which carried at least one of robustness-improving mutations (F47L, F65L, K102E, V164A or I172V) in each replicate population and in each generation 
Robust=['F47L', 'F65L', 'K102E', 'V164A','I172V']
def replicate_list(seqlist):
	numSeq=len(seqlist)
	if numSeq>0:
		n=0
		mut=0
		nRobust=0
		for line in seqlist:
			if line[47-1]=="L" or line[65-1]=='L' or line[102-1]=="E" or line[164-1]=='A' or line[172-1]=='V':
				n=n+1
				for j in range(len(ref)):
					if line[j]!=ref[j]:
						mut=mut+1
				for robu in Robust:
					pos= int(robu[1:-1])
					AA=robu[-1:]
					if line[pos-1]==AA:
						nRobust=nRobust+1					
		if n>0:				
			snpNum='%.4f'%(float(mut-nRobust)/n)
		else:
			snpNum=0

	else:
		snpNum=''
	return snpNum

#the following codes are used for calculating the mean number of amino-acid-changing mutations in YFP molecules which didn't carry any of robustness-improving mutations (F47L, F65L, K102E, V164A, I172V) in each replicate population and in each generation 
def replicate_list2(seqlist):
	numSeq=len(seqlist)
	if numSeq>0:
		n=0
		mut=0
		for line in seqlist:
			if line[47-1]!="L" and line[65-1]!="L" and line[102-1]!="E" and line[164-1]!="A" and line[172-1]!="V":
				n=n+1
				for j in range(len(ref)):
					if line[j]!=ref[j]:
						mut=mut+1
		if n>0:	
			snpNum='%.4f'%(float(mut)/n)
		else:
			snpNum=0
	else:
		snpNum=''
	return snpNum

SNP_all=[]
def allSNP_File(filepath):
	nowDir = os.path.split(filepath)[0]   
	fileName = os.path.split(filepath)[1] 
	Seqlist_each=repSeq_File(filepath)
	SNP_each1=[fileName[5:-7],fileName[5:-6],'With-F']
	SNP_each2=[fileName[5:-7],fileName[5:-6],'Without-F']
	for g in range (len(GList)):
		SNP_each1.append(replicate_list(Seqlist_each[g]))
		SNP_each2.append(replicate_list2(Seqlist_each[g]))
	SNP_all.append(SNP_each1) 
	SNP_all.append(SNP_each2) 

#the following codes are used for reading all input files (in the folder "~/ProSeq") which contain protein sequences of each evolving population sequenced by SMRT sequencing
def eachFile(filepath):
	os.chdir(filepath)
	pathDir = os.listdir(filepath)      
	for s in pathDir:
		newDir=os.path.join(filepath,s)     
		if os.path.isfile(newDir) :         
			if os.path.splitext(newDir)[1]==".fasta":  
				allSNP_File(newDir)                     

eachFile("~/ProSeq")


#the following codes are used for writing the result into a csv file
with open("~/7_number_of_robustness.csv", 'wb') as csvfile:
	Wri = csv.writer(csvfile)
	Wri.writerow(['Population','Replicate','Type']+GList)
	for each in SNP_all:
			Wri.writerow(each)