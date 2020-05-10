#This program is used for calculating the frequency of sequences carrying both G66S and Y204 during phase II evolution
import os
import csv
import sys
import re

GList=['II-1','II-2','II-3','II-4']
GSea=['phase-II_1st','phase-II_2nd','phase-II_3rd','phase-II_4th']

#the following codes are used for grouping sequences from each generation of phase II evolution 

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

#the following codes are used for calculating frequency of YFP molecules which carried G66S and Y204 in each replicate population during each generation of phase II evolution 
def SNP_list(seqlist):
	numSeq=len(seqlist)
	freList=[]

	if numSeq>0:
		n=0
		for line in seqlist:
			if line[65]=="S" and line[203]=='C':
				n=n+1
		fre="%.4f" %(float(100.0*n)/numSeq)
	else:  
		fre=''

	return fre



SNP_all=[]
def allSNP_File(filepath):
	nowDir = os.path.split(filepath)[0]   
	fileName = os.path.split(filepath)[1] 
	Seqlist_each=repSeq_File(filepath)
	SNP_each=[fileName[5:-7],fileName[5:-6],'U']

	for g in range (len(GList)):
		SNP_each.append(SNP_list(Seqlist_each[g]))

	SNP_all.append(SNP_each) 


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
with open("~/7_frequency_of_U.csv", 'wb') as csvfile:
	Wri = csv.writer(csvfile)
	Wri.writerow(['Population','Replicate','Mutation']+GList)
	for each in SNP_all:
			Wri.writerow(each)