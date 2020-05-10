#This program is for determining the frequency of nonsense mutations in each replicate population during phase I evolution
import os
import csv
import sys
import re

#the following codes are used for grouping sequences based on generation 
GList=['I-1','I-2','I-3','I-4']
GSea=['phase-I_1st','phase-I_2nd','phase-I_3rd','phase-I_4th']

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

#the following codes are used for determining the frequency of nonsense mutations in each replicate population during each generation of phase I evolution
def roundNonGenoFre_list(seqlist):
	numseq=len(seqlist)
	Fre=''
	if numseq>0:
		Num=0
		for seq in seqlist:
			for j in range(len(ref)-1):
				if seq[j]=="*":
					Num=Num+1	
		fre=float(Num)/numseq*100
		Fre="%.4f"%fre
	return Fre

Rep_all=[]
def allRep_File(filepath):
	fileName = os.path.split(filepath)[1]
	Seqlist_each=repSeq_File(filepath)
	Fre_14=[]
	for g in range (len(GList)):
		Rep_gx=roundNonGenoFre_list(Seqlist_each[g])
		Fre_14.append(Rep_gx)
	Rep_each=[fileName[5:-7],fileName[5:-6],'Nonsense']+Fre_14
	Rep_all.append(Rep_each)

#the following codes are used for reading all input files (in the folder "~/ProSeq") which contain protein sequences of each evolving population sequenced by SMRT sequencing
def eachFile(filepath):
	os.chdir(filepath)
	pathDir = os.listdir(filepath)
	for s in pathDir:
		newDir=os.path.join(filepath,s)    
		if os.path.isfile(newDir) :         
			if os.path.splitext(newDir)[1]==".fasta":  
				allRep_File(newDir)                                    
eachFile("~/ProSeq")


#the following codes are used for writing the result into a csv file
with open("~/10_frequencyOfNonsenseMutations.csv", 'wb') as csvfile:
	Wri = csv.writer(csvfile)
	Wri.writerow(['Population','Replicate','Mutation']+GList)
	for each in Rep_all:
			Wri.writerow(each)