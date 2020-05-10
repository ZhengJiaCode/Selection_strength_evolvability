#This program is used for calculating genetic diversity (Average pairwise sequence distance) of variants in each replicate population and in each generation
import os
import csv
import sys
import re

#the following codes are used for grouping sequences based on generation 
GList=['I-1','I-2','I-3','I-4','II-1','II-2','II-3','II-4','Anc1','Anc2','Con1','Con2']
GSea=['phase-I_1st','phase-I_2nd','phase-I_3rd','phase-I_4th','phase-II_1st','phase-II_2nd','phase-II_3rd','phase-II_4th',
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
#the following codes are used for calculating genetic diversity (Average pairwise sequence distance) of variants in each replicate population and in each generation
def diversity_list(seqlist):
	numSeq=len(seqlist)
	diversity=[]
	if numSeq!=0:
		TPD=0
		for i in range (0,720):
			DiMatrix=[]
			for line in seqlist:
				if not re.search(">",line):
					DiMatrix.append(list(line)[i])
			y=list(DiMatrix)
			NumA=y.count("A")
			NumT=y.count("T")
			NumC=y.count("C")
			NumG=y.count("G")
			Dj = NumA * (NumT+NumC+NumG)+NumC * (NumT+NumG)+NumT * NumG
			PDj = Dj/(numSeq*(numSeq-1)/2.0)
			TPD=TPD+PDj
		APD=float(TPD)/720
		DIR='%.6f'%APD
	else:
		DIR=''
	return DIR

Pop_all=[]
def Population_File(filepath):
	nowDir = os.path.split(filepath)[0]   
	fileName = os.path.split(filepath)[1] 
	Seqlist_each=repSeq_File(filepath)
	Pop_each=[fileName[0:-7],fileName[0:-6]]
	for g in range (len(GList)):
		Pop_each.append(diversity_list(Seqlist_each[g]))
	Pop_all.append(Pop_each) 

#the following codes are used for reading all input files (in the folder "~/2020sel_smrt_seq") which contain cDNA sequences of each evolving population sequenced by SMRT sequencing
def eachFile(filepath):
	os.chdir(filepath)
	pathDir = os.listdir(filepath)      
	for s in pathDir:
		newDir=os.path.join(filepath,s)     
		if os.path.isfile(newDir) :         
			if os.path.splitext(newDir)[1]==".fasta":  
				Population_File(newDir)                     

eachFile("~/2020sel_smrt_seq")


#the following codes are used for writing the result into a csv file
with open("~/5_genetic_diversity.csv", 'wb') as csvfile:
	Wri = csv.writer(csvfile)
	Wri.writerow(['Population','Replicate']+GList)
	for each in Pop_all:
			Wri.writerow(each)