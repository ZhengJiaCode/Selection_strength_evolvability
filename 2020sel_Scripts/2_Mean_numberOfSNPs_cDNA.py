#This program is used for calculating the mean number of SNPs per YFP molecule sequenced by SMRT in each replicate population and in each generation
import os
import csv
import sys
import re

#ref indicates YFP (ancestor) cDNA sequence
ref="ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCTTCGGCTACGGCCTGCAATGCTTCGCCCGCTACCCCGACCACATGAAGCTGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCTACCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTGA"

#the following codes are used for grouping sequences in a replicate population based on generation 
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

#the following codes are used for getting the mean number of SNPs per YFP molecule in each replicate population and in each generation 
def SNP_list(seqlist):
	numSeq=len(seqlist)
	if numSeq>0:
		diff=0
		for seq in seqlist:
			for i in range (len(ref)):
				if seq[i]!=ref[i]:
					diff=diff+1
		snpNum='%.4f'%(float(diff)/numSeq)
	else:  
		snpNum=''

	return snpNum

SNP_all=[]
def allSNP_File(filepath):
	nowDir = os.path.split(filepath)[0]   
	fileName = os.path.split(filepath)[1] 
	Seqlist_each=repSeq_File(filepath)
	SNP_each=[fileName[0:-7],fileName[0:-6]]
	for g in range (len(GList)):
		SNP_each.append(SNP_list(Seqlist_each[g]))
	SNP_all.append(SNP_each) 

#the following codes are used for reading all input files (in the folder "~/2020sel_smrt_seq") which contain cDNA sequences of each evolving population sequenced by SMRT sequencing
def eachFile(filepath):
	os.chdir(filepath)
	pathDir = os.listdir(filepath)      
	for s in pathDir:
		newDir=os.path.join(filepath,s)     
		if os.path.isfile(newDir) :         
			if os.path.splitext(newDir)[1]==".fasta":  
				allSNP_File(newDir)                     

eachFile("~/2020sel_smrt_seq")

#the following codes are used for writing the result into a csv file
with open("~/2_Mean_numberOfSNPs_cDNA.csv", 'wb') as csvfile:
	Wri = csv.writer(csvfile)
	Wri.writerow(['Population','Replicate']+GList)
	for each in SNP_all:
			Wri.writerow(each)