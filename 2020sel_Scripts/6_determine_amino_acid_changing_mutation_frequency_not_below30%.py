#This program is used for identifying amino-acid changing mutations with frequency >=30% at the end of phase II evolution and for culculating their frequencies in each replicate population during phase II evolution
import os
import csv
import sys
import re

#ref indicates YFP (ancestor) protein sequence
ref="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKLHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK*"

#the following codes are used for grouping sequences from each generation of phase II evolution 
GList=['II-1','II-2','II-3','II-4']
GSea=['phase-II_1st','phase-II_2nd','phase-II_3rd','phase-II_4th']

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


#the followig codes are used for idenfitying high-frequency amino-acid changing mutations which reached >=30% in at least one replicate population at the end of phase II (generation 4)
def roundMut_list(seqlist):
	mutSNPx_gx=[]
	numseq=len(seqlist)
	if numseq>0:
		mut_pos=[]	
		for j in range (len(ref)):

			for i in range (len(seqlist)):
				if seqlist[i][j]!=ref[j]:
					mut_pos.append(ref[j]+str(j+1)+seqlist[i][j])

		for mut in set(mut_pos):
			Num=mut_pos.count(mut) 
			fre=float(100.0*Num)/numseq
			if fre>=30:
				mutSNPx_gx.append(mut)
	return mutSNPx_gx


def PosMutsort(mutlist):
	PosMut=[]
	PosMutsort=[]
	for mut in set(mutlist):
		PosMut.append([mut[1:-1],mut])
	for i in range (len(ref)):
		for j in range (len(PosMut)):
			if PosMut[j][0]==str(i):
				PosMutsort.append(PosMut[j])
	return PosMutsort


Mut_all=[]
def allMut_File(filepath):
	Seqlist_each=repSeq_File(filepath)
	Mut_each=[]

	Mut_gx=roundMut_list(Seqlist_each[4-1])
	for i in range (len(Mut_gx)):
		Mut_each.append(Mut_gx[i])
	for i in range (len(Mut_each)):
		Mut_all.append(Mut_each[i])

#the following codes are used for reading all input files (in the folder "~/ProSeq") which contain protein sequences of each evolving population sequenced by SMRT sequencing
def eachFile_Mut(filepath):
	os.chdir(filepath)
	pathDir = os.listdir(filepath)
	for name in pathDir:
		newDir=os.path.join(filepath,name)    
		if os.path.isfile(newDir) :         
			if os.path.splitext(newDir)[1]==".fasta":  
				allMut_File(newDir) 				            

eachFile_Mut("~/ProSeq")


#####################################
#the followig codes are used for calculating frequencies of high-frequency (>=30%) amino-acid changing mutations in each replicate population during phase II evolution
def roundFreSel_list(seqlist):
	Mut_sel=PosMutsort(Mut_all)
	mutsel_fre=[]
	mut_gx=[]
	numseq=len(seqlist)
	if numseq>0:
		for i in range (len(Mut_sel)):
			mutsel_fre.append(Mut_sel[i]+['0'])
		mut_pos=[]	
		for j in range (len(ref)):

			for seq in seqlist:
				if seq[j]!=ref[j]:
					mut_pos.append(ref[j]+str(j+1)+seq[j])
		for mut in set(mut_pos):
			Num=mut_pos.count(mut)
			fre="%.4f"%(float(100.0*Num)/numseq)
			mutFre=[mut[1:-1],mut,fre]
			mut_gx.append(mutFre)
	else:
		for i in range (len(Mut_sel)):
			mutsel_fre.append(Mut_sel[i]+[''])

	for sel in mutsel_fre:
		for real in mut_gx:
			if sel[1]==real[1]:
				sel[2]=real[2]
	return mutsel_fre


SNP_all=[]
def allSNP_File(filepath):
	nowDir = os.path.split(filepath)[0]  
	fileName = os.path.split(filepath)[1] 
	Seqlist_each=repSeq_File(filepath)
	SNP_each=[]
	for g in range (len(GList)):
		SNP_gx=roundFreSel_list(Seqlist_each[g])
		for num in range (len(SNP_gx)):
 			SNP_gx[num]=[fileName[5:-7],fileName[5:-6],GList[g]]+SNP_gx[num]  
		SNP_each.append(SNP_gx)
	SNP_each0New=[]
	for j in range (len(SNP_each[0])):
		newMut0=SNP_each[0][j][0:2]+SNP_each[0][j][3:5] 
		fre0_14=[]
		for Rg in SNP_each:
			fre0_14.append(Rg[j][5])
		Mut_fre0_14=newMut0+fre0_14
		SNP_each0New.append(Mut_fre0_14)
	SNP_all.append(SNP_each0New) 

#the following codes are used for reading all input files (in the folder "~/ProSeq") which contain protein sequences of each evolving population sequenced by SMRT sequencing
def eachFile_SNP(folderpath):
	os.chdir(folderpath)
	pathDir = os.listdir(folderpath) 
	for name in pathDir:
		newDir=os.path.join(folderpath,name)
		if os.path.isfile(newDir) :         
			if os.path.splitext(newDir)[1]==".fasta":  
				allSNP_File(newDir)                             


eachFile_SNP("~/ProSeq")

#the following codes are used for writing the results into a csv file

with open("~/6_indentify_amino_acid_changing_mutation_frequency30II4.csv", 'wb') as csvfile:
	Wri = csv.writer(csvfile)
	Wri.writerow(['Population','Replicate','Position','Mutation']+GList)
	for each in SNP_all:
		for gx in each:
			Wri.writerow(gx)