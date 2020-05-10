#This program is used for calculating frequencies of genotypes E18G, E18G+F and the ratio between them (r) in the replicate population S2, of genotypes M79L, M79L+F and the ratio between them (r) in the replicate population N4, of genotypes R110S, R110S+F and the ratio between them (r) in the replicate population S1, and of genotypes N145S, N145S+F and the ratio between them (r) in the replicate population S1 during phase II evolution. Here F refers to any of the five foldability-improving mutations F47L, F65L, K102E, V164A and I172V.
import csv
import sys
import re
import os
import ntpath

#the following codes are used for grouping sequences based on generation
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
		
#the following codes are used for calculating frequencies of genotypes X(E18G,M79L,R110S or N145S), X+F and the ratio between them (r) in the corresponding replicate population
GenoType=['E18G','M79L','R110S','N145S']
Rep_list=['S2','N4','S1','S1']
Inpath="~/ProSeq/"

def replicate_list(Inx,seqlist):
	numSeq=len(seqlist)
	frequency_list=[]	
	if numSeq>0:		
		pos=int(GenoType[Inx][1:-1])
		AA=GenoType[Inx][-1:]
		X=0
		X_F=0
		for line in seqlist:
			if line[pos-1]==AA:
				X=X+1
				if line[47-1]=="L" or line[65-1]=='L' or line[102-1]=="E" or line[164-1]=='A' or line[172-1]=='V':
					X_F=X_F+1

		fre_X=100.0*X/numSeq
		fre_XF=100.0*X_F/numSeq
		if  fre_X>0:
			r=100.0*fre_XF/fre_X
		else:
			r=0

		frequency_list.append([GenoType[Inx],"%.4f"%fre_X])
		frequency_list.append([GenoType[Inx]+'+F',"%.4f"%fre_XF])
		frequency_list.append(['r',"%.4f"%r])
	else: 
		frequency_list=[[GenoType[Inx],''],[GenoType[Inx]+'+F',''],['r','']]
	return frequency_list


SNP_all=[]
def allSNP_File(Inx):
	filepath=os.path.join(Inpath,'2_aa_'+Rep_list[Inx]+'.fasta')     
	nowDir = os.path.split(filepath)[0]   
	fileName = os.path.split(filepath)[1] 
	Seqlist_each=repSeq_File(filepath)
	SNP_each=[]
	for g in range (len(GList)):
		SNP_gx=replicate_list(Inx,Seqlist_each[g])
		for n in range (len(SNP_gx)):
 			SNP_gx[n]=[fileName[5:-7],fileName[5:-6],GList[g]]+SNP_gx[n]
		SNP_each.append(SNP_gx)
	SNP_all.append(SNP_each)

for i in range(len(Rep_list)):
        allSNP_File(i)

#the following codes are used for writing the result into a csv file
with open("~/11_frequencyOf_X.csv", 'wb') as csvfile:
	Wri = csv.writer(csvfile)
	Wri.writerow(['Population','Replicate','Generation','Genotype','Frequency'])
	for each in SNP_all:
		for g in each:
			for j in g:
				Wri.writerow(j)