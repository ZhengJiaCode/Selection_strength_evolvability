#This program is used for converting cDNA sequences of YFP molecules into protein sequences
import re
import sys
import os
import ntpath

codon={'GCT':'A','GCC':'A','GCA':'A', 'GCG':'A','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R','AAT':'N','AAC':'N','GAT':'D','GAC':'D','TGT':'C','TGC':'C','CAA':'Q','CAG':'Q','GAA':'E', 'GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G','CAT':'H','CAC':'H','ATT':'I','ATC':'I','ATA':'I','ATG':'M','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L','AAA':'K','AAG':'K','TTT':'F', 'TTC':'F','CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P','TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S','ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 	'TGG':'W','TAT':'Y', 'TAC':'Y','GTT':'V','GTC':'V','GTA':'V','GTG':'V','TAA':'*','TGA':'*','TAG':'*'}

OutpathI='~/ProSeq/'
def Protein_File(filepath):
	f1 = open(filepath, "r")
	SNPlines = f1.readlines()
	
	FileName=ntpath.basename(filepath) 
	newFile=os.path.join(OutpathI,'3_aa_'+FileName)     
	openfile=open(newFile,"w")
	openfile=open(newFile,"w")
	for line in SNPlines:
		aa_seq=[]
		if re.search(">",line):
			openfile.write("%s"%line)
		elif len(line)==721:
			for i in range (0,720,3):
				x=line[i:i+3]
				aa_seq.append(codon[x])
				join_aa="".join(aa_seq)
			openfile.write("%s\n"%join_aa)

def eachFile(filepath):
	os.chdir(filepath)
	pathDir = os.listdir(filepath)      
	for s in pathDir:
		newDir=os.path.join(filepath,s)     
		if os.path.isfile(newDir) :         
			if os.path.splitext(newDir)[1]==".fasta":  
				Protein_File(newDir)                     

eachFile("~/2020sel_smrt_seq")







