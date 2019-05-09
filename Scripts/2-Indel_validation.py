# -*- coding: us-ascii -*-

import sys

### Usage ###

#2-Indel_validation.py - First script of the Part iii: Indel validation
#This script will identify insertion/deletion regions and read depth around the
#insertion/deletion

#input file: a fasta file and the path of the read depth files 

#By Julien Boutte, April 2019
#Copyright (c) 2019 Julien Boutte.
#Version 1.0.0

#This program is free software: you can redistribute it and/or modify it under
#the terms of the GNU General Public License as published by the Free Software
#Foundation, either version 3 of the License, or (at your option) any later
#version. A copy of this license is available at <http://www.gnu.org/licenses/>.
#Great effort has been taken to make this software perform its said
#task, however, this software comes with ABSOLUTELY NO WARRANTY,
#not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    
#start program
#Parsing of the fasta file
input_file=open(sys.argv[1],'r')

liste_name=[]
liste_seq=[]

seq_temp=''

for line in input_file:
    if line[0]!='>':
        temp=line.split()
        if len(temp)>0:
            seq_temp=seq_temp+temp[0]
    else:
        liste_name.append(line[1:].replace('\n','').replace('\r',''))
        if seq_temp!='':
            liste_seq.append(seq_temp)
            seq_temp=''

liste_seq.append(seq_temp) 

liste_deb=[]
liste_fin=[]

#print liste_name
#print liste_seq

#beginning and end of sequences
i=0
while i!=len(liste_seq):
    j=0
    while liste_seq[i][j]=='-' or liste_seq[i][j:j+9].count('-')!=0:
        j=j+1
    liste_deb.append(j)
    j=len(liste_seq[0])-1
    if j!='-':
        liste_fin.append(j+1)
    else:
        while liste_seq[i][j]=='-' or liste_seq[i][j:j-9].count('-')!=0:
            j=j-1
        liste_fin.append(j+1)
    i=i+1

#print liste_deb
#print liste_fin

nb_SNPs_un=0
nb_SNPs_clean=0

liste_SNPs_clean=[]


if len(liste_seq)==0:
   print "error with this sample:", sys.argv[1]
   print "script exit"
   exit()
else:
    i=0
    while i!=len(liste_seq[0]): #each position
        temp=[]
        j=0
        while j!=len(liste_seq): #each sample
            if liste_deb[j]<=i and liste_fin[j]>=i:
                #print liste_name[j], i
                temp.append(liste_seq[j][i]) 
            j=j+1
        #print temp
        if len(list(set(temp)))>1: #SNPs
            proba=list(set(temp))
            temp2=[]
            for n in proba:
                if temp.count(n)>1:
                    temp2.append(n)
            if len(temp2)>1:
                nb_SNPs_clean=nb_SNPs_clean+1 #SNPs with nucleotide present in more than 1 seq
                if temp.count('-')>=1:
                    liste_SNPs_clean.append(i)
            else:
                nb_SNPs_un=nb_SNPs_un+1 # SNP specific to one sequence
        i=i+1

    #print liste_SNPs_clean

    mini=len(liste_seq[0])
    maxi=0

    i=1
    boundaries=[]
    while i!=len(liste_SNPs_clean):
        deb=liste_SNPs_clean[i-1]
        while liste_SNPs_clean[i]-1==liste_SNPs_clean[i-1] and i!=len(liste_SNPs_clean)-1:
            i=i+1
        fin=liste_SNPs_clean[i-1]
        temp=[]
        if fin-deb>1:
            temp.append(deb)
            temp.append(fin)
            boundaries.append(temp)
            if fin-deb+1<mini:
                mini=fin-deb+1
            if fin-deb+1>maxi:
                maxi=fin-deb+1
        #print i, liste_SNPs_clean[i], liste_SNPs_clean[i-1]
        i=i+1
    #print boundaries

    #PART 2: look read depth
    output=open('temp_part1.txt','w')
    i=0
    while i!=len(boundaries): #for each insertion/deletion
        output.write('>'+sys.argv[1]+'\t'+str(boundaries[i][0])+'\t'+str(boundaries[i][1])+'\n')
        j=0
        while j!=len(liste_seq): #each sample
            if liste_deb[j]<=boundaries[i][0] and liste_fin[j]>=boundaries[i][1]: #the sequence present the insertion/deletion
                #open read depth speicific to the sequence
                file_depth=open(str(sys.argv[2])+str(liste_name[j])+'.txt','r')
                liste_depth=[]
                h=0
                for line in file_depth:
                    #print liste_name[i], h, len(liste_seq[i])
                    temp=line.replace('\n','').split('\t')
                    if liste_seq[j][h]!='N' and liste_seq[j][h]!='n':
                        liste_depth.append(temp[2])
                    else:
                        liste_depth.append('NA')
                        h=h+1
                file_depth.close()

                ###
                
                #identification of the real position of begining and end of the insertion/deletion
                pos=0
                k=0
                while k!=boundaries[i][0]: #identification of the real start position pos
                    if liste_seq[j][k]!='-':
                        pos=pos+1
                    k=k+1
                    pos2=0
                k=0
                while k!=boundaries[i][1]+1: #identification of the real end position pos2
                    if liste_seq[j][k]!='-':
                        pos2=pos2+1
                    k=k+1                
                output.write(liste_name[j]+'\t'+str(pos)+'\t'+str(pos2)+'\t'+str(liste_depth[pos:pos2])+'\t') #test to validate position
                #next step: indicate N or 0 if deletion and start 10 bp before and end 10 bp after
                if pos-10>0 and pos2+10<len(liste_seq[j].replace('-','')):
                    output.write(str(pos-10)+'\t'+str(pos2+10)+'\t'+str(liste_depth[pos-10:pos2+10])+'\n') 
                elif pos-10>0:
                    #output.write(str(pos-10)+'\t'+str(len(liste_seq[j].replace('-','')))+'\t'+str(liste_depth[pos-10:len(liste_seq[j].replace('-',''))])+'\n')
                    output.write(str(pos-10)+'\t'+str(len(liste_seq[j].replace('-',''))-1)+'\t'+str(liste_depth[pos-10:])+'\n')
                else:
                    output.write(str(0)+'\t'+str(pos2+10)+'\t'+str(liste_depth[0:pos2+10])+'\n')
                ###
                
            j=j+1
        i=i+1
    output.close()



