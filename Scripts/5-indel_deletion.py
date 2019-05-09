# -*- coding: us-ascii -*-

import sys

### Usage ###

#5-indel_deletion.py - Generate temporary fasta file
#delete insertion/deletion regions that must not be coded as binary characters

#input file: MRD2_MRD3_TX.txt, command line: python 5-indel_deletion.py MRD2_MRD3_TX.txt path_to_fasta_file

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

#Function to ordonate windows for each alignment fasta file
def tri_dico(liste):
    if len(liste)==1:
        return liste
    else:
        tempo=[]
        for i in liste:
            tempo.append(i[0])
        order=sorted(tempo, reverse=True) ###
        if len(list(sorted(order)))!=len(liste):
            print 'ERROR FUNCTION 1'
        final=[]
        i=0
        while i!=len(order):
            for j in liste:
                if j[0]==order[i]:
                    final.append(j)
            i=i+1
        if len(final)==len(liste):
            return final
        else:
            print 'ERROR FUNCTION 2'
#Part 1: extraction in .txt file of the different windows
       
inputfile=open(sys.argv[1],'r') #mean-max-close_D2_GAP_NoGAP_[with,without].txt
#folder=sys.argv[2] #words with or without only

dico_seq={}

for line in inputfile:
    temp=line.split('\t')
    name=temp[0].replace('\n','').replace('\r','').split('/')
    if name[-1] not in dico_seq:
        dico_seq[name[-1]]=[]
    liste_tempo=[]
    liste_tempo.append(int(temp[1]))
    liste_tempo.append(int(temp[2]))
    dico_seq[name[-1]].append(liste_tempo)
    liste_tempo=[]

inputfile.close()

#print len(dico_seq)
#print dico_seq

#Part 2: for each element of dico_seq:
#open fasta file and modify it

#Warning: necessary to ordonate windows by descendant values

#print dico_seq

for i in dico_seq:
    i2=i.replace('\n','').replace('\r','').split('/')
    #organise dico_seq[i]
    tempo=[]
    for j in dico_seq[i]:
        tempo.append(j)
    tempo2=tri_dico(tempo)
    #tempo2 is the descending order of dico_seq
    output=open(i.replace('.fasta','_temp.fasta'),'w')
    liste_name=[]
    liste_seq=[]
    seq_temp=''
    inputfile=open(str(sys.argv[2])+str(i2[-1]),'r')
    for line in inputfile:
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
    #outputfile
    j=0
    while j!=len(liste_seq):
        seq_temp=liste_seq[j]
        for k in tempo2:
            seq_temp=seq_temp[:k[0]]+seq_temp[k[1]+1:]
        output.write('>'+liste_name[j]+'\n'+seq_temp+'\n')
        seq_temp=''
        j=j+1
    output.close()



