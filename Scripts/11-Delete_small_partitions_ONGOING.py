# -*- coding: us-ascii -*-

### Usage ###

#11-Delete_small_partitions.py - Delete small partition and generate a new fasta file

#input file: list of exon, fasta file, length of small partition and name of the new fasta
#command line example: 11-Delete_small_partitions_ONGOING.py Example/II-temporary_files/exon_positions.txt Example/I-input_files/fasta_files/Example.fasta Example_partitioned.fasta 100

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


import sys

#1- Inputfile 1
input1=open(sys.argv[1],"r")
dico_exon={}

test=0
for line in input1:
    if test==0:
        test=1
        temp=line.replace('\n','')
        temp=temp.replace('.txt','')
        dico_exon[temp]=[]
        name1=temp
    else:
        test=0
        temp2=line.replace('\n','')
        temp2=temp2.replace('[','')
        temp2=temp2.replace(']','')
        temp2=temp2.replace(' ','')
        temp2=temp2.split(',')
        list_pos=[]
        plop=[]
        j=0
        test2=0
        while j!=len(temp2):
            if len(plop)==2:
                list_pos.append(plop)
                plop=[]
            plop.append(int(temp2[j]))
            j=j+1
        list_pos.append(plop)
        dico_exon[temp]=list_pos
        
input1.close()

input_file=open(sys.argv[2],'r')

list_name=[]
list_seq=[]

seq_temp=''

for line in input_file:
    if line[0]!='>':
        temp=line.split()
        if len(temp)>0:
            seq_temp=seq_temp+temp[0]
    else:
        list_name.append(line[1:].replace('\n','').replace('\r',''))
        if seq_temp!='':
            list_seq.append(seq_temp)
            seq_temp=''

list_seq.append(seq_temp) 
    
lg=len(list_seq[0])

#2- Identification list_intron and list_exon
        
#check start

list_intron=[]
list_exon=[]
if dico_exon[name1][0][0]!=1:
    temp=[]
    temp.append(1)
    temp.append(dico_exon[name1][0][0]-1)
    list_intron.append(temp)
temp2=[]
j=0
while j!=len(dico_exon[name1]):
    temp=[]
    temp.append(dico_exon[name1][j][0])
    temp.append(dico_exon[name1][j][1])
    list_exon.append(temp)
    j=j+1
    
end=dico_exon[name1][-1][1]

if len(dico_exon[name1])>1:
    j=0
    tempo=[]
    while j!=len(dico_exon[name1]):
        tempo.append(dico_exon[name1][j][0])
        tempo.append(dico_exon[name1][j][1])
        j=j+1
    del(tempo[-1])
    del(tempo[0])
    j=0
    while j!=len(tempo):
        temp=[]
        temp.append(tempo[j]+1)
        temp.append(tempo[j+1]-1)
        list_intron.append(temp)
        j=j+2
if end<lg:
    temp=[]
    temp.append(end+1)
    temp.append(lg)
    list_intron.append(temp)

#3- identify small introns and small exons

list_tot=[]
ordo=[]
for i in list_exon:
    if int(i[1])<int(i[0]): plop=1
    else:
        temp=[]
        temp.append(int(i[0]))
        temp.append(int(i[1]))
        ordo.append(int(i[0]))
        temp.append('E')
        list_tot.append(temp)
    
for i in list_intron:
    if int(i[1])<int(i[0]): plop=1
    else:
        temp=[]
        temp.append(int(i[0]))
        temp.append(int(i[1]))
        ordo.append(int(i[0]))
        temp.append('I')
        list_tot.append(temp)

#Ordonate list_tot
ordo_ordo=sorted(ordo)

list_tot_ordo=[]
for i in ordo_ordo:
    for j in list_tot:
        if j[0]==i:
            list_tot_ordo.append(j)

if len(list_tot)!=len(list_tot_ordo):
    print 'error length, exit'
    print list_tot
    print list_tot_ordo
    exit()

#delete small regions in list_tot_ordo
list_delete_fasta=[] #for fasta file

save=[]
for i in list_tot_ordo:
    save.append(i)

try:
    length=int(sys.argv[4])
except:
    length=100

i=0
while i!=len(list_tot_ordo):
    if list_tot_ordo[i][1]-list_tot_ordo[i][0]+1<length:
        #for fasta file
        temp=[]
        temp.append(list_tot_ordo[i][0]-1) ##### ATTENTION -1 #####
        temp.append(list_tot_ordo[i][1]-1) ##### ATTENTION -1 #####
        list_delete_fasta.append(temp) 
        #change values
        supp=list_tot_ordo[i][1]-list_tot_ordo[i][0]+1
        j=i+1
        while j<len(list_tot_ordo):
            temp=[]
            temp.append(list_tot_ordo[j][0]-supp)
            temp.append(list_tot_ordo[j][1]-supp)
            temp.append(list_tot_ordo[j][2])
            list_tot_ordo[j]=temp
            j=j+1
        del(list_tot_ordo[i])
    else:
        i=i+1

#create list_final_exon containing new positions of exons
list_final_exon=[]

for i in list_tot_ordo:
    if i[2]=='E':
        temp=[]
        temp.append(i[0])
        temp.append(i[1])
        list_final_exon.append(temp)

#Using list_delete_fasta, list_name and list_seq
#create a new fasta file

list_retained=[]

if len(list_delete_fasta)>0:
    if list_delete_fasta[0][0]!=0:
        temp=[]
        temp.append(0)
        temp.append(list_delete_fasta[0][0]-1)
        list_retained.append(temp)
    i=0
    while i!=len(list_delete_fasta)-1:
        temp=[]
        temp.append(list_delete_fasta[i][1]+1)
        temp.append(list_delete_fasta[i+1][0]-1)
        list_retained.append(temp)
        i=i+1

    if list_delete_fasta[i][1]!=lg-1:
        temp=[]
        temp.append(list_delete_fasta[i][1]+1)
        temp.append(lg-1)
        list_retained.append(temp)

    output=open(sys.argv[3],'w')

    i=0
    while i!=len(list_name):
        output.write('>'+list_name[i]+'\n')
        j=0
        while j!=len(list_retained):
            output.write(list_seq[i][list_retained[j][0]:list_retained[j][1]])
            j=j+1
        output.write('\n')
        i=i+1
    output.close()

else:
    
    output=open(sys.argv[3],'w')

    i=0
    while i!=len(list_name):
        output.write('>'+list_name[i]+'\n')
        output.write(list_seq[i]+'\n')
        i=i+1
    output.close()































