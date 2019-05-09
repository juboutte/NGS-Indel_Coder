# -*- coding: us-ascii -*-

### Usage ###

#10-partitioned_nexus_files_creation.py - Generate nexus file to run IQTREE with exon/intron partitions

#input file: aligned fasta file exons_positions.txt path
#command line example: python Scripts/10-partitioned_nexus_files_creation.py Example/II-temporary_files/exon_positions.txt Example/I-input_files/fasta_files/Example.fasta Example/II-temporary_files/Example2.nex

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
#Goal of this program:

import sys

#Inputfile 1: exons_positions.txt

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

#for i in dico_exon:
#    print i, dico_exon[i]
#exit()

#Open each fasta, and identify exon and intron

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
if end<len(list_seq[0]):
    temp=[]
    temp.append(end+1)
    temp.append(len(list_seq[0]))
    list_intron.append(temp)

###

input3=open(sys.argv[3],"r")
temp=''
for line in input3:
    if "part1" in line and "charset" in line:
        temp=line.split(': ')
        #print line
#print temp
input3.close()

output=open(name1+"_partition.nex","w")
output.write("#nexus"+'\n'+'\n')
output.write("begin sets;"+'\n')
output.write(temp[0]+': ')
#DNA-Exon
j=0
while j!=len(list_exon)-1:
    output.write(str(list_exon[j][0])+'-'+str(list_exon[j][1])+'\\3,')
    output.write(str(list_exon[j][0]+1)+'-'+str(list_exon[j][1])+'\\3,')
    output.write(str(list_exon[j][0]+2)+'-'+str(list_exon[j][1])+'\\3,')
    j=j+1
output.write(str(list_exon[j][0])+'-'+str(list_exon[j][1])+'\\3,')
output.write(str(list_exon[j][0]+1)+'-'+str(list_exon[j][1])+'\\3,')
output.write(str(list_exon[j][0]+2)+'-'+str(list_exon[j][1])+'\\3;'+'\n')
#DNA-intron
intron=0
if len(list_intron)>0:
    output.write(temp[0].replace("part1","part2")+': ')
    j=0
    while j!=len(list_intron)-1:
        output.write(str(list_intron[j][0])+'-'+str(list_intron[j][1])+',')
        j=j+1
    output.write(str(list_intron[j][0])+'-'+str(list_intron[j][1])+';'+'\n')
    intron=1
#Binary
input3=open(sys.argv[3],"r")
temp=''
for line in input3:
    if "indel" in line and "charset" in line:
        temp=line.split(': ')
        #print line
#print temp
input3.close()
if len(temp)==2:
    if intron==1:
        output.write(temp[0].replace("part2","part3")+': ')
    else:
        output.write(temp[0]+': ')
    output.write(temp[1].replace('\n','')+'\n'+'end;'+'\n')
output.close()
#exit()
#create output file with name of output files created
