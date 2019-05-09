# -*- coding: us-ascii -*-

### Usage ###

#9-identification_boundaries.py - Generate txt file with exon positions in the aligned fasta file

#input file: bastn result
#command line example: python ../Scripts/9-identification_boundaries_ONGOING.py res_blastn.txt name_sample

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

#Goal of this program: using blast result, extract position of the introns
#and create a list of introns
#delete redundancy and create new for overlap

inputfile=open(sys.argv[1],'r')
output=open('exon_positions.txt','w')

list_introns=[]

for line in inputfile:
    if line[0]!='#':
        temp=line.split('\t')
        tempo=[]
        if int(temp[6])<int(temp[7]):
            tempo.append(int(temp[6]))
            tempo.append(int(temp[7]))
        else:
            tempo.append(int(temp[7]))
            tempo.append(int(temp[6]))
        list_introns.append(tempo)

#print list_introns

test=0
while test==0:
    test=1
    i=0
    while i<=len(list_introns)-2:
        j=0
        while j<=len(list_introns)-1:
            if i!=j:
                #print i,j, len(list_introns)
                #print list_introns[i],list_introns[j]
                if list_introns[i][0]<=list_introns[j][0] and list_introns[i][1]>=list_introns[j][1]:
                    del(list_introns[j])
                    test=0
                elif list_introns[j][0]<=list_introns[i][0] and list_introns[j][1]>=list_introns[i][1]:
                    del(list_introns[i])
                    test=0
                    j=0
                elif list_introns[i][0]<=list_introns[j][0] and list_introns[j][0]<=list_introns[i][1]:
                    tempo=[]
                    tempo.append(min(list_introns[i][0],list_introns[j][0]))
                    tempo.append(max(list_introns[i][1],list_introns[j][1]))
                    #print list_introns[i], list_introns[j], tempo
                    if i<j:
                        del(list_introns[j])
                        del(list_introns[i])
                    else:
                        del(list_introns[i])
                        del(list_introns[j])
                    j=0
                    i=0
                    list_introns.append(tempo)
                    test=0
                elif list_introns[j][0]<=list_introns[i][0] and list_introns[i][0]<=list_introns[j][1]:
                    tempo=[]
                    tempo.append(min(list_introns[i][0],list_introns[j][0]))
                    tempo.append(max(list_introns[i][1],list_introns[j][1]))
                    #print list_introns[i], list_introns[j], tempo
                    if i<j:
                        del(list_introns[j])
                        del(list_introns[i])
                    else:
                        del(list_introns[i])
                        del(list_introns[j])
                    j=0
                    i=0
                    list_introns.append(tempo)
                    test=0
                else:
                    j=j+1
            else: j=j+1
        i=i+1


#reorder list_introns
list_temp=[]
for i in list_introns:
    list_temp.append(i[0])
list_temp.sort()

#print list_temp

list_introns_finale=[]

i=0
while i!=len(list_temp):
    j=0
    while list_introns[j][0]!=list_temp[i]:
        j=j+1
    list_introns_finale.append(list_introns[j])
    i=i+1

output.write(str(sys.argv[2])+'\n')
output.write(str(list_introns_finale)+'\n')
output.close()
