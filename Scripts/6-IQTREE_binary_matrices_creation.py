# -*- coding: us-ascii -*-

### Usage ###

#6-IQTREE_binary_matrices_creation.py - Generate IQTREE file for binary characters

#input file: Two files created by 2MATRIX software
#command line example: python 6-IQTREE_binary_matrices_creation.py Example.phy Example.part 

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

#using files created by 2MATRIX create 2 .phy for IQ-tree
input_phy=open(sys.argv[1],'r')
input_part=open(sys.argv[2],'r')

#output_phy1=open(sys.argv[1].replace('.phy','')+'_dna.phy','w')
output_phy2=open(sys.argv[1].replace('.phy','')+'_indel.phy','w')

start_dna=0
end_dna=0

start_indel=0
end_indel=0

test=0
for line in input_part:
    if test==0:
        temp=line.replace('\n','').split('=')
        temp2=temp[1].split('-')
        start_dna=int(temp2[0])
        end_dna=int(temp2[1])        
        test=1
    else:
        temp=line.replace('\n','').split('=')
        temp2=temp[1].split('-')
        start_indel=int(temp2[0])
        end_indel=int(temp2[1])

input_part.close()

list_supp=[]

test=0
for line in input_phy:
    if test==0:
        test=1
    else:
        temp=line.split(' ')
        name=temp[0].split('_')
        #if temp[-1][start_indel-1:].count('0')==0 and temp[-1][start_indel-1:].count('1')==0 and temp[-1][start_indel-1:].count('-')==0:
        if len(temp[-1][start_indel-1:].replace('\n',''))==temp[-1][start_indel-1:].replace('\n','').count('?'):
            list_supp.append(name[0])

print list_supp

input_phy.close()
input_phy=open(sys.argv[1],'r')

test=0
for line in input_phy:
    if test==0:
        test=1
        temp=line.split(' ')
        #output_phy1.write(str(int(temp[0])-len(list_supp))+' '+str(end_dna-start_dna+1)+'\n')
        output_phy2.write(str(int(temp[0])-len(list_supp))+' '+str(end_indel-start_indel+1)+'\n')
    else:
        temp=line.split(' ')
        name=temp[0].split('_')
        if name[0] not in list_supp:
            #output_phy1.write(name[0]+' ')
            #output_phy1.write(temp[-1][:start_indel-1].replace('\n','')+'\n')
            output_phy2.write(name[0]+' ')
            output_phy2.write(temp[-1][start_indel-1:].replace('\n','')+'\n')

input_phy.close()
#output_phy1.close()
output_phy2.close()
