# -*- coding: us-ascii -*-

### Usage ###

#8-nexus_files_creation.py - Generate IQTREE file for DNA and indel characters

#input file: Two files created by 2MATRIX software
#command line example: python Scripts/8-nexus_files_creation.py Example/II-temporary_files/Example2_dna.phy Example/II-temporary_files/Example_indel.phy MyFolder_T20/

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

#Using Xaligned-coded.part create nexus file requiert for iq-tree

inputfile1=open(sys.argv[1],'r') # _dna.phy
inputfile2=open(sys.argv[2],'r') # _indel.phy

output=open(sys.argv[1].replace('_dna.phy','.nex'),'w')

try:
    path=sys.argv[3]
except:
    path=''
    
end_dna=0
end_indel=0

test=0
for line in inputfile1:
    if test==0:
        temp=line.replace('\n','').split(' ')
        end_dna=int(temp[1])        
        break
inputfile1.close()

test=0
for line in inputfile2:
    if test==0:
        temp=line.replace('\n','').split(' ')
        end_indel=int(temp[1])        
        break
inputfile2.close()

temp=sys.argv[1].split('/')
name1=temp[-1]
temp=sys.argv[2].split('/')
name2=temp[-1]

#output write
output.write('#nexus'+'\n'+'\n')
output.write('begin sets;'+'\n')
output.write('	charset part1 = '+str(path)+name1+': '+str(1)+'-'+str(end_dna)+';'+'\n')
output.write('	charset part2 = '+str(path)+name2+': '+str(1)+'-'+str(end_indel)+';'+'\n'+'\n')
#output.write('	charpartition mine = GTR+G:part1, MK+ASC:part2;'+'\n')
output.write('end;'+'\n')
output.close()
