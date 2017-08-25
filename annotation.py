#!/usr/bin/env python
# encoding: utf-8
import sys

args = sys.argv[1:]
if len(args)<3:
    print("You should enter 3 files!")
    sys.exit()
site = open(args[0],'r')
annotation = open(arg[1],'r')
final = open(arg[2],'w')
genes=[]
for line in annotation.readlines():
    genes.append(line.strip().split('\t'))
genenames=[]
geneids=[]
for line in site.readlines():
    genename = []
    geneid = []
    temp = line.strip().split('\t')
    for gene in genes:
        if (temp[1]<=gene[2] and temp[1]>=gene[1]) or (temp[2]<=gene[2] and temp[2]>=gene[1]):
            genename.append(gene[4])
            geneid.append(gene[3])
        if genename==[]:
            genename.append('-')
        if geneid==[]:
            geneid.append('-')
        genenames.append(genename)
        geneid[len(geneid)]
        geneids.append(geneid)
site.close()
annotation.close()
p=[]
for i in range(0,len(genenames)):
    p.append(''.join(genenames[i])+'\t'+''.join(geneid[i])+'\n')
final.writelines(p)
final.close()



