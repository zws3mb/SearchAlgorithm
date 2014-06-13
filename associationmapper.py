__author__ = 'seidz'
#! /usr/bin/env python
#Import the string functions from python
import string

firsts = dict()
lasts = dict()
# Open the input text file for reading
dataFile = open('sample.txt','r')#open('C:\Users\Agilex01\Projects\Project 2-Namematching\DataSet\queries.txt', 'r') #sys.stdin.readlines()

# Loop through each line of the input data file
for eachLine in dataFile:
# setup a temporay variable 1|ANN P|DANIELS
    index, first, last = eachLine.split('|')
    if (first.find(' ')!=-1):
        firsts[index]=first.split(' ')[0].strip().replace('-','').replace('\'','') #Just takes first names, for now. We'll see if there's a difference in the distributions.
        #firsts[index]=temp.strip()
    else:
        firsts[index]=first.strip().replace('-','').replace('\'','').replace('\\','').replace('\/','').replace('"','').replace(',','')
    lasts[index]=last.strip().replace('-','').replace('\'','').replace('\\','').replace('\/','').replace('"','').replace(',','')
dataFile.close()
outfile = open('sampleassociation0.txt','w')
outfile2= open('sampleassociation0last.txt','w')
for name in firsts.values():
    for i in range(len(name)-1):
        if(name[i+1:i+3]!=''):
            outfile.write('%s:%s\t%s\n' % (name[i:i+2],name[i+1:i+3],"1"))
outfile.close()
for lnam in lasts.values():
    for i in range(len(lnam)-1):
        if(lnam[i+1:i+3]!=''):
            outfile2.write('%s:%s\t%s\n' % (lnam[i:i+2],lnam[i+1:i+3],"1"))

