__author__ = 'seidz'
names=open('sample.txt','r')
firsts=dict()
# for eachLine in names:
#     index, first, last = eachLine.split('|')
#     if (first.find(' ')!=-1):
#         firsts[index]=first.split(' ')[0].strip().replace('-','').replace('\'','') #Just takes first names, for now. We'll see if there's a difference in the distributions.
#         #firsts[index]=temp.strip()
#     else:
#         firsts[index]=first.strip().replace('-','').replace('\'','')
# outfile=open('namelengths.txt','w')
# for name in firsts.values():
#     outfile.write('%s\n'%(len(name)))
#print 1001211311 |1001211411
def onecount(value):
    count=0;
    while(value!=0):
        value = value&(value-1);
        count=count+1;
        print bin(value)
    return count
# print '\t%b'%(0b101)
# print '\t%b'%(0b01011)
# print bin((int('01011',2)))
# print bin(int(0b101)^int(0b01011))
# combo=dict()
# firsts=open('sampleoutput2.txt')
# lasts=open('sampleoutput2lasts.txt')
# for line in firsts:
#     kmer,freq=line.split('\t')
#     try:
#         freq=int(freq)
#         if kmer in combo.keys():
#             combo[kmer]=combo[kmer]+freq
#         else:
#             combo[kmer]=freq
#     except ValueError:
#         continue
# firsts.close()
# for line in lasts:
#     kmer,freq=line.split('\t')
#     try:
#         freq=int(freq)
#         if kmer in combo.keys():
#             combo[kmer]=combo[kmer]+freq
#         else:
#             combo[kmer]=freq
#     except ValueError:
#         continue
# lasts.close()
# outfile = open('frequencies.txt','w')
# for key,value in combo.iteritems():
#     outfile.write('%s\t%s\n'%(key,value))
# outfile.close()
# infile = open('queries.txt','r')
# o1 = open('q1.txt','w')
# o2 = open('q2.txt','w')
# o3 = open('q3.txt','w')
# ln =0
# temp =dict()
# temp['124312']=0.843828
# temp['14124']=0.3782
# temp['13211']=0.10120
# for line in infile:
#     if ln < 3000:
#         o1.write(line)
#     if ln >=3000 and ln < 6000:
#         o2.write(line)
#     if ln >=6000:
#         o3.write(line)
#     ln+=1
# o1.close()
# o2.close()
# o3.close()
# infile.close()
#
# def cutdict(dic):
#     for g in (sorted(dic,key=dic.get,reverse=False)[50:]):
#         dic.pop(g)
# cutdict(temp)
# def equalize(bstring1,bstring2):
#     for i in range(0,abs((bstring1.bit_length())-bstring2.bit_length())):
#         bstring1=bstring1*2
#     return bstring1
# def equalize(bstring1,bstring2):
#     if bstring1.bit_length() <bstring2.bit_length():
#         for i in range(0,abs((bstring1.bit_length())-bstring2.bit_length())):
#             bstring1=bstring1*2
#         return bstring1
#     elif bstring2.bit_length() <bstring1.bit_length():
#         for i in range(0,abs((bstring1.bit_length())-bstring2.bit_length())):
#             bstring2=bstring2*2
#         return bstring2
# a=225
# b=132
# print bin(a)
# print bin(b)
#
# if a.bit_length <b.bit_length:
#             a=equalize(a,b)
# elif a.bit_length >b.bit_length:
#             kcomp=equalize(qcode[0],kmer[0])
#
#
#
# if a.bit_length()<b.bit_length():
#     a=equalize(a,b)
# elif a.bit_length()>b.bit_length():
#     b=equalize(a,b)
#
# print bin(a)
# print bin(b)
# # print bin(a^b)
#print bin(equalize(a,b))
#print bin(equalize(a,b))
# # print bin(b)
# # print bin(equalize(a,b)^b)
#
#
# print bin(a&b)
infile=open('frequencies.txt','r') #input the frequencies
symb2freq=dict()
for line in infile:
    kmer,freq=line.split('\t')
    try:
        freq=int(freq)
    except ValueError:
        continue
    symb2freq[kmer]=freq
alf=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
for letter in alf:
    for letter2 in alf:
        if letter+letter2 not in symb2freq:
            symb2freq[letter+letter2]=1
            print 'insert made'
outfile=open('masterfrequencies.txt','w')
for kmer,freq in symb2freq.iteritems():
    outfile.write(kmer+'\t%s\n'%freq)