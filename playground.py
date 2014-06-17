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
    return count
# print '\t%b'%(0b101)
# print '\t%b'%(0b01011)
# print bin((int('01011',2)))
# print bin(int(0b101)^int(0b01011))
print int('110001000',2).bit_length()-int('11',2).bit_length()