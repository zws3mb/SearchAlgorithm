__author__ = 'seidz'
from heapq import heappush, heappop, heapify
from collections import defaultdict

def encode(symb2freq):
    """Huffman encode the given dict mapping symbols to weights"""
    heap = [[wt, [sym, ""]] for sym, wt in symb2freq.items()]
    heapify(heap)
    while len(heap) > 1:
        lo = heappop(heap)
        hi = heappop(heap)
        for pair in lo[1:]:
            pair[1] = '0' + pair[1]
        for pair in hi[1:]:
            pair[1] = '1' + pair[1]
        heappush(heap, [lo[0] + hi[0]] + lo[1:] + hi[1:])
    return sorted(heappop(heap)[1:], key=lambda p: (len(p[-1]), p))
def onecount(value):
    count=0;
    while(value!=0):
        value = value&(value-1);
        count=count+1;
    return count
KLENGTH=2
#infile=open('custom2.txt','r')
infile=open('sampleoutput2.txt','r')
symb2freq=defaultdict()
for line in infile:
    kmer,freq=line.split('\t')
    try:
        freq=int(freq)
    except ValueError:
        continue
    symb2freq[kmer]=freq
huff = encode(symb2freq)
#print "Symbol\tWeight\tHuffman Code"
code=dict()
for p in huff:
    #print "%s\t%s\t%s" % (p[0], symb2freq[p[0]], p[1])
    code[p[0]]=p[1]

names=open('sample.txt','r')
firsts=dict()
namesize=0
for eachLine in names:
    namesize=namesize+1
    index, first, last = eachLine.split('|')
    if (first.find(' ')!=-1):
        firsts[index]=first.split(' ')[0].strip().replace('-','').replace('\'','') #Just takes first names, for now. We'll see if there's a difference in the distributions.
        #firsts[index]=temp.strip()
    else:
        firsts[index]=first.strip().replace('-','').replace('\'','')
temp=0
#print firsts
Q="ZACHARY"
qcode=list()
for i in range(len(Q)):
    if Q[i:i+KLENGTH] in code:
        #query=query+" "+code["ANN"[i:i+2]]
        qcode.append(int(code[Q[i:i+KLENGTH]],2))
        #qlen=qlen+int(code["ANN"[i:i+2]]).bit_length()
qscore=0
for item in qcode:
    qscore=qscore+(item&item)
#print "QUERY:%s "%(qscore),
print qcode
#1|ANN P|DANIELS
outfile = open('results.txt','w')
encode=dict()
for name in firsts.values():
    encode[name]=list()
    for i in range(len(name)):
        if name[i:i+KLENGTH] in code:
            encode[name].append(int(code[name[i:i+KLENGTH]],2))
            #temp=temp+int(code[name[i:i+2]]).bit_length()
    #outfile.write('%s\n'%(temp))
    temp=0
# for key in sorted(code, key=code.get, reverse=True):
#     print key+' '+code[key]
#print qcode
val=0
min=0
results=dict()
#print encode['ANN'][0].bit_length()
CUTOFF=.5
for chars,kmer in encode.iteritems(): #For every name
    if len(kmer)>=len(qcode):
        g=0 #query index
        numpass=0
        prop=0.0
        i =0
        syll=1
        #for i in range(len(kmer)): #coded name kmer index
        while i <len(kmer):
            klen=kmer[i].bit_length()
            if klen==0:
                klen=7
            if(i!=0 and i==len(qcode)+numpass):#i%(len(qcode)+1)==0):
                if(syll==klen/2):
                     syll=0
                syll+=1
                numpass+=1
                g=0
                i=numpass
                if(val/prop<min and val<=prop*CUTOFF):
                    min=val/prop
                    # print 'MIN SET:%s'%(min)
                val=0
                prop=0.0
            #if(kmer[i].bit_length()>prop):
            prop=prop+klen+(i)*(g+1)#+abs(kmer[i].bit_length()-qcode[g].bit_length())#qcode[g].bit_length()
            val=val+onecount((qcode[g] ^ kmer[i]))+abs(klen-qcode[g].bit_length())+(i)*(g)#*(syll)#+(1+syll)
            # print chars+':'+chars[i:i+2]+' I:%s G:%s'%(i,g)
            # print bin(qcode[g])
            # print bin(kmer[i])
            # print bin((qcode[g] ^ kmer[i]))
            # print onecount((qcode[g] ^ kmer[i]))
            # print 'PROP %s'%(prop)
            # print 'VAL %s'%(val)
            g=g+1
            i=i+1
            if(numpass==0 and i==g==len(qcode)):
                    min=val/prop
                    print 'FIRST MIN SET:%s'%(min)
            if(len(kmer)==len(qcode)):
                #if(val<min):
                min=val/prop
            #print val,
    elif len(kmer)<=2:
        if kmer[0]==qcode[0]:
            min=onecount(kmer[len(kmer)-1]^qcode[len(qcode)-1])/(kmer[len(kmer)-1].bit_length())
        numpass=1
        #else:
            #min=len(qcode)*namesize
    elif len(qcode)>len(kmer):
        g=0 #query index
        numpass=0
        prop=0.0
        i=0
        syll=1
        #for i in range(len(qcode)):
        while i < len(qcode):
            if(i!=0 and i==len(kmer)+numpass):#i%(len(kmer))==0):
                if(syll==len(kmer)/2):
                     syll=0
                syll+=1

                if(i==(len(kmer)+numpass)):
                    min=val/prop
                numpass=numpass+1
                g=0
                i=numpass
                #if(val<min):
                if(val/prop<min and val<=prop*CUTOFF):
                    min=val/prop
                val=0
                prop=0.0
            #if(kmer[i].bit_length()>prop):
            prop=prop+qcode[i].bit_length()+(i)*(g)#+abs(kmer[g].bit_length()-qcode[i].bit_length()))#kmer[g].bit_length()

            val=val+onecount((kmer[g] ^ qcode[i]))+abs(kmer[g].bit_length()-qcode[i].bit_length())+(i)*(g)#*(syll)+abs(len(kmer)-len(qcode))#(g/(1+syll))#+(1+syll)# count of differences between codes, penalized for position in each record

            g=g+1
            #val=val|(qcode[i] ^ kmer[i])/(g+1)
            #min=val
            #print val,
            i=i+1
    results[chars]=min#/float(min(len(kmer),len(qcode)))#+(abs(len(kmer)-len(qcode))+1)#/(len(kmer)*len(qcode))
    val=0
    min=99999999999999999
for w in sorted(results, key=results.get, reverse=True):
    #outfile.write(w+'\t%s\n'%(results[w]))
    print(w+'\t%s'%(results[w]))