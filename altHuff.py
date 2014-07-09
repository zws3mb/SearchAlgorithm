__author__ = 'seidz'
import Queue
import threading
import time
import numpy as numpy
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

def onecount(value): #count the number of 1's in a binary sequence.
    count=0;
    while(value!=0):
        value = value&(value-1);
        count=count+1;
    return count


vowels=['A','E','I','O','U','Y','AE','AI','AO','AU','AY','EA','EI','EO','EU','EY','IO','IU','IY','OA','OE','OI','OU','OY','UA','UE','UI','UY']
def syllable(name):         #crude syllabification
    first=False
    syll=1
    if name[0] in vowels:
        first=True
    temp=list()
    temp.append(1)
    for part in range(1,len(name)):
        if name[part] in vowels:
            if not name[part-1]+name[part] in vowels: #and name[part-1]+name[part] !='IA':
                if first:
                    syll+=1
                    #syll=1
                else:
                    first=True
        #else:
            #syll+=1
        temp.append(syll)
    return temp


def freqToTree(filename):       #utility to convert lines of element\tfreq into a dictionary
    #Building the Huffman tree
    infile=open(filename,'r') #input the frequencies
    symb2freq=defaultdict()
    for line in infile:
        kmer,freq=line.split('\t')
        try:
            freq=int(freq)
        except ValueError:
            continue
        symb2freq[kmer]=freq                #store the frequency with the kmer
    return symb2freq


def translate(huffman):                 #converts array to dict form
    huff = huffman#encode(freqToTree('sampleoutput2lasts.txt'))                #build the tree
    code=dict()
    for p in huff:
        #print "%s\t%s\t%s" % (p[0], symb2freq[p[0]], p[1])
        code[p[0]]=p[1]                     #put kmers and codes into dictionary form
    return code
def equalize(bstring1,bstring2):
    if bstring1.bit_length() <bstring2.bit_length():
        for i in range(0,abs((bstring1.bit_length())-bstring2.bit_length())):
            bstring1=bstring1*2
        return bstring1
    elif bstring2.bit_length() <bstring1.bit_length():
        for i in range(0,abs((bstring1.bit_length())-bstring2.bit_length())):
            bstring2=bstring2*2
        return bstring2
    else:
        print 'ERROR IN EQUALIZATION %s %s' %(bstring1,bstring2)
        return bstring1
def clean(input):                       #string cleaning
    return input.strip().replace('-','').replace('\'','').replace('\\','').replace('/','').replace('"','').replace(',','').replace(' ','').replace('.','').replace('?','').replace('1','').replace('2','').replace('3','').replace('4','').replace('5','').replace('6','').replace('7','').replace('8','').replace('9','').replace('0','').replace('(',"").replace(')','')


def treeFactory(code):
    try:
        treesrc=open('tree.txt','r')
        for line in treesrc:
            kmer,freq=line.split('\t')
            code[kmer]=freq
        treesrc.close()
    except:
        code = translate(encode(freqToTree('frequencies.txt')))         #Use the utilities to actually build the tree
        '''OUTPUT THE TREE FOR FUTURE USE'''
        outtree=open('tree.txt','w')
        for id, bstring in code.iteritems():
            outtree.write('%s\t%s\n'%(id,bstring))
        #print firsts
        outtree.close()
code = dict()
treeFactory(code)
#Read in the 'watchlist' or db to search against
KLENGTH=2                                   #length of the kmer
firsts=dict()
lasts=dict()
lookupf=dict()
lookupl=dict()
def loadNames(lookupf,lookupl):
    names=open('index.txt','r')
    for eachLine in names:
        index, first, last = eachLine.split('|')
        first=clean(first)
        last=clean(last).replace('\n','')
        if first not in lookupf:
            lookupf[first]=list()
        lookupf[first].append(index)
        if last not in lookupl:
            lookupl[last]=list()
        lookupl[last].append(index)
        if (first.find(' ')!=-1):
            temp0=clean(first.split(' ')[0]) #Just takes first names, for now. We'll see if there's a difference in the distributions.
            firsts[temp0]=temp0
        else:
            temp0=clean(first)
            firsts[temp0]=temp0
        templ=clean(last)
        lasts[templ]=templ
    names.close()

loadNames(lookupf, lookupl)
# lookupf=[None]*(900000)
# lookupl=[None]*(900000)
#Encode each of the names according to the huffman tree
encodefirsts=dict()

for index,name in firsts.iteritems():
    #val=int(index)
    #lookupf[val]=name
    if len(name)>0:
        encodefirsts[index]=list()
    for i in range(len(name)):
        if name[i:i+KLENGTH] in code:
            encodefirsts[index].append(int(code[name[i:i+KLENGTH]],2))
            #temp=temp+int(code[name[i:i+2]]).bit_length()
        else:
            print 'Error encoding first name'+name
            encodefirsts.pop(index)
            break
encodelasts=dict()
for index,lname in lasts.iteritems():
    #val=int(index)

    #lookupl[val]=lname
    if len(lname)>0:
        encodelasts[index]=list()
    for i in range(len(lname)):
        if lname[i:i+KLENGTH] in code:
            encodelasts[index].append(int(code[lname[i:i+KLENGTH]],2))
            #temp=temp+int(code[name[i:i+2]]).bit_length()
        else:
            print 'Error encoding last name'+lname
            encodelasts.pop(index)
            break


#print "QUERY:%s "%(qscore),
# print qcode
#1|ANN P|DANIELS
# outfile = open('results.txt','w')

#outfile.write('%s\n'%(temp))
temp=0
# for key in sorted(code, key=code.get, reverse=True):
#     print key+' '+code[key]
#print qcode

#print encode['ANN'][0].bit_length()
CUTOFF=1
sylldict=dict()
for index, kmer in encodefirsts.iteritems():
    sylldict[index]=syllable(index)
for index,kmer in encodelasts.iteritems():
    if index not in sylldict:
        sylldict[index]=syllable(index)
# '''EXPERIMENTAL NUMPY'''
# encodefirsts=numpy.array([(k,)+numpy.asarray(v) for k,v in encodefirsts.iteritems()])
# encodelasts=numpy.array([(k,)+numpy.asarray(v) for k,v in encodelasts.iteritems()])


outfile=open('r3.txt','a')
def reconcile(iq,firstres,lastres):
    output=0
    attempt=0
    for g,v in firstres.iteritems():#(sorted(firstres,key=firstres.get,reverse=False)[:50]):
        # if attempt<1000:
        #     for h in sorted(lastres,key=lastres.get,reverse=False):
        #         if output <100:
        #             if g==h:
        #                 outfile.write('%s|%s|%s\n'%(i,g,firstres[g]+lastres[g]))
        #                 print '%s|%s|%s'%(i,g,firstres[g]+lastres[g])
        #             output+=1
        #         else:
        #             break
        #     output=0
        #     attempt+=1
        # else:
        #     break
        #print '%s|%s|%s\n'%(iq,g,v)
        for id in lookupf[g]:
            for h,y in lastres.iteritems():
                if id in lookupl[h]:
                    outfile.write('%s|%s|%s\n'%(iq,id,v+y))
                    #print '%s|%s|%s\n'%(iq,id,v+y)
    # for h,v in lastres.iteritems():#(sorted(lastres,key=lastres.get,reverse=False)[:50]):
    #     outfile.write('%s|%s|%s\n'%(iq,h,v))
    #     print '%s|%s|%s\n'%(iq,h,v)
    print iq
def algorithm (query, huffdict, lookup):
    encode=huffdict
    Q=query
    qcode=list()
    for i in range(len(Q)):
        if Q[i:i+KLENGTH] in code:
            qcode.append(int(code[Q[i:i+KLENGTH]],2))

    val=0
    min=9999999999
    firstval=999999
    results=dict()
    #print encode['ANN'][0].bit_length()
    CUTOFF=1
    for chars,kmer in encode.iteritems(): #For every name
        if len(kmer)>=len(qcode):
            g=0 #query index
            numpass=0
            prop=0.0
            i =0
            syll=sylldict[chars]
            #for i in range(len(kmer)): #coded name kmer index
            while i <=len(kmer):
                if i==len(kmer):
                    if(val/prop<min and val<=prop*CUTOFF):
                        min=val/prop
                    break
                klen=kmer[i].bit_length()
                if klen==0:
                    klen=15
                if(i!=0 and i==len(qcode)+numpass):#i%(len(qcode)+1)==0):
                    # if(syll==klen/2):
                    #      syll=0
                    # syll+=1
                    numpass+=1
                    g=0
                    i=numpass
                    if(val/prop<min and val<=prop*CUTOFF):
                        min=val/prop
                        #print 'MIN SET:%s'%(min)
                    val=0
                    prop=0.0
                #if(kmer[i].bit_length()>prop):
                qcomp=qcode[g]
                kcomp=kmer[i]
                if qcode[g].bit_length() <kmer[i].bit_length():
                    qcomp=equalize(qcode[g],kmer[i])
                elif qcode[g].bit_length() >kmer[i].bit_length():
                    kcomp=equalize(qcode[g],kmer[i])

                prop=prop+max(klen,qcode[g].bit_length())+syll[i]*val
                val=val+onecount((qcomp^ kcomp))+abs(klen-qcode[g].bit_length())+(syll[i]*val)
                if len(chars)==-1:
                    print chars+':'+chars[i:i+2]+'::'+Q+':'+Q[g:g+2]+' I:%s G:%s'%(i,g)
                    print bin(qcomp)
                    #print bin(equalize(qcode[g],kmer[i]))
                    print bin(kcomp)
                    print '---------------'
                    print bin((qcomp^kcomp))
                    print 'SYL %s DIFF %s'%(syll[i],abs(kmer[i].bit_length()-qcode[g].bit_length()))
                    # # print onecount((qcode[g] ^ kmer[i]))
                    print 'VAL %s'%(val)
                    print 'PROP %s'%(prop)
                if i==0 and g==0:
                    firstval=val/prop
                g=g+1
                i=i+1
                if(numpass==0 and i==g==len(qcode)):
                    min=val/prop
                #                firstval=val/prop
                # print 'FIRST MIN SET:%s'%(min)
                if(len(kmer)==len(qcode)):
                    #if(val<min):
                    min=val/prop

                #print val,
        elif len(kmer)<=2:
            qcomp=qcode[0]
            kcomp=kmer[0]
            if qcode[0].bit_length() <kmer[0].bit_length():
                qcomp=equalize(qcode[0],kmer[0])
            elif qcode[0].bit_length() >kmer[0].bit_length():
                kcomp=equalize(qcode[0],kmer[0])
            if kmer[0]==qcode[0]:
                min=(onecount((kcomp^qcomp))+abs(kmer[0].bit_length()-qcode[0].bit_length()))/(max(kmer[0].bit_length(),qcode[0].bit_length()))
                firstval=min
            numpass=1
            #else:
            #min=len(qcode)*namesize
        elif len(qcode)>len(kmer):
            g=0 #query index
            numpass=0
            prop=0.0
            i=0
            syll=syllable(Q)
            #for i in range(len(qcode)):
            while i <=len(qcode):
                if i==len(qcode):
                    if(val/prop<min and val<=prop*CUTOFF):
                        min=val/prop
                    break
                if(i!=0 and i==len(kmer)+numpass):#i%(len(kmer))==0):
                    # if(syll==len(kmer)/2):
                    #      syll=0
                    # syll+=1

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
                qcomp=qcode[i]
                kcomp=kmer[g]
                if qcode[i].bit_length() <kmer[g].bit_length():
                    qcomp=equalize(qcode[i],kmer[g])
                elif qcode[i].bit_length() >kmer[g].bit_length():
                    kcomp=equalize(qcode[i],kmer[g])

                #if(kmer[i].bit_length()>prop):
                prop=prop+max(qcode[i].bit_length(),kmer[g].bit_length())+syll[i]*val
                val=val+onecount((kcomp ^ qcomp))+abs(kmer[g].bit_length()-qcode[i].bit_length())+(syll[i]*val)
                if len(kmer)==-1:
                    print Q+':'+Q[i:i+2]+'::'+chars+':'+chars[g:g+2]+' I:%s G:%s'%(i,g)
                    print bin(kcomp)
                    # print bin(equalize(kmer[g],qcode[i]))
                    print bin(qcomp)
                    print '---------------'
                    print bin((kcomp ^ qcomp))
                    print 'SYL %s DIFF %s'%(syll[i],abs(kmer[g].bit_length()-qcode[i].bit_length()))
                    # # print onecount((qcode[g] ^ kmer[i]))
                    print 'VAL %s'%(val)
                    print 'PROP %s'%(prop)
                if i==0 and g==0:
                    firstval=val/prop
                if(numpass==0 and i==g==len(qcode)):
                    min=val/prop
                #    firstval=val/prop
                g=g+1
                i=i+1
        fweight=max(abs((len(kmer)-len(qcode)))/10.0 ,1)
        results[chars]=(min+fweight*firstval)/2.0#+(.1*onecount(kmer[0]^qcode[0])*abs(max(len(kmer),len(qcode))-6))#/float(min(len(kmer),len(qcode)))#+(abs(len(kmer)-len(qcode))+1)#/(len(kmer)*len(qcode))
        #if min <.95:
        #print '%s %s'%(onecount(kmer[0]^qcode[0]),abs(max(len(kmer),len(qcode))-6))
        val=0
        min=99999999
        firstval=9999999
    return results
print 'PREPROCESSING FINISHED'
firstanswers=dict()
lastanswers=dict()
def cutdict(dic):
    for g in (sorted(dic,key=dic.get,reverse=False)[500:]):
        dic.pop(g)

'''PULL QUEUE CONTENTS FROM FILE'''
dataFile = open('q3.txt','r')
firstqs=dict()
lastqs=dict()
for eachLine in dataFile:
# setup a temporay variable 1|ANN P|DANIELS
    index, first, last = eachLine.split('|')
    index=int(index)
    if (first.find(' ')!=-1):
        firstqs[index]=clean(first.split(' ')[0]) #Just takes first names, for now. We'll see if there's a difference in the distributions.
        #firsts[index]=temp.strip()
    else:
        firstqs[index]=clean(first)
    lastqs[index]=clean(last)
dataFile.close()
'''END QUEUE PULL'''
imin=numpy.min(firstqs.keys())
imax=numpy.max(firstqs.keys())
# Fill the queue
for i in range(imin,imax):
    dataf=firstqs[i]
    datag=lastqs[i]
    ind =i
    stime=time.time()
    print "PROCESSING %s %s" % (dataf,datag)
    if dataf in firstanswers:
        firstresults=firstanswers[dataf]
    else:
        firstresults=algorithm(dataf,encodefirsts,lookupf)
        cutdict(firstresults)
        firstanswers[dataf]=firstresults
    print 'FIRST-TIME: %s'%(time.time()-stime)
    mtime=time.time()
    if datag in lastanswers:
        lastresults=lastanswers[datag]
    else:
        lastresults=algorithm(datag,encodelasts,lookupl)
        cutdict(lastresults)
        lastanswers[datag]=lastresults
    print 'LAST-TIME: %s'%(time.time()-mtime)
    rtime=time.time()
    print "RECONCILING RESULTS..."
    reconcile(ind,firstresults,lastresults)
    etime=time.time()
    print 'R-TIME: %s'%(etime-rtime)
    print "%s %s FINISHED: %s" % (dataf,datag,etime-stime)

# Notify threads it's time to exit




# for i in range(0,len(firstqs)):
#     print 'Running query:'+firstqs['%s'%i]+' '+lastqs['%s'%i]
#     starttime=time.time()
#     firstresults=algorithm(firstqs['%s'%i],encodefirsts,lookupf)
#     endtime=time.time()
#     print 'Firstname completed'
#     print endtime-starttime
#     starttime=time.time()
#     lastresults=algorithm(lastqs['%s'%i],encodelasts,lookupl)
#     endttime=time.time()
#     print 'Lastname completed'
#     print endttime-starttime
#     thread.start_new_thread(reconcile,(firstresults,lastresults))
outfile.close()
# print len(firstresults)
# print len(lastresults)
#firstresults=sorted(firstresults,key=firstresults.get,reverse=True)
#lastresults=sorted(lastresults,key=lastresults.get,reverse=False)

# for w in sorted(results, key=results.get, reverse=False):
     #outfile.write(w+'\t%s\n'%(results[w]))
#     print(w+'\t%s'%(results[w]))