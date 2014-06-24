import Queue
import threading
import time
import thread
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
def syllable(name):
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
#print syllable('ANNIE')
def freqToTree(filename):
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
def translate(huffman):
    huff = huffman#encode(freqToTree('sampleoutput2lasts.txt'))                #build the tree
    code=dict()
    for p in huff:
        #print "%s\t%s\t%s" % (p[0], symb2freq[p[0]], p[1])
        code[p[0]]=p[1]                     #put kmers and codes into dictionary form
    return code
def clean(input):
    return input.strip().replace('-','').replace('\'','').replace('\\','').replace('/','').replace('"','').replace(',','').replace(' ','').replace('.','').replace('?','').replace('1','').replace('2','').replace('3','').replace('4','').replace('5','').replace('6','').replace('7','').replace('8','').replace('9','').replace('0','')



code = translate(encode(freqToTree('frequencies.txt')))
#Read in the 'watchlist' or db to search against
KLENGTH=2                                   #length of the kmer
names=open('sample.txt','r')
firsts=dict()
lasts=dict()
for eachLine in names:
    index, first, last = eachLine.split('|')
    if (first.find(' ')!=-1):
        firsts[index]=clean(first.split(' ')[0]) #Just takes first names, for now. We'll see if there's a difference in the distributions.
        #firsts[index]=temp.strip()
    else:
        firsts[index]=clean(first)
    lasts[index]=clean(last)
names.close()
lookupf=[None]*(900000)
lookupl=[None]*(900000)
#Encode each of the names according to the huffman tree
encodefirsts=dict()
for index,name in firsts.iteritems():
    val=int(index)
    lookupf[val]=name
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
    val=int(index)

    lookupl[val]=lname
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

temp=0
#print firsts

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
    sylldict[lookupf[int(index)]]=syllable(kmer)
for index,kmer in encodelasts.iteritems():
    if lookupl[int(index)] not in sylldict:
        sylldict[lookupl[int(index)]]=syllable(kmer)
outfile=open('results.txt','a')
def reconcile(firstres,lastres):
    output=0
    attempt=0
    for g in sorted(firstres,key=firstres.get,reverse=False)[:50]:
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
        print '%s|%s|%s\n'%(i,g,firstres[g])
        outfile.write('%s|%s|%s\n'%(i,g,firstres[g]))
    for h in sorted(lastres,key=lastres.get,reverse=False)[:50]:
        outfile.write('%s|%s|%s\n'%(i,h,lastres[h]))
        print '%s|%s|%s\n'%(i,h,lastres[h])
def algorithm (query, huffdict, lookup):
    encode=huffdict
    Q=query
    qcode=list()
    for i in range(len(Q)):
        if Q[i:i+KLENGTH] in code:
            qcode.append(int(code[Q[i:i+KLENGTH]],2))

    val=0
    min=0
    results=dict()
    for chars,kmer in encode.iteritems(): #For every name
        val=0
        min=0
        if len(kmer)>=len(qcode):
            g=0 #query index
            numpass=0
            prop=0.0
            i =0
            syll=sylldict[lookup[int(chars)]]
            while i <len(kmer):
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
                prop=prop+klen+val*syll[i]#+abs(kmer[i].bit_length()-qcode[g].bit_length())#qcode[g].bit_length()
                val=val+onecount((qcode[g] ^ kmer[i]))+abs(klen-qcode[g].bit_length())+(val)*(syll[i])#+(1+syll)
                # print chars+':'+chars[i:i+2]+' I:%s G:%s'%(i,g)
                # print bin(qcode[g])
                # print bin(kmer[i])
                # print bin((qcode[g] ^ kmer[i]))
                # print onecount((qcode[g] ^ kmer[i]))
                # print 'SYLL %s'%(syll[i])
                # print 'PROP %s'%(prop)
                # print 'VAL %s = %s + %s'%(val,abs(klen-qcode[g].bit_length()),val*syll[i])
                g=g+1
                i=i+1
                if(numpass==0 and i==g==len(qcode)):
                    min=val/prop
                    #print 'FIRST MIN SET:%s'%(min)
                if(len(kmer)==len(qcode)):
                    #if(val<min):
                    min=val/prop
                    #print val,
        elif len(kmer)<=2:
            if kmer[0]==qcode[0]:
                min=onecount(kmer[len(kmer)-1]^qcode[len(qcode)-1])+abs(kmer[len(kmer)-1].bit_length()-qcode[len(qcode)-1].bit_length())
            else:
                min=2.0 #arbitrary high proportion of min differences
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
            while i < len(qcode):
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
                #if(kmer[i].bit_length()>prop):
                prop=prop+qcode[i].bit_length()+val*syll[i]#+abs(kmer[g].bit_length()-qcode[i].bit_length()))#kmer[g].bit_length()

                val=val+onecount((kmer[g] ^ qcode[i]))+abs(kmer[g].bit_length()-qcode[i].bit_length())+val*(syll[i])#+abs(len(kmer)-len(qcode))#(g/(1+syll))#+(1+syll)# count of differences between codes, penalized for position in each record

                g=g+1
                #val=val|(qcode[i] ^ kmer[i])/(g+1)
                #min=val
                #print val,
                i=i+1
        results[chars]=min#/float(min(len(kmer),len(qcode)))#+(abs(len(kmer)-len(qcode))+1)#/(len(kmer)*len(qcode))
    val=0
    min=99999999999999999
    return results
print 'PREPROCESSING FINISHED'
exitFlag = 0
class myThread (threading.Thread):
    def __init__(self, threadID, name, q):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.q = q
    def run(self):
        #print "Starting " + self.name
        process_data(self.name, self.q)
        #print "Exiting " + self.name

def process_data(threadName, q):
    while not exitFlag:
        queueLock.acquire()
        if not workQueue.empty():
            dataf,datag = q.get()
            queueLock.release()
            firstresults=algorithm(dataf,encodefirsts,lookupf)
            lastresults=algorithm(datag,encodelasts,lookupl)
            reconcile(firstresults,lastresults)
            print "processing %s %s" % (dataf,datag)
        else:
            queueLock.release()
        time.sleep(1)

threadList = ["Thread-1", "Thread-2", "Thread-3","Thread-4"]
nameList = ["One", "Two", "Three", "Four", "Five"]
queueLock = threading.Lock()
workQueue = Queue.Queue()

threads = []
threadID = 1

# Create new threads
for tName in threadList:
    thread = myThread(threadID, tName, workQueue)
    thread.start()
    threads.append(thread)
    threadID += 1

'''PULL QUEUE CONTENTS FROM FILE'''
dataFile = open('diagnostic.txt','r')
firstqs=dict()
lastqs=dict()
for eachLine in dataFile:
# setup a temporay variable 1|ANN P|DANIELS
    index, first, last = eachLine.split('|')
    if (first.find(' ')!=-1):
        firstqs[index]=clean(first.split(' ')[0]) #Just takes first names, for now. We'll see if there's a difference in the distributions.
        #firsts[index]=temp.strip()
    else:
        firstqs[index]=clean(first)
    lastqs[index]=clean(last)
dataFile.close()
'''END QUEUE PULL'''


# Fill the queue
queueLock.acquire()
for i in range(1,3):
    workQueue.put((firstqs['%s'%i],lastqs['%s'%i]))
queueLock.release()

# Wait for queue to empty
while not workQueue.empty():
    pass

# Notify threads it's time to exit
exitFlag = 1

# Wait for all threads to complete
for t in threads:
    t.join()


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