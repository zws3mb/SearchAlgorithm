__author__ = 'seidz'
__author__ = 'seidz'
#!/usr/bin/env python

import sys

current_kmer = None
kmer = None
infile=open('sampleassociation0last.txt','r')
infile=sorted(infile)
print infile[1]
#outfile = open('sampleassociation2last.txt','w')
# input comes from STDIN
for line in infile:
    # remove leading and trailing whitespace
    line = line.strip()
    #print line
    # parse the input we got from mapper.py
    kmer, v1 = line.split('\t')

    # convert count (currently a string) to int
    try:
        v1 = int(v1)
    except ValueError:
        # count was not a number, so silently
        # ignore/discard this line
        print "error"
        continue

    # this IF-switch only works because Hadoop sorts map output
    # by key (here: word) before it is passed to the reducer
    if current_kmer == kmer:
        current_v1 += v1
    else:
        if current_kmer:
            # write result to STDOUT
            print('%s\t%s\n' % (current_kmer, current_v1))
        current_kmer = kmer
        current_v1 = v1



# do not forget to output the last word if needed!
if current_kmer == kmer:
    print('%s\t%s\n' % (current_kmer, current_v1))
